#!/usr/bin/env python3
"""
infer_kos.py — Infer KEGG Orthology (KO) assignments for ASVs via vsearch.

Improvements over v1 (single-best-hit):
  1. Weighted top-N averaging — top N hits per ASV weighted by identity×coverage.
  2. Batch matrix indexing — collects all genome row indices first, then does a
     single matrix slice instead of one getrow() per genome per ASV. For large
     runs this replaces an O(ASVs × top_N) sparse access loop with a single
     scipy fancy-index operation, giving substantial speedup at scale.
  3. Header map handles both GCF_xxx|... and GCF_xxx;taxonomy=... formats.
  4. Auto-raises max_hits if below top_n_hits.
"""

import argparse
import os
import sys
import subprocess
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.sparse import load_npz


def parse_args():
    p = argparse.ArgumentParser(
        description="Infer KO abundances from ASV sequences via vsearch."
    )
    p.add_argument("--asv_fasta",           required=True)
    p.add_argument("--ko_fasta",            required=True)
    p.add_argument("--ko_npz",              required=True)
    p.add_argument("--ko_genomes",          required=True)
    p.add_argument("--ko_list",             required=True)
    p.add_argument("--output",              required=True)
    p.add_argument("--threads",             type=int,   default=1)
    p.add_argument("--identity_threshold",  type=float, default=0.4)
    p.add_argument("--min_coverage",        type=float, default=0.5)
    p.add_argument("--top_n_hits",          type=int,   default=10,
                   help="Top N hits to average KO profiles over per ASV (default 10). "
                        "Set 1 for single-best-hit.")
    p.add_argument("--max_hits",            type=int,   default=1000)
    p.add_argument("--write_zeros",         action="store_true")
    return p.parse_args()


def log(msg):
    print(f"[infer_kos] {msg}", flush=True)


def get_query_lengths(fasta_path):
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def load_genome_ko_matrix(ko_npz, ko_genomes_path, ko_list_path):
    log("Loading KO sparse matrix...")
    matrix  = load_npz(ko_npz)
    genomes = pd.read_csv(ko_genomes_path, header=None)[0].tolist()
    kos     = pd.read_csv(ko_list_path,    header=None)[0].tolist()
    genome_index = {g: i for i, g in enumerate(genomes)}
    log(f"  Matrix: {len(genomes):,} genomes × {len(kos):,} KOs  "
        f"(density: {matrix.nnz / (len(genomes) * len(kos)):.4%})")
    return matrix, genome_index, kos


def build_ko_header_map(ko_fasta_path):
    log("Indexing reference FASTA headers...")
    mapping = {}
    for rec in SeqIO.parse(ko_fasta_path, "fasta"):
        gid = rec.id.split("|")[0].split(";")[0]
        mapping[rec.id] = gid
    log(f"  {len(mapping):,} reference sequences indexed")
    return mapping


def run_vsearch(asv_fasta, ko_fasta, identity, max_hits, threads, out_path):
    log(f"Running vsearch (identity >= {identity:.0%}, max_hits={max_hits}, "
        f"threads={threads})...")
    cmd = [
        "vsearch", "--usearch_global", asv_fasta,
        "--db",       ko_fasta,
        "--id",       str(identity),
        "--blast6out", out_path,
        "--threads",  str(threads),
        "--maxhits",  str(max_hits),
        "--query_cov", "0.0",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"vsearch failed with exit code {result.returncode}")
    log("vsearch completed.")


def parse_blast6(blast6_path, ko_header_map, query_lengths, min_coverage):
    log(f"Parsing vsearch hits (min_coverage >= {min_coverage:.0%})...")
    hits          = defaultdict(list)
    total_hits    = 0
    filtered_hits = 0
    missing_refs  = set()

    with open(blast6_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            qid, sid = parts[0], parts[1]
            pident   = float(parts[2]) / 100.0
            aln_len  = int(parts[3])
            total_hits += 1

            qlen     = query_lengths.get(qid, 0)
            coverage = aln_len / qlen if qlen > 0 else 0.0
            if coverage < min_coverage:
                filtered_hits += 1
                continue

            gid = ko_header_map.get(sid)
            if gid is None:
                missing_refs.add(sid)
                continue

            hits[qid].append((gid, pident * coverage))

    if missing_refs:
        log(f"  WARNING: {len(missing_refs)} ref IDs not in header map "
            f"(first 5: {list(missing_refs)[:5]})")
    kept = total_hits - filtered_hits
    log(f"  Hits: {total_hits:,} total | {kept:,} passed coverage filter | "
        f"{filtered_hits:,} filtered")
    log(f"  ASVs with >= 1 hit: {len(hits):,}")
    return hits


def compute_ko_scores(hits, genome_index, matrix, kos, top_n_hits):
    """
    Weighted top-N KO averaging with batch matrix indexing.

    Batch approach:
      Instead of calling matrix.getrow(idx) inside the per-ASV loop
      (which creates a new sparse row object each time), we:
        1. Collect the union of all genome row indices needed across ALL ASVs.
        2. Slice them out of the matrix in one operation:
             submatrix = matrix[sorted_indices, :].toarray()
           This uses scipy's efficient CSR fancy indexing — one memory
           allocation for all needed rows.
        3. Build a lookup dict {genome_id: dense_row_vector} once.
        4. Per-ASV KO scoring is then pure numpy dot products over dense arrays.

    This replaces O(ASVs × top_N) sparse operations with 1 sparse slice +
    O(ASVs × top_N) dense array ops — dramatically faster for large runs.
    """
    log(f"Computing weighted top-{top_n_hits} KO assignments (batch matrix indexing)...")

    n_kos = len(kos)

    # --- Pass 1: compute top-N genome weights per ASV ---
    asv_genome_weights = {}   # {asv: {genome_id: weight}}
    needed_genomes     = set()

    for asv, genome_hits in hits.items():
        top_hits = sorted(genome_hits, key=lambda x: x[1], reverse=True)[:top_n_hits]

        gw = defaultdict(float)
        for gid, score in top_hits:
            if gid in genome_index:
                gw[gid] += score
        if not gw:
            continue

        total = sum(gw.values())
        asv_genome_weights[asv] = {g: w / total for g, w in gw.items()}
        needed_genomes.update(gw.keys())

    if not needed_genomes:
        log("  No valid genomes found — all ASVs have no KO hits")
        return {}

    # --- Pass 2: batch-fetch all needed genome rows in one matrix slice ---
    log(f"  Fetching {len(needed_genomes):,} genome rows from matrix (batch)...")
    sorted_genomes = sorted(needed_genomes)
    row_indices    = [genome_index[g] for g in sorted_genomes]

    # Single CSR fancy-index: shape (n_needed_genomes × n_kos), dense float32
    genome_matrix  = matrix[row_indices, :].toarray().astype(np.float32)
    genome_row_map = {g: genome_matrix[i] for i, g in enumerate(sorted_genomes)}

    # --- Pass 3: per-ASV weighted average using pre-fetched dense rows ---
    result  = {}
    skipped = 0

    for asv, gw_norm in asv_genome_weights.items():
        ko_scores = np.zeros(n_kos, dtype=np.float32)
        for gid, weight in gw_norm.items():
            ko_scores += weight * genome_row_map[gid]

        nonzero_idx = np.nonzero(ko_scores)[0]
        if len(nonzero_idx) == 0:
            skipped += 1
            continue
        result[asv] = {kos[i]: float(ko_scores[i]) for i in nonzero_idx}

    if skipped:
        log(f"  {skipped} ASVs produced all-zero KO vectors (no annotated KOs in top hits)")
    log(f"  KO scores computed for {len(result):,} ASVs")
    return result


def write_output(result, kos, all_asvs, output_path, write_zeros):
    log("Writing output...")
    active_kos = sorted({ko for scores in result.values() for ko in scores})
    log(f"  Active KOs: {len(active_kos):,} / {len(kos):,} total")

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with open(output_path, "w") as out:
        out.write("\t".join(["ASV"] + active_kos) + "\n")
        written = 0
        for asv in all_asvs:
            if asv not in result:
                if write_zeros:
                    out.write("\t".join([asv] + ["0.000"] * len(active_kos)) + "\n")
                    written += 1
                continue
            scores = result[asv]
            row = [asv] + [f"{scores.get(ko, 0.0):.3f}" for ko in active_kos]
            out.write("\t".join(row) + "\n")
            written += 1

    log(f"  Wrote {written:,} rows × {len(active_kos):,} columns → {output_path}")


def main():
    args = parse_args()

    for path, name in [
        (args.asv_fasta,  "asv_fasta"),
        (args.ko_fasta,   "ko_fasta"),
        (args.ko_npz,     "ko_npz"),
        (args.ko_genomes, "ko_genomes"),
        (args.ko_list,    "ko_list"),
    ]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: {name} not found: {path}")

    if args.max_hits < args.top_n_hits:
        new_max = args.top_n_hits * 10
        log(f"WARNING: max_hits ({args.max_hits}) < top_n_hits ({args.top_n_hits}). "
            f"Raising to {new_max}.")
        args.max_hits = new_max

    matrix, genome_index, kos = load_genome_ko_matrix(
        args.ko_npz, args.ko_genomes, args.ko_list
    )
    ko_header_map = build_ko_header_map(args.ko_fasta)
    query_lengths = get_query_lengths(args.asv_fasta)
    all_asvs      = list(query_lengths.keys())
    log(f"ASVs to process: {len(all_asvs):,} | Top-N: {args.top_n_hits}")

    tmp_fd, tmp_path = tempfile.mkstemp(
        suffix=".tsv",
        prefix=f"vsearch_hits_{os.getpid()}_",
        dir=os.path.dirname(os.path.abspath(args.output)) or "."
    )
    os.close(tmp_fd)

    try:
        run_vsearch(args.asv_fasta, args.ko_fasta,
                    args.identity_threshold, args.max_hits,
                    args.threads, tmp_path)
        hits   = parse_blast6(tmp_path, ko_header_map, query_lengths, args.min_coverage)
        result = compute_ko_scores(hits, genome_index, matrix, kos, args.top_n_hits)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
            log(f"Temp file removed: {tmp_path}")

    write_output(result, kos, all_asvs, args.output, args.write_zeros)

    n_with = len(result)
    n_tot  = len(all_asvs)
    mean_k = np.mean([len(v) for v in result.values()]) if result else 0
    log("")
    log("=== KO Inference Summary ===")
    log(f"  ASVs total          : {n_tot:,}")
    log(f"  ASVs with KO hits   : {n_with:,} ({n_with/max(n_tot,1):.1%})")
    log(f"  Mean KOs per ASV    : {mean_k:.1f}")
    log(f"  Identity threshold  : >= {args.identity_threshold:.0%}")
    log(f"  Coverage threshold  : >= {args.min_coverage:.0%}")
    log(f"  Top-N hits averaged : {args.top_n_hits}")
    log("============================")


if __name__ == "__main__":
    main()
