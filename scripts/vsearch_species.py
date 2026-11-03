#!/usr/bin/env python3
"""
vsearch_species.py — Species-level taxonomy via lowest-common-ancestor of the
best vsearch hits.

Runs vsearch against a reference FASTA with headers in CUSTARD format:
  >GCF_016756535.1;taxonomy=d__Bacteria;p__Pseudomonadota;...;s__Mesorhizobium sp016756535

For each ASV the search is exhaustive, every reference tied at the best score is
collected, and the reported lineage is the lowest common ancestor of that tied
set. A sole best hit yields its full lineage; N tied references yield only what
they agree on. Reported columns:
  - vsearch_species  : species, or empty when the best hits disagree there
  - vsearch_genus    : genus, likewise
  - vsearch_identity : identity of the best (tied) hits (%)
  - vsearch_top_hit  : best single hit genome ID + identity
  - vsearch_n_hits   : number of hits that passed filters
  - vsearch_n_tied   : how many references tied for best — >1 means the call
                       was capped at their LCA
  - vsearch_lca_rank : deepest rank actually resolved

Note for anyone comparing against older output: this used to run vsearch with
its default --maxaccepts 1, which returns the first reference above --id and
stops. Every ASV therefore had exactly one hit, the advertised top-N weighted
consensus never ran, and the reported species was whichever equally-good
reference vsearch reached first. On the ZymoBIOMICS standard that silently
assigned E. coli reads to s__Enterobacter_flexneri_A, one of seven conflicting
labels carried by that byte-identical amplicon. Species calls are fewer now and
correspondingly more trustworthy.

This is run independently and produces a separate TSV merged into the
annotated ASV table by the annotate_asv_table rule.

Usage (called from Snakemake rule):
  python scripts/vsearch_species.py
      --asv_fasta       results/asv_sequences.fasta
      --ref_fasta       data/db/dada2_formatted_rpoc_amplicon.fasta
      --output          results/asv_vsearch_taxonomy.tsv
      --threads         20
      --identity        0.97      # species-level threshold (97% = classic species)
      --min_coverage    0.9       # high coverage required for reliable species call
      --top_n           5
"""

import argparse
import os
import subprocess
import re
import sys
import tempfile
from collections import defaultdict

import pandas as pd
from Bio import SeqIO


RANK_ORDER = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX = {
    "domain":  "d__",
    "phylum":  "p__",
    "class":   "c__",
    "order":   "o__",
    "family":  "f__",
    "genus":   "g__",
    "species": "s__",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--asv_fasta",    required=True)
    p.add_argument("--ref_fasta",    required=True,
                   help="DADA2 reference FASTA with ;taxonomy= headers")
    p.add_argument("--output",       required=True)
    p.add_argument("--threads",      type=int,   default=1)
    p.add_argument("--identity",     type=float, default=0.97,
                   help="Min identity for species call (default 0.97 = 97%%)")
    p.add_argument("--min_coverage", type=float, default=0.9,
                   help="Min query coverage fraction (default 0.9)")
    p.add_argument("--top_n",        type=int,   default=5,
                   help="Top N hits to use for consensus (default 5)")
    p.add_argument("--max_hits",     type=int,   default=100)
    return p.parse_args()


def log(msg):
    print(f"[vsearch_species] {msg}", flush=True)


def parse_custard_header(header):
    """
    Parse a CUSTARD FASTA header into {rank: value} dict.

    Accepts rec.description (full header after '>') so species names
    containing spaces are preserved intact. For example:
      's__Mesorhizobium sp016756535'
    is a single semicolon-delimited field — splitting on ';' keeps the
    space within the value rather than truncating it.

    Header format (rec.description):
      GCF_016756535.1;taxonomy=d__Bacteria;p__Pseudomonadota;...;s__Mesorhizobium sp016756535

    Returns (genome_id, tax_dict).
    """
    header = header.lstrip(">").strip()
    parts  = header.split(";")
    genome = parts[0].split()[0]   # genome ID = first whitespace-delimited token

    fields = [p.strip() for p in parts[1:]]
    fields = [p[len("taxonomy="):] if p.startswith("taxonomy=") else p for p in fields]

    tax_dict = {}
    for part in fields:
        for rank, prefix in RANK_PREFIX.items():
            if part.startswith(prefix):
                # Capture everything after prefix - spaces within the value are fine
                val = part[len(prefix):].strip()
                if val:
                    tax_dict[rank] = val
                break   # one rank per field

    # FIX: accept references that carry NO d__/p__/.../s__ rank prefixes.
    # rpocDB writes the prefixes. The GTDB SSU release does not. Without this
    # fallback every rank column came back empty for the 16S arm (0 of 556 ASVs)
    # while rpoC parsed normally (92%). One arm silently lost its vsearch
    # taxonomy, which is worse than a hard failure because nothing reports it.
    # A bare header holds exactly the 7 ranks in order, for example:
    #   ID;Bacteria;Pseudomonadota;...;Escherichia;Escherichia coli
    if not tax_dict and len(fields) == len(RANK_ORDER):
        for rank, val in zip(RANK_ORDER, fields):
            val = val.strip()
            if val:
                tax_dict[rank] = val

    return genome, tax_dict


def index_ref_headers(ref_fasta):
    """
    Return {seq_id: tax_dict} for every record in the reference FASTA.

    BioPython's rec.id stops at the first whitespace, which truncates species
    names that contain a space (e.g. 's__Mesorhizobium sp016756535' becomes
    's__Mesorhizobium'). We use rec.description (the full header after '>')
    for taxonomy parsing, but key the dict by rec.id so it matches the
    sequence IDs returned by vsearch blast6 output.
    """
    log("Indexing reference taxonomy headers...")
    index = {}
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        # rec.description = full header line after '>'
        # rec.id          = everything up to the first space (used as blast6 sseqid)
        _, tax = parse_custard_header(rec.description)
        index[rec.id] = tax
    log(f"  {len(index):,} reference sequences indexed")
    return index


def get_query_lengths(fasta_path):
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}


def run_vsearch(asv_fasta, ref_fasta, identity, max_hits, threads, out_path):
    log(f"Running vsearch (identity >= {identity:.0%}, max_hits={max_hits})...")
    cmd = [
        "vsearch", "--usearch_global", asv_fasta,
        "--db",        ref_fasta,
        "--id",        str(identity),
        "--blast6out", out_path,
        "--threads",   str(threads),
        "--maxhits",   str(max_hits),
        "--query_cov", "0.0",
        # vsearch defaults to --maxaccepts 1 / --maxrejects 32, i.e. it stops at
        # the FIRST reference above --id and returns that one. --maxhits alone
        # does not change this: it caps output, not the search. Without these
        # two flags every ASV came back with exactly one hit, the top-N
        # weighted consensus below never actually ran, and the single hit
        # reported was whichever equally-good reference vsearch happened to
        # reach first. 0/0 makes the search exhaustive, so the top-N vote and
        # the co-optimal tie check see the real candidate pool.
        "--maxaccepts", "0",
        "--maxrejects", "0",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"vsearch failed: {r.returncode}")
    log("vsearch done.")


# Scores within this relative distance of the best are treated as tied. Exact
# ties are the norm, not the exception: rpocDB holds ~1,700 amplicons that are
# byte-identical yet carry conflicting labels (one Collinsella amplicon carries
# 41 species labels), so identity and coverage — and therefore score — come out
# bit-for-bit equal.
TIE_MARGIN = 1e-6



# --- GTDB label normalisation -------------------------------------------------
# Two label artefacts otherwise force a genuine single-species match down to genus:
#
#  1. Infraspecific suffixes. GTDB splits a species into ANI subclusters labelled
#     s__X_y_A, _B, ... Those are the SAME species, so a tie between
#     s__Streptococcus_parasanguinis_M and _V must resolve to
#     s__Streptococcus_parasanguinis, not to g__Streptococcus.
#     Genus suffixes (g__Blautia_A, and the genus token inside the species
#     binomial) are real GTDB genera and are left untouched.
#
#  2. Placeholder species. s__Genus_sp012345678 is the ABSENCE of a name, not a
#     competing name. When a tie mixes named binomials with placeholders, the
#     placeholders are dropped; they are used only if nothing named remains.
_SP_SUFFIX = re.compile(r"_[A-Z]{1,2}$")
_PLACEHOLDER = re.compile(r"_sp\d{4,}$")


def _strip_species_suffix(species):
    """s__Streptococcus_parasanguinis_M -> s__Streptococcus_parasanguinis"""
    return _SP_SUFFIX.sub("", species) if species else species


def _is_placeholder(species):
    return bool(species) and bool(_PLACEHOLDER.search(species))


def normalise_species_values(vals):
    """Collapse a set of tied species labels to what they actually disagree on."""
    vals = {v for v in vals if v}
    if not vals:
        return vals
    named = {v for v in vals if not _is_placeholder(v)}
    if named:                      # placeholders never veto a real name
        vals = named
    return {_strip_species_suffix(v) for v in vals}


def consensus_lineage(cooptimal_tax):
    """
    Resolve a lineage from domain down to species as the lowest common ancestor
    of the co-optimal hits — every reference tied at the best score.

    Those references explain the ASV equally well, so a rank they disagree on
    is not resolvable, and a rank they agree on is. Ranks are walked top-down
    and the first disagreement (or first unannotated rank) terminates the
    lineage; everything below is left as None.

    When one reference is the sole best hit this returns its full lineage, which
    is the common case. When N references tie it returns what they share — for
    the Zymo E. coli amplicon, which is byte-identical to references labelled
    s__Escherichia_coli, s__Enterobacter_flexneri_A, s__Escherichia_boydii and
    four more, that is f__Enterobacteriaceae.

    This deliberately replaces a weighted majority vote over the top-N hits.
    That vote could not work: vsearch was being run with its default
    --maxaccepts 1, so exactly one hit was ever returned and the "consensus"
    was that single arbitrary hit. Re-running it exhaustively exposed the other
    half of the problem — a vote lets several weaker hits overturn a unique
    perfect match, which is worse than not voting at all.

    Returns (ranks_dict, lca_rank); lca_rank names the deepest resolved rank,
    or None if nothing resolved.
    """
    ranks    = {}
    lca_rank = None
    stopped  = False

    for rank in RANK_ORDER:
        if stopped:
            ranks[rank] = None
            continue

        vals = {tax.get(rank) for _, tax in cooptimal_tax}
        if rank == "species":
            # Suffix and placeholder artefacts are not real disagreements.
            norm = normalise_species_values(vals)
            if len(norm) == 1 and None not in vals or (len(norm) == 1 and vals - {None, ""}):
                vals = norm
        if len(vals) != 1 or None in vals or "" in vals:
            ranks[rank] = None
            stopped = True
            continue

        ranks[rank] = vals.pop()
        lca_rank = rank

    return ranks, lca_rank


def parse_blast6_and_score(blast6_path, ref_index, query_lengths, min_coverage, top_n):
    """
    Parse blast6, apply coverage filter, collect hits per ASV.

    Returns (top_hits, cooptimal_hits):
      top_hits       {asv: [(score, tax, pident, hit_id), ...]}  best top_n
      cooptimal_hits {asv: [(score, tax, pident, hit_id), ...]}  every hit tied
                     at the best score — not capped at top_n, because the size
                     of the tied pool is exactly what the caller needs to know.
    """
    log(f"Parsing hits (min_coverage >= {min_coverage:.0%})...")
    all_hits   = defaultdict(list)
    n_total    = 0
    n_filtered = 0

    with open(blast6_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            qid, sid = parts[0], parts[1]
            pident   = float(parts[2]) / 100.0
            aln_len  = int(parts[3])
            n_total += 1

            qlen     = query_lengths.get(qid, 0)
            coverage = aln_len / qlen if qlen > 0 else 0.0
            if coverage < min_coverage:
                n_filtered += 1
                continue

            tax   = ref_index.get(sid, {})
            score = pident * coverage
            all_hits[qid].append((score, tax, pident * 100, sid))

    log(f"  {n_total:,} hits | {n_total - n_filtered:,} passed coverage filter")
    log(f"  ASVs with >= 1 hit: {len(all_hits):,}")

    top_hits       = {}
    cooptimal_hits = {}
    n_tied         = 0
    for asv, hits in all_hits.items():
        ordered = sorted(hits, key=lambda x: x[0], reverse=True)
        top_hits[asv] = ordered[:top_n]

        best   = ordered[0][0]
        co_opt = [h for h in ordered if h[0] >= best * (1.0 - TIE_MARGIN)]
        cooptimal_hits[asv] = co_opt
        if len(co_opt) > 1:
            n_tied += 1

    log(f"  ASVs whose best hit is tied by >= 1 other reference: {n_tied:,} "
        f"({100*n_tied/max(len(all_hits),1):.1f}%)")
    return top_hits, cooptimal_hits


def build_taxonomy_table(top_hits, cooptimal_hits, all_asvs, identity_thresh):
    """
    Build output DataFrame. One row per ASV.
    Columns: ASV | vsearch_domain | ... | vsearch_species |
             vsearch_identity | vsearch_top_hit | vsearch_n_hits
    """
    rows      = []
    n_lca     = 0   # ASVs where a tie forced a shallower call
    for asv in all_asvs:
        hits = top_hits.get(asv)
        if not hits:
            row = {"ASV": asv}
            for rank in RANK_ORDER:
                row[f"vsearch_{rank}"] = None
            row["vsearch_identity"]  = None
            row["vsearch_top_hit"]   = None
            row["vsearch_n_hits"]    = 0
            row["vsearch_lca_rank"]  = None
            row["vsearch_n_tied"]    = 0
            rows.append(row)
            continue

        co_opt   = cooptimal_hits.get(asv, hits[:1])
        co_tax   = [(h[0], h[1]) for h in co_opt]
        # Identity of the co-optimal hits — they are tied, so this is just the
        # best identity. Reported as vsearch_identity and used for the species
        # threshold. Previously this averaged the top-N, which let weaker hits
        # drag a perfect match below the cutoff.
        best_id  = sum(h[2] for h in co_opt) / len(co_opt)
        top_hit  = f"{hits[0][3]} ({hits[0][2]:.1f}%)"

        ranks, lca_rank = consensus_lineage(co_tax)

        # For species: only report if identity >= threshold
        if best_id < identity_thresh * 100:
            ranks["species"] = None
            if lca_rank == "species":
                lca_rank = "genus" if ranks.get("genus") else lca_rank

        row = {"ASV": asv}
        for rank in RANK_ORDER:
            row[f"vsearch_{rank}"] = ranks.get(rank)

        if lca_rank != "species":
            n_lca += 1

        row["vsearch_identity"] = round(best_id, 2)
        row["vsearch_top_hit"]  = top_hit
        row["vsearch_n_hits"]   = len(hits)
        row["vsearch_lca_rank"] = lca_rank
        row["vsearch_n_tied"]   = len(co_opt)
        rows.append(row)

    df = pd.DataFrame(rows)

    n_with_species = df["vsearch_species"].notna().sum()
    n_total        = len(df)
    log(f"  Species-level calls: {n_with_species:,} / {n_total:,} ASVs "
        f"({100*n_with_species/max(n_total,1):.1f}%)")
    log(f"  Resolved below species (tie or missing rank): {n_lca:,} ASVs "
        f"— see vsearch_lca_rank")
    return df


def main():
    args = parse_args()

    for path, name in [(args.asv_fasta, "asv_fasta"), (args.ref_fasta, "ref_fasta")]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: {name} not found: {path}")

    if args.max_hits < args.top_n:
        args.max_hits = args.top_n * 20

    ref_index     = index_ref_headers(args.ref_fasta)
    query_lengths = get_query_lengths(args.asv_fasta)
    all_asvs      = list(query_lengths.keys())
    log(f"ASVs: {len(all_asvs):,} | Identity threshold: {args.identity:.0%} | Top-N: {args.top_n}")

    tmp_fd, tmp_path = tempfile.mkstemp(
        suffix=".tsv", prefix=f"vsearch_sp_{os.getpid()}_",
        dir=os.path.dirname(os.path.abspath(args.output)) or "."
    )
    os.close(tmp_fd)

    try:
        run_vsearch(args.asv_fasta, args.ref_fasta,
                    args.identity, args.max_hits, args.threads, tmp_path)
        top_hits, cooptimal_hits = parse_blast6_and_score(
            tmp_path, ref_index, query_lengths, args.min_coverage, args.top_n
        )
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    df = build_taxonomy_table(top_hits, cooptimal_hits, all_asvs, args.identity)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    log(f"Output: {args.output}")

    # Summary
    log("")
    log("=== Vsearch Species Summary ===")
    for rank in RANK_ORDER:
        col = f"vsearch_{rank}"
        if col in df.columns:
            n = df[col].notna().sum()
            log(f"  {rank:<10}: {n:,} / {len(df):,} ASVs ({100*n/max(len(df),1):.1f}%)")
    log("================================")


if __name__ == "__main__":
    main()
