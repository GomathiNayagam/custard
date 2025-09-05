#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import subprocess
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv_fasta", required=True)
    parser.add_argument("--ko_fasta", required=True)
    parser.add_argument("--ko_npz", required=True)
    parser.add_argument("--ko_genomes", required=True)
    parser.add_argument("--ko_list", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--identity_threshold", type=float, default=0.4)
    args = parser.parse_args()

    # Load KO sparse matrix
    matrix = load_npz(args.ko_npz)
    genomes = pd.read_csv(args.ko_genomes, header=None)[0].tolist()
    kos = pd.read_csv(args.ko_list, header=None)[0].tolist()

    genome_to_kos = {}
    for i, genome in enumerate(genomes):
        genome_to_kos[genome] = matrix.getrow(i).toarray().flatten()

    # Build genome header mapping from KO fasta
    ko_headers = {}
    for rec in SeqIO.parse(args.ko_fasta, "fasta"):
        gid = rec.id.split("|")[0]  # remove suffixes, take GCF_/GCA_
        ko_headers[rec.id] = gid

    # Perform vsearch search
    vsearch_out = "tmp_vsearch_hits.tsv"
    subprocess.run([
        "vsearch", "--usearch_global", args.asv_fasta,
        "--db", args.ko_fasta,
        "--id", str(args.identity_threshold),
        "--blast6out", vsearch_out,
        "--threads", str(args.threads),
        "--maxhits", "1000"
    ], check=True)

    # Parse vsearch output
    # blast6 format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    hits = defaultdict(list)
    with open(vsearch_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            qseqid, sseqid, pident, align_len = parts[0], parts[1], float(parts[2])/100, int(parts[3])
            weight = pident * align_len  # weighted score
            hits[qseqid].append((ko_headers[sseqid], weight))

    # Compute fractional KO assignments
    result = {}
    for asv, genome_hits in hits.items():
        # sum weights per genome
        genome_weights = defaultdict(float)
        for gid, w in genome_hits:
            genome_weights[gid] += w
        total_weight = sum(genome_weights.values())

        ko_scores = defaultdict(float)
        for gid, gw in genome_weights.items():
            if gid not in genome_to_kos:
                continue
            presence = genome_to_kos[gid]
            for idx, val in enumerate(presence):
                if val == 1:
                    ko_scores[kos[idx]] += gw / total_weight

        result[asv] = ko_scores

    # Write output TSV
    all_kos = kos
    with open(args.output, "w") as out:
        header = ["ASV"] + all_kos
        out.write("\t".join(header) + "\n")
        for asv in result:
            row = [asv] + [f"{result[asv].get(k,0):.3f}" for k in all_kos]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
