#!/usr/bin/env python3
"""
extract_discarded_reads.py — Extract reads discarded by DADA2 filterAndTrim.

Strategy (no R changes required):
  DADA2's filterAndTrim reads from the cutadapt-trimmed FASTQs and writes
  only the passing reads. By collecting the read IDs that ARE in the passing
  filtered FASTQs and subtracting from the trimmed FASTQs, we recover the
  discarded reads.

  Discarded = trimmed_reads - filtered_reads

  Output: per-sample FASTA of discarded reads, then merged across samples
  into a single FASTA for Kraken2 classification.
"""

import argparse
import gzip
import os
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--trimmed_dir",  required=True,
                   help="Directory containing cutadapt-trimmed FASTQs")
    p.add_argument("--filtered_dir", required=True,
                   help="Directory containing DADA2 filterAndTrim output FASTQs")
    p.add_argument("--output_fasta", required=True,
                   help="Output merged FASTA of all discarded reads")
    p.add_argument("--read1_suffix", default="_R1.fastq.gz")
    p.add_argument("--samples",      nargs="+", required=True,
                   help="Sample names to process")
    return p.parse_args()


def read_ids_in_fastq(fastq_path):
    """Return set of read IDs present in a FASTQ (gzipped or plain)."""
    ids = set()
    if not os.path.exists(fastq_path):
        return ids
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            fh.readline()   # seq
            fh.readline()   # +
            fh.readline()   # qual
            # Read ID is everything after '@' up to the first space
            rid = header.strip().lstrip("@").split()[0]
            ids.add(rid)
    return ids


def extract_discarded(trimmed_fastq, passing_ids, out_fh, sample):
    """Write reads from trimmed_fastq NOT in passing_ids to out_fh as FASTA."""
    n_total    = 0
    n_written  = 0
    opener = gzip.open if str(trimmed_fastq).endswith(".gz") else open
    with opener(trimmed_fastq, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().strip()
            fh.readline()   # +
            fh.readline()   # qual
            n_total += 1
            rid = header.strip().lstrip("@").split()[0]
            if rid not in passing_ids:
                # Write as FASTA with sample prefix
                out_fh.write(f">{sample}|{rid}\n{seq}\n")
                n_written += 1
    return n_total, n_written


def main():
    args = parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output_fasta)), exist_ok=True)

    total_trimmed   = 0
    total_passing   = 0
    total_discarded = 0

    with open(args.output_fasta, "w") as out_fh:
        for sample in args.samples:
            trimmed_r1  = os.path.join(args.trimmed_dir,  f"{sample}{args.read1_suffix}")
            filtered_r1 = os.path.join(args.filtered_dir, f"{sample}{args.read1_suffix}")

            if not os.path.exists(trimmed_r1):
                print(f"[extract_discarded] WARNING: trimmed R1 not found for {sample}, skipping",
                      flush=True)
                continue

            passing_ids = read_ids_in_fastq(filtered_r1)
            n_trim, n_disc = extract_discarded(trimmed_r1, passing_ids, out_fh, sample)

            n_pass = n_trim - n_disc
            total_trimmed   += n_trim
            total_passing   += n_pass
            total_discarded += n_disc
            pct = 100.0 * n_disc / n_trim if n_trim > 0 else 0
            print(f"[extract_discarded] {sample}: {n_trim:,} trimmed | "
                  f"{n_pass:,} passed DADA2 | {n_disc:,} discarded ({pct:.1f}%)",
                  flush=True)

    print(f"\n[extract_discarded] Total: {total_trimmed:,} trimmed reads | "
          f"{total_discarded:,} discarded ({100*total_discarded/max(total_trimmed,1):.1f}%)",
          flush=True)
    print(f"[extract_discarded] Discarded FASTA: {args.output_fasta}", flush=True)


if __name__ == "__main__":
    main()
