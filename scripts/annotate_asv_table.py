#!/usr/bin/env python3
"""
annotate_asv_table.py — Merge taxonomy sources into a single annotated ASV table.

Output column order:
  ASV
  | primary taxonomy (kraken2 or dada2):
      domain  phylum  class  order  family  genus  species
  | vsearch taxonomy (higher-resolution, separate column set):
      vsearch_domain  ...  vsearch_species
      vsearch_identity  vsearch_top_hit  vsearch_n_hits
  | per-sample read counts: s1  s2  s3  ...
  | total

This produces asv_table_annotated.tsv alongside the original asv_table.tsv

"""

import argparse
import os
import sys

import pandas as pd


TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--asv_table",       required=True,
                   help="asv_table.tsv from DADA2 (ASV × samples, counts)")
    p.add_argument("--primary_taxonomy", required=True,
                   help="Taxonomy TSV from kraken2 or dada2 rule")
    p.add_argument("--vsearch_taxonomy", required=True,
                   help="Taxonomy TSV from vsearch_species rule")
    p.add_argument("--output",           required=True,
                   help="Output annotated TSV path")
    p.add_argument("--method",           default="kraken2",
                   help="'kraken2' or 'dada2' — used for primary column prefix")
    return p.parse_args()


def log(msg):
    print(f"[annotate_asv_table] {msg}", flush=True)


def normalise_asv_col(df, name="input"):
    """Rename whatever the ASV ID column is called to 'ASV' (string)."""
    candidates = [c for c in df.columns
                  if c.upper() in ("ASV", "ASVID", "ASV_ID", "#OTU ID", "SEQUENCE")]
    if candidates:
        df = df.rename(columns={candidates[0]: "ASV"})
    else:
        # Fall back to first column
        df = df.rename(columns={df.columns[0]: "ASV"})
    df["ASV"] = df["ASV"].astype(str)
    return df


def normalise_tax_columns(df, prefix=""):
    """
    Rename taxonomy rank columns to a consistent set.
    Handles GTDB-style (d__Bacteria) and plain (Bacteria) values.
    Renames columns to prefix_rank or rank depending on prefix.
    """
    col_map = {}
    for col in df.columns:
        col_lower = col.lower()
        for rank in TAX_RANKS:
            if rank in col_lower and col_lower not in ("asv", "asvid"):
                target = f"{prefix}{rank}" if prefix else rank
                col_map[col] = target
                break
    return df.rename(columns=col_map)


def main():
    args = parse_args()

    for path, name in [
        (args.asv_table,        "asv_table"),
        (args.primary_taxonomy, "primary_taxonomy"),
        (args.vsearch_taxonomy, "vsearch_taxonomy"),
    ]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: {name} not found: {path}")

    # -------------------------------------------------------------------------
    # Load ASV count table
    # -------------------------------------------------------------------------
    log("Loading ASV count table...")
    counts = pd.read_csv(args.asv_table, sep="\t", dtype=str)
    counts = normalise_asv_col(counts, "asv_table")
    sample_cols = [c for c in counts.columns if c != "ASV"]
    # Coerce counts to numeric
    for col in sample_cols:
        counts[col] = pd.to_numeric(counts[col], errors="coerce").fillna(0).astype(int)
    counts["total"] = counts[sample_cols].sum(axis=1)
    log(f"  {len(counts):,} ASVs × {len(sample_cols)} samples")

    # -------------------------------------------------------------------------
    # Load primary taxonomy (kraken2 or dada2)
    # -------------------------------------------------------------------------
    log(f"Loading primary taxonomy ({args.method})...")
    primary = pd.read_csv(args.primary_taxonomy, sep="\t", dtype=str)
    primary = normalise_asv_col(primary, "primary_taxonomy")
    primary = normalise_tax_columns(primary, prefix="")  # plain rank names
    # Keep only ASV + rank columns
    keep_primary = ["ASV"] + [r for r in TAX_RANKS if r in primary.columns]
    primary = primary[keep_primary]
    log(f"  {len(primary):,} rows | ranks found: "
        f"{[r for r in TAX_RANKS if r in primary.columns]}")

    # -------------------------------------------------------------------------
    # Load vsearch taxonomy
    # -------------------------------------------------------------------------
    log("Loading vsearch taxonomy...")
    vsearch = pd.read_csv(args.vsearch_taxonomy, sep="\t", dtype=str)
    vsearch = normalise_asv_col(vsearch, "vsearch_taxonomy")
    # vsearch columns already prefixed vsearch_* — keep all of them
    vsearch_cols = ["ASV"] + [c for c in vsearch.columns if c.startswith("vsearch_")]
    vsearch = vsearch[vsearch_cols]
    log(f"  {len(vsearch):,} rows | columns: {vsearch_cols[1:]}")

    # -------------------------------------------------------------------------
    # Merge: ASV → primary taxonomy → vsearch taxonomy → sample counts → total
    # -------------------------------------------------------------------------
    log("Merging tables...")
    merged = counts[["ASV"]].copy()

    # Primary taxonomy columns first
    merged = merged.merge(primary, on="ASV", how="left")

    # vsearch taxonomy columns second
    merged = merged.merge(vsearch, on="ASV", how="left")

    # Sample count columns + total last
    merged = merged.merge(counts[["ASV"] + sample_cols + ["total"]], on="ASV", how="left")

    # Sort by total reads descending
    merged["total"] = pd.to_numeric(merged["total"], errors="coerce").fillna(0)
    merged = merged.sort_values("total", ascending=False).reset_index(drop=True)

    # -------------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)

    log(f"Annotated table: {len(merged):,} ASVs × {len(merged.columns)} columns")
    log(f"Column order: ASV | primary taxonomy | vsearch taxonomy | sample counts | total")

    # Species call rate summary
    # FIX: a duplicate column name (for example two "Species" columns) makes
    # merged[col] return a DataFrame, not a Series. n was then a Series and
    # f"{n:,}" raised:
    #   TypeError: unsupported format string passed to Series.__format__
    # That killed the rule after the merge had already succeeded. Take the
    # first matching column so the count is always a scalar.
    for col, label in [("species", "Primary species"), ("vsearch_species", "Vsearch species")]:
        if col in merged.columns:
            sub = merged[col]
            if isinstance(sub, pd.DataFrame):
                sub = sub.iloc[:, 0]
            n = int(sub.notna().sum())
            log(f"  {label}: {n:,} / {len(merged):,} ASVs "
                f"({100*n/max(len(merged),1):.1f}%)")

    log(f"Output: {args.output}")


if __name__ == "__main__":
    main()
