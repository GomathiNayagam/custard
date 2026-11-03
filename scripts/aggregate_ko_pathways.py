#!/usr/bin/env python3
"""
aggregate_ko_pathways.py — Roll per-ASV KO scores up into named KEGG pathways.

asv_predicted_kos.tsv (from infer_kos.py) is a per-ASV matrix of ~thousands
of bare KO identifiers (K00001, K00002, ...) with no biological names
attached — not directly useful to a reader. This script:

  1. Computes per-sample KO abundance: (samples x ASVs) . (ASVs x KOs),
     the same weighting infer_kos/ko_pca already use.
  2. Maps each KO to its KEGG pathway(s) via a static reference table
     (data/ko_pathway_map.tsv, pulled once from the KEGG REST API).
  3. Sums KO abundance per pathway. A KO belonging to multiple pathways
     contributes its full abundance to each — the standard convention
     (same approach PICRUSt2/HUMAnN use), since a gene's measured copy
     number isn't divided across the pathways it participates in.

Output: pathway_abundance.tsv, samples x named pathways. This is the
directly-useful table — "Glycolysis / Gluconeogenesis: 0.42" instead of
forty bare KO columns nobody can read without a KEGG lookup.
"""

import argparse
import sys

import pandas as pd


def log(msg):
    print(f"[aggregate_ko_pathways] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--asv_table", required=True)
    p.add_argument("--ko_table", required=True)
    p.add_argument("--ko_pathway_map", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    asv_table = pd.read_csv(args.asv_table, sep="\t", index_col=0)
    ko_table = pd.read_csv(args.ko_table, sep="\t", index_col=0)

    common_asvs = asv_table.index.intersection(ko_table.index)
    if len(common_asvs) == 0:
        sys.exit(
            "ERROR: no ASVs shared between asv_table and ko_table. "
            "Check that ASV IDs match in both files."
        )
    asv_table = asv_table.loc[common_asvs]
    ko_table = ko_table.loc[common_asvs]

    log(f"Computing per-sample KO abundance ({len(asv_table.columns)} samples "
        f"x {len(ko_table.columns)} KOs)...")
    ko_abundance = asv_table.T.dot(ko_table)   # samples x KOs

    log("Loading KO -> pathway map...")
    ko_map = pd.read_csv(args.ko_pathway_map, sep="\t")
    # One KO may map to several pathways; keep every (ko, pathway) pair.
    ko_map = ko_map[ko_map["ko"].isin(ko_abundance.columns)]
    if ko_map.empty:
        sys.exit(
            "ERROR: none of the predicted KOs matched the ko_pathway_map "
            "reference. Check --ko_pathway_map points at data/ko_pathway_map.tsv."
        )
    log(f"  {ko_map['ko'].nunique():,} / {ko_abundance.shape[1]:,} predicted "
        f"KOs map to a named pathway ({ko_map['pathway_id'].nunique():,} pathways)")

    # Column label combines id + name so two differently-numbered KEGG
    # pathways can never collide under one label.
    ko_map["label"] = ko_map["pathway_id"] + ": " + ko_map["pathway_name"]

    pathway_abundance = pd.DataFrame(
        0.0, index=ko_abundance.index, columns=sorted(ko_map["label"].unique())
    )
    for label, grp in ko_map.groupby("label"):
        kos_in_pathway = grp["ko"].tolist()
        pathway_abundance[label] = ko_abundance[kos_in_pathway].sum(axis=1)

    # Drop pathways with zero abundance across every sample (uninformative).
    pathway_abundance = pathway_abundance.loc[:, (pathway_abundance > 0).any(axis=0)]
    log(f"  {pathway_abundance.shape[1]:,} pathways with nonzero abundance "
        f"in at least one sample")

    pathway_abundance.index.name = "Sample"
    pathway_abundance.to_csv(args.output, sep="\t")
    log(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
