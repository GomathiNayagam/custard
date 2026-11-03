#!/usr/bin/env python3
"""
pathway_diff_abundance.py — Which pathways actually differ between groups.

ko_pca.py shows whether sample groups separate in functional space but
never says what's driving that separation. This script tests every
pathway in pathway_abundance.tsv for a difference between the metadata
groups and reports a ranked table — the thing a reader actually wants:
"which functions are enriched where."

Test: Mann-Whitney U for exactly two groups, Kruskal-Wallis for more than
two (both non-parametric — no distributional assumption on abundance
scores, and robust with the small sample sizes typical of amplicon
studies). P-values are FDR-corrected (Benjamini-Hochberg) across all
pathways tested, since testing hundreds of pathways at once needs
multiple-testing correction or most "hits" are noise.
"""

import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests


def log(msg):
    print(f"[pathway_diff_abundance] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pathway_abundance", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--alpha", type=float, default=0.05,
                    help="FDR significance threshold for the 'significant' flag (default 0.05)")
    return p.parse_args()


def main():
    args = parse_args()

    abund = pd.read_csv(args.pathway_abundance, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata, sep="\t", index_col=0)

    common = abund.index.intersection(meta.index)
    if len(common) < 2:
        sys.exit(
            f"ERROR: only {len(common)} sample(s) overlap between "
            "pathway_abundance and metadata — need at least 2."
        )
    abund = abund.loc[common]
    group_col = "Group" if "Group" in meta.columns else meta.columns[0]
    groups = meta.loc[common, group_col]

    group_names = sorted(groups.unique())
    n_groups = len(group_names)
    if n_groups < 2:
        sys.exit(
            f"ERROR: only one group ({group_names}) present — need at least "
            "2 groups to test for a difference."
        )
    log(f"{len(common)} samples, {n_groups} groups: {group_names}")

    test_name = "Mann-Whitney U" if n_groups == 2 else "Kruskal-Wallis"
    log(f"Testing {abund.shape[1]:,} pathways ({test_name})...")

    rows = []
    for pathway in abund.columns:
        values = abund[pathway]
        by_group = [values[groups == g].values for g in group_names]

        # A pathway that is all-zero, or constant, in every group carries no
        # information and a rank-based test on it is undefined - skip rather
        # than report a meaningless p-value.
        if all(np.all(g == g[0]) for g in by_group) and len(set(np.concatenate(by_group))) <= 1:
            continue

        try:
            if n_groups == 2:
                stat, pval = mannwhitneyu(by_group[0], by_group[1], alternative="two-sided")
            else:
                stat, pval = kruskal(*by_group)
        except ValueError:
            # All values identical across every group - kruskal/mannwhitneyu
            # raise rather than return a p-value; genuinely uninformative.
            continue

        row = {"pathway": pathway, "test": test_name, "statistic": stat, "p_value": pval}
        for g, vals in zip(group_names, by_group):
            row[f"mean[{g}]"] = vals.mean()

        if n_groups == 2:
            m0, m1 = row[f"mean[{group_names[0]}]"], row[f"mean[{group_names[1]}]"]
            pseudo = max(m0, m1, 1e-6) * 1e-3  # avoid log2(0) without distorting real signal
            row["log2fc"] = np.log2((m1 + pseudo) / (m0 + pseudo))

        rows.append(row)

    if not rows:
        sys.exit("ERROR: no pathway had enough variation across groups to test.")

    result = pd.DataFrame(rows)
    result["p_adj"] = multipletests(result["p_value"], method="fdr_bh")[1]
    result["significant"] = result["p_adj"] < args.alpha
    result = result.sort_values("p_adj")

    result.to_csv(args.output, sep="\t", index=False)
    n_sig = result["significant"].sum()
    log(f"Wrote {args.output}  ({n_sig}/{len(result)} pathways significant at FDR < {args.alpha})")


if __name__ == "__main__":
    main()
