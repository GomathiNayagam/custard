#!/usr/bin/env python3
"""
ko_pca.py — PCA of weighted KO abundance per sample.

Key improvements over previous version:
  - No colored scatter dots. Groups shown via filled ellipses (95% CI) with
    group label centred inside each ellipse boundary.
  - Sample names placed inside the plot with font size that scales with N.
  - Adaptive figure size and layout so large experiments don't crowd labels.
  - Correct 95% confidence ellipse scaling via chi-squared quantile (was 2σ).
  - Single-sample and two-sample groups handled gracefully (no ellipse drawn,
    label placed at point).
  - Colour-blind-friendly palette (Wong 2011 / Okabe-Ito).
  - Tight layout with no legend clutter — group identity conveyed by ellipse
    colour + in-plot label only.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from matplotlib.colors import to_rgba
from sklearn.decomposition import PCA
from scipy.stats import chi2

# ---------------------------------------------------------------------------
# Snakemake I/O
# ---------------------------------------------------------------------------
asv_table_path = snakemake.input["asv_table"]
ko_table_path  = snakemake.input["ko_table"]
metadata_path  = snakemake.input["metadata"]
output_path    = snakemake.output[0]

os.makedirs(os.path.dirname(output_path), exist_ok=True)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
asv_table = pd.read_csv(asv_table_path, sep="\t", index_col=0)
ko_table  = pd.read_csv(ko_table_path,  sep="\t", index_col=0)
meta      = pd.read_csv(metadata_path,  sep="\t", index_col=0)

# Align ASVs
common_asvs = asv_table.index.intersection(ko_table.index)
if len(common_asvs) == 0:
    raise ValueError(
        "No ASVs shared between asv_table and ko_table. "
        "Check that ASV IDs match (e.g. ASV001 in both files)."
    )
asv_table = asv_table.loc[common_asvs]
ko_table  = ko_table.loc[common_asvs]

# Weighted KO abundance per sample: (samples × ASVs) · (ASVs × KOs)
ko_abundance = asv_table.T.dot(ko_table)

# Align samples
common_samples = ko_abundance.index.intersection(meta.index)
if len(common_samples) < 2:
    raise ValueError(
        f"Only {len(common_samples)} sample(s) overlap between KO table and "
        "metadata — need at least 2 for PCA."
    )
ko_abundance = ko_abundance.loc[common_samples]
meta         = meta.loc[common_samples]

# Drop all-zero KO columns (uninformative for PCA)
ko_abundance = ko_abundance.loc[:, (ko_abundance > 0).any(axis=0)]

n_samples = len(common_samples)

# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------
n_components = min(2, ko_abundance.shape[1], n_samples)
if n_components < 2:
    raise ValueError(
        "Not enough variance dimensions for 2-component PCA. "
        "Check that ko_abundance has at least 2 non-zero KO columns."
    )

pca   = PCA(n_components=2)
X_pca = pca.fit_transform(ko_abundance.values)

coords = pd.DataFrame(
    X_pca, index=ko_abundance.index, columns=["PC1", "PC2"]
)
# Attach group column — default to "Group" or first metadata column
group_col = "Group" if "Group" in meta.columns else meta.columns[0]
coords = coords.join(meta[[group_col]])
coords.rename(columns={group_col: "Group"}, inplace=True)

groups = coords["Group"].unique()
n_groups = len(groups)

# ---------------------------------------------------------------------------
# Colour palette — Okabe-Ito / Wong (colour-blind safe, up to 8 groups;
# cycles with slight lightness offset for larger sets)
# ---------------------------------------------------------------------------
OKABE_ITO = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
]

def make_palette(n):
    base = OKABE_ITO * ((n // len(OKABE_ITO)) + 1)
    return base[:n]

palette = make_palette(n_groups)
group_colors = dict(zip(groups, palette))

# ---------------------------------------------------------------------------
# Adaptive figure geometry
# ---------------------------------------------------------------------------
# Width and height grow with sample count but are capped to avoid
# extremely large files. Font size shrinks as labels crowd the space.
fig_w      = max(7, min(14, 6 + n_samples * 0.25))
fig_h      = max(6, min(11, 5 + n_samples * 0.20))
# Label font: starts at 9pt for few samples, shrinks to 5pt for many
label_fs   = max(5.0, min(9.0, 9.0 - (n_samples - 5) * 0.15))
# Group centroid label font
group_fs   = max(7.0, min(13.0, label_fs * 1.4))

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# ---------------------------------------------------------------------------
# 95% confidence ellipse helper
# ---------------------------------------------------------------------------
CHI2_95 = np.sqrt(chi2.ppf(0.95, df=2))   # ≈ 2.448; correct 95% CI scaling

def confidence_ellipse_95(x, y, ax, color, alpha_fill=0.13, alpha_edge=0.85, lw=2.0):
    """
    Draw a 95% confidence ellipse for (x, y) data.
    Returns the centroid (mean_x, mean_y).
    Requires >= 3 points; returns centroid without drawing for n < 3.
    """
    mean_x, mean_y = np.mean(x), np.mean(y)
    if len(x) < 3:
        return mean_x, mean_y

    cov = np.cov(x, y)
    # Guard against degenerate covariance (all points collinear)
    if np.linalg.matrix_rank(cov) < 2:
        return mean_x, mean_y

    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]

    angle  = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width  = 2 * CHI2_95 * np.sqrt(vals[0])
    height = 2 * CHI2_95 * np.sqrt(vals[1])

    fill_color = to_rgba(color, alpha=alpha_fill)
    edge_color = to_rgba(color, alpha=alpha_edge)

    # Filled ellipse
    ell_fill = Ellipse(
        xy=(mean_x, mean_y), width=width, height=height, angle=angle,
        facecolor=fill_color, edgecolor="none", zorder=2
    )
    # Ellipse border
    ell_edge = Ellipse(
        xy=(mean_x, mean_y), width=width, height=height, angle=angle,
        facecolor="none", edgecolor=edge_color, linewidth=lw,
        linestyle="--", zorder=3
    )
    ax.add_patch(ell_fill)
    ax.add_patch(ell_edge)
    return mean_x, mean_y

# ---------------------------------------------------------------------------
# Draw ellipses + group labels + sample name labels
# ---------------------------------------------------------------------------
label_padding = dict(
    boxstyle="round,pad=0.15", facecolor="white",
    edgecolor="none", alpha=0.65
)

for g, grp_data in coords.groupby("Group"):
    color = group_colors[g]
    x     = grp_data["PC1"].values
    y     = grp_data["PC2"].values

    cx, cy = confidence_ellipse_95(x, y, ax, color)

    # Group centroid label — drawn on top of ellipse, no background box,
    # white outline so it's readable against any ellipse fill colour
    ax.text(
        cx, cy, str(g),
        ha="center", va="center",
        fontsize=group_fs, fontweight="bold",
        color=color, zorder=6,
        path_effects=[
            pe.withStroke(linewidth=2.5, foreground="white")
        ]
    )

    # Individual sample labels — placed at their PCA coordinates
    for sample, row in grp_data.iterrows():
        ax.text(
            row["PC1"], row["PC2"],
            str(sample),
            ha="center", va="bottom",
            fontsize=label_fs,
            color="#222222",
            zorder=7,
            bbox=label_padding,
            # Offset label slightly above the implicit dot position
            transform=ax.transData,
        )

# Draw thin sample position markers (small +) so positions are clear
# even when labels overlap; no filled dots to avoid visual clutter
for g, grp_data in coords.groupby("Group"):
    color = group_colors[g]
    ax.scatter(
        grp_data["PC1"], grp_data["PC2"],
        marker="+", s=40, color=color,
        linewidths=1.2, zorder=5, alpha=0.8
    )

# ---------------------------------------------------------------------------
# Axes decoration
# ---------------------------------------------------------------------------
ax.set_xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}% variance)",
    fontsize=10, labelpad=8
)
ax.set_ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}% variance)",
    fontsize=10, labelpad=8
)
ax.set_title(
    f"PCA of Weighted KO Abundance  |  {n_samples} samples, {n_groups} groups",
    fontsize=11, pad=12, fontweight="bold"
)

# Light grid — helps read coordinates without cluttering
ax.grid(True, linestyle=":", linewidth=0.5, color="#cccccc", zorder=0)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)

# Add a small colour legend for groups in the corner — compact, no dots,
# just coloured text lines so it takes minimal space
legend_lines = [
    plt.Line2D([0], [0], color=group_colors[g], linewidth=3, label=str(g))
    for g in groups
]
ax.legend(
    handles=legend_lines,
    title=group_col,
    title_fontsize=8,
    fontsize=7.5,
    framealpha=0.7,
    edgecolor="#dddddd",
    loc="best",
    ncol=max(1, n_groups // 8),   # multi-column for many groups
)

# Expand axis limits slightly so ellipses and labels don't clip
x_pad = (coords["PC1"].max() - coords["PC1"].min()) * 0.18 + 0.5
y_pad = (coords["PC2"].max() - coords["PC2"].min()) * 0.18 + 0.5
ax.set_xlim(coords["PC1"].min() - x_pad, coords["PC1"].max() + x_pad)
ax.set_ylim(coords["PC2"].min() - y_pad, coords["PC2"].max() + y_pad)

plt.tight_layout(pad=1.5)
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()

print(f"PCA plot saved -> {output_path}")
print(f"  Samples: {n_samples} | Groups: {n_groups} | KOs used: {ko_abundance.shape[1]}")
print(f"  PC1: {pca.explained_variance_ratio_[0]*100:.1f}%  "
      f"PC2: {pca.explained_variance_ratio_[1]*100:.1f}%")
