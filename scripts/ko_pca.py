import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import numpy as np
import os

# Load ASV table (rows = ASVs, columns = samples)
asv_table = pd.read_csv(snakemake.input["asv_table"], sep="\t", index_col=0)

# Load KO predictions (rows = ASVs, columns = KOs)
ko_table = pd.read_csv(snakemake.input["ko_table"], sep="\t", index_col=0)

# Load metadata
meta = pd.read_csv(snakemake.input["metadata"], sep="\t", index_col=0)

# Keep only ASVs present in both tables
common_asvs = asv_table.index.intersection(ko_table.index)
asv_table = asv_table.loc[common_asvs]
ko_table = ko_table.loc[common_asvs]

# Compute KO abundance per sample: matrix multiplication (ASV abundance × KO presence)
ko_abundance = asv_table.T.dot(ko_table)  # shape: samples × KOs

# Keep only samples present in metadata
common_samples = ko_abundance.index.intersection(meta.index)
ko_abundance = ko_abundance.loc[common_samples]
meta = meta.loc[common_samples]

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(ko_abundance.values)
coords = pd.DataFrame(X_pca, index=ko_abundance.index, columns=["PC1","PC2"])
coords = coords.merge(meta, left_index=True, right_index=True)

# Plot
groups = coords["Group"].unique()
palette = sns.color_palette("tab10", n_colors=len(groups))
group_colors = dict(zip(groups, palette))

fig, ax = plt.subplots(figsize=(8,6))
for g, data in coords.groupby("Group"):
    ax.scatter(data["PC1"], data["PC2"], label=g, color=group_colors[g], s=50)

    # Optional: 95% confidence ellipse
    if len(data) > 2:
        cov = np.cov(data[["PC1","PC2"]].T)
        mean = data[["PC1","PC2"]].mean().values
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        width, height = 2 * np.sqrt(vals) * 2
        ell = Ellipse(xy=mean, width=width, height=height, angle=theta,
                      edgecolor=group_colors[g], facecolor='none', lw=2, ls='--')
        ax.add_patch(ell)

ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
ax.set_title("PCA of Weighted KO Abundance")
ax.legend(title="Group")

os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)
plt.close()

