import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from skbio.diversity import alpha_diversity

# Load ASV table and metadata
asv = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
meta = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)

# Filter out empty samples (columns with zero sum)
non_empty_samples = asv.columns[asv.sum(axis=0) > 0]
if len(non_empty_samples) < asv.shape[1]:
    removed_samples = set(asv.columns) - set(non_empty_samples)
    print(f"Removing empty samples: {removed_samples}")
asv = asv[non_empty_samples]

# Filter ASVs with zero counts across all retained samples
asv = asv.loc[asv.sum(axis=1) > 0, :]

# Ensure metadata matches the retained samples
meta = meta.loc[meta.index.isin(asv.columns), :]

# Alpha diversity (Observed and Shannon)
alpha = pd.DataFrame({
    "Observed": (asv > 0).sum(axis=0),
    "Shannon": alpha_diversity("shannon", asv.T.values, ids=asv.columns)
})
alpha.to_csv(snakemake.output[0], sep="\t")

# Rarefaction curves per sample
def rarefaction_counts(counts, max_depth=None, steps=20):
    """Generate rarefaction points for a single sample"""
    counts = counts[counts > 0]
    n_reads = counts.sum()
    if max_depth is None:
        max_depth = n_reads
    depths = np.linspace(1, max_depth, steps, dtype=int)
    obs_asvs = []
    for d in depths:
        if d >= n_reads:
            obs_asvs.append(len(counts))
            continue
        sampled = np.random.choice(counts.index.repeat(counts.values), size=d, replace=False)
        obs_asvs.append(len(np.unique(sampled)))
    return depths, obs_asvs

fig, ax = plt.subplots(figsize=(10,6))
for sample in asv.columns:
    depths, obs = rarefaction_counts(asv[sample], steps=50)
    ax.plot(depths, obs, label=sample)
ax.set_xlabel("Sequencing depth (reads)")
ax.set_ylabel("Observed ASVs")
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
plt.tight_layout()
os.makedirs(os.path.dirname(snakemake.output[2]), exist_ok=True)
plt.savefig(snakemake.output[2], dpi=300)

# Save beta diversity matrix for completeness (optional)
from skbio.diversity import beta_diversity
bc_dist = beta_diversity("braycurtis", asv.T.values, ids=asv.columns)
bc_dist.to_data_frame().to_csv(snakemake.output[1], sep="\t")

