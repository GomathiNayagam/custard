import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Ensure output directory exists
os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)

# Load ASV table (samples as rows)
asv = pd.read_csv(snakemake.input["asv_table"], sep="\t", index_col=0)
asv = asv.T  # transpose if ASVs are rows

# Function to compute rarefaction curve
def rarefaction_curve(counts, step=100):
    """counts: array of integer counts of ASVs in a sample"""
    total = counts.sum()
    depths = np.arange(step, total+1, step)
    richness = []
    for d in depths:
        if d > total:
            richness.append(len(counts[counts>0]))
        else:
            # subsample without replacement
            subsample = np.random.choice(np.repeat(np.arange(len(counts)), counts), d, replace=False)
            richness.append(len(np.unique(subsample)))
    return depths, richness

# Plot rarefaction for each sample
plt.figure(figsize=(12,6))
for sample in asv.columns:
    counts = asv[sample].values
    if counts.sum() == 0:
        continue
    depths, richness = rarefaction_curve(counts)
    plt.plot(depths, richness, label=sample)

plt.xlabel("Sequencing depth")
plt.ylabel("Observed ASVs")
plt.title("Rarefaction curves")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)

