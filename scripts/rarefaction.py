import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)

asv = pd.read_csv(snakemake.input["asv_table"], sep="\t", index_col=0)

# ASV table: rows=ASVs, cols=samples. Transpose to rows=samples, cols=ASVs
asv = asv.T

def rarefaction_curve(counts, n_iters=10, step=100):
    total = int(counts.sum())
    if total == 0:
        return [], []
    depths = list(range(step, total + 1, step))
    if depths[-1] != total:
        depths.append(total)
    richness = []
    pool = np.repeat(np.arange(len(counts)), counts.astype(int))
    for d in depths:
        observed = []
        for _ in range(n_iters):
            subsample = np.random.choice(pool, min(d, len(pool)), replace=False)
            observed.append(len(np.unique(subsample)))
        richness.append(np.mean(observed))
    return depths, richness

fig, ax = plt.subplots(figsize=(10, 6))

colors = plt.cm.tab10(np.linspace(0, 1, len(asv.index)))

for i, sample in enumerate(asv.index):
    counts = asv.loc[sample].values
    if counts.sum() == 0:
        continue
    depths, richness = rarefaction_curve(counts)
    ax.plot(depths, richness, label=sample, color=colors[i], linewidth=2)

ax.set_xlabel("Sequencing depth", fontsize=12)
ax.set_ylabel("Observed ASVs", fontsize=12)
ax.set_title("Rarefaction curves", fontsize=14)
ax.legend(loc='lower right', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=150, bbox_inches='tight')
plt.close()
