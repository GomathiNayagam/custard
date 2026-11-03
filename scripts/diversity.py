import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa

asv = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
tax = pd.read_csv(snakemake.input[1], sep="\t")
meta = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)

non_empty = asv.columns[asv.sum(axis=0) > 0]
asv = asv[non_empty]
asv = asv.loc[asv.sum(axis=1) > 0, :]
meta = meta.loc[meta.index.isin(asv.columns), :]

os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)

# Alpha diversity
alpha = pd.DataFrame({
    "Observed": (asv > 0).sum(axis=0),
    "Shannon": alpha_diversity("shannon", asv.T.values, ids=asv.columns)
})
alpha.to_csv(snakemake.output[0], sep="\t")

# Beta diversity (Bray-Curtis)
bc_dm = beta_diversity("braycurtis", asv.T.values, ids=list(asv.columns))
bc_dm.to_data_frame().to_csv(snakemake.output[1], sep="\t")

# PCoA ordination plot
if len(asv.columns) >= 3:
    pc = pcoa(bc_dm)
    props = pc.proportion_explained

    fig, ax = plt.subplots(figsize=(8, 6))

    coords = pc.samples[['PC1', 'PC2']]
    coords.index = list(asv.columns)

    groups = {}
    for sample in coords.index:
        if sample in meta.index and 'Group' in meta.columns:
            g = meta.loc[sample, 'Group']
        else:
            g = 'Unknown'
        groups.setdefault(g, []).append(sample)

    colors = plt.cm.Set2(np.linspace(0, 1, max(len(groups), 3)))

    for i, (group, samples) in enumerate(sorted(groups.items())):
        sub = coords.loc[[s for s in samples if s in coords.index]]
        ax.scatter(sub['PC1'], sub['PC2'], color=colors[i], s=120, label=group,
                   edgecolors='black', linewidth=0.5, zorder=3)
        for s in sub.index:
            label = s.split('_L001')[0]
            ax.annotate(label, (sub.loc[s, 'PC1'], sub.loc[s, 'PC2']),
                        fontsize=8, ha='left', va='bottom', xytext=(5, 5),
                        textcoords='offset points')

    ax.set_xlabel(f"PC1 ({100*props.iloc[0]:.1f}%)", fontsize=12)
    ax.set_ylabel(f"PC2 ({100*props.iloc[1]:.1f}%)", fontsize=12)
    ax.set_title("PCoA (Bray-Curtis)", fontsize=14)
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(snakemake.output[2], dpi=150, bbox_inches='tight')
    plt.close()
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "PCoA requires ≥3 samples", ha='center', va='center',
            transform=ax.transAxes, fontsize=14)
    ax.set_axis_off()
    plt.savefig(snakemake.output[2], dpi=150)
    plt.close()
