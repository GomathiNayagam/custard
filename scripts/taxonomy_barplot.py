import pandas as pd
import matplotlib.pyplot as plt
import os

# Load ASV table, taxonomy, and metadata
asv = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
tax = pd.read_csv(snakemake.input[1], sep="\t", index_col=0)
meta = pd.read_csv(snakemake.input[2], sep="\t", index_col=0)

# Remove empty samples
non_empty_samples = asv.columns[asv.sum(axis=0) > 0]
if len(non_empty_samples) < asv.shape[1]:
    removed_samples = set(asv.columns) - set(non_empty_samples)
    print(f"Removing empty samples: {removed_samples}")
asv = asv[non_empty_samples]

# Align metadata
meta = meta.loc[meta.index.isin(asv.columns), :]

# Merge ASV counts with taxonomy
asv_tax = asv.merge(tax, left_index=True, right_index=True)

# Group by taxonomy level (e.g., Genus)
tax_level = "Genus"  # change as needed
numeric_cols = asv_tax.select_dtypes(include=['number']).columns
grouped = asv_tax.groupby(tax_level)[numeric_cols].sum()

# Remove taxa with zero counts across all samples
grouped = grouped.loc[(grouped.sum(axis=1) > 0), :]

# Collapse rare/Unclassified taxa
threshold = 1.0  # percentage threshold to keep individual taxa
rel_abundance = grouped.div(grouped.sum(axis=0), axis=1) * 100
rare_taxa = rel_abundance.index[(rel_abundance.max(axis=1) < threshold)]
if len(rare_taxa) > 0:
    rel_abundance.loc["Other/Unclassified"] = rel_abundance.loc[rare_taxa].sum()
    rel_abundance = rel_abundance.drop(rare_taxa)

# Reorder columns by metadata
rel_abundance = rel_abundance[meta.index]

# Plot stacked barplot
fig, ax = plt.subplots(figsize=(12,6))
rel_abundance.T.plot(kind="bar", stacked=True, ax=ax, width=0.8, edgecolor="none", color=plt.cm.tab20.colors)

# Optional: dashed lines between groups
if "Group" in meta.columns:
    group_labels = meta["Group"].values
    for i in range(len(group_labels)-1):
        if group_labels[i] != group_labels[i+1]:
            ax.axvline(i + 0.5, color="black", linestyle="--", linewidth=1)

ax.set_ylabel("Relative abundance (%)")
ax.set_xlabel("Samples")
ax.legend(bbox_to_anchor=(1.05,1), loc='upper left', title=tax_level)
plt.tight_layout()

# Save figure
os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
plt.savefig(snakemake.output[0], dpi=300)

