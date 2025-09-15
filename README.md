# CUSTARD: DADA2-based Amplicon Processing Tool

**CUSTARD** is a reproducible pipeline for processing paired-end amplicon sequencing reads using the DADA2 workflow, specifically tailored to assign taxonomy using a custom `rpoC` gene-based reference database. This tool constitutes the first module in a two-stage system, with future provision for functional state prediction based on representative genome annotations.

---

custard/
├── Snakefile
├── config.yaml
├── scripts/
│   ├── run_dada2.R
│   ├── taxonomy_barplot.py
│   ├── ko_pca.py
│   ├── diversity_metrics.py
│   └── ...
├── data/
│   └── fastq/                  # Place your raw FASTQ files here
├── db/
│   └── rpoc_dada2_taxonomy.fa  # Custom taxonomy file (for DADA2 mode)
└── results/
    ├── asv_table.tsv
    ├── asv_sequences.fasta
    ├── asv_kraken2_taxonomy.tsv
    ├── taxonomy/
    │   └── taxonomy_barplot.png
    ├── ko/
    │   └── ko_pca.png
    └── diversity/
        ├── alpha_diversity.tsv
        ├── beta_diversity.tsv
        └── rarefaction_curves.png

---

🛠️ Installation  
**System dependencies (Ubuntu)**  
```bash
sudo apt update && sudo apt install -y     build-essential     r-base     libxml2-dev     libcurl4-openssl-dev     libssl-dev     libfontconfig1-dev     libfreetype6-dev     libpng-dev     libharfbuzz-dev     libfribidi-dev     libtiff5-dev     libjpeg-dev     zlib1g-dev     libbz2-dev     pkg-config kraken2
```

**Conda environment**  
```bash
conda create -n custard snakemake -c bioconda -c conda-forge
conda activate custard
pip install ete3 scikit-bio
```

**R dependencies**  
Launch R and install:  
```R
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("dada2", "Biostrings"))
```

---

## Installation Help
Run the following commands to set up the environment:

### Download Database
Download and unpack the database tar file from Zenodo and move the extracted contents into /data/db:

```bash
wget https://zenodo.org/records/17065095/files/db.tar.gz?download=1 -O db.tar.gz
tar -xzvf db.tar.gz

```

---


⚙️ **Configuration (`config.yaml`)**  
```yaml
input_dir: "data/fastq"
output_dir: "results"
taxa_db: "db/rpoc_dada2_taxonomy.fa"
threads: 8
method: "kraken2"   # Options: kraken2 (default) or dada2
metadata: "metadata.tsv"   # Optional, required for intermediate use
```

---

🚀 **Running the Pipeline**  

**Beginner (default, Kraken2 mode)**  
Place raw paired-end FASTQ files (with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes) in `data/fastq/`, then run:  
```bash
cd custard
snakemake --cores 8
```
This will run the default Kraken2 workflow.  
No metadata or custom taxonomy is required.  

**Intermediate (DADA2 + metadata)**  
For more control and reproducibility, use the DADA2 workflow.  

Provide a custom taxonomy reference (e.g., `db/rpoc_dada2_taxonomy.fa`).  

Add a metadata file (`metadata.tsv`) with the following format:  
```
SampleID    Group
Sample1     Control
Sample2     Treatment
Sample3     Treatment
```

Run the pipeline in DADA2 mode:  
```bash
cd custard
snakemake --cores 4 --config method=dada2 metadata=metadata.tsv
```

---

🧾 **Output**  

All results are saved in `results/`:  

**Shared outputs**  
- `asv_table.tsv` – ASV abundance table  
- `asv_sequences.fasta` – Representative ASV sequences  

**Kraken2 mode**  
- `asv_kraken2_taxonomy.tsv` – Taxonomic assignment via Kraken2  

**DADA2 mode**  
- `taxa_assignments.tsv` – Taxonomic assignment using custom database  
- `seqtab.rds` – Sequence table (RDS format)  
- `taxa.rds` – Raw taxonomy assignments (RDS format)  

---

🔮 **Next Steps (Planned in Part 2)**  
- Match ASVs to representative genomes using rpoC hits  
- Retrieve precomputed EGGNOG annotations from representative genomes  
- Infer unobserved functional states and generate KO/COG profiles for each sample  

---

📄 **License**  
This repository is distributed for academic and research use. No warranty is provided.  

✉️ **Contact**  
For queries or contributions, please contact the pipeline author.  
