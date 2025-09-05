import pandas as pd
import os

# Load config
configfile: "config.yaml"

# Set default method
if "method" not in config:
    config["method"] = "kraken2"

if config["method"] not in ["kraken2", "dada2"]:
    raise ValueError("config['method'] must be 'kraken2' or 'dada2'")

# Taxonomy output path
TAXONOMY_OUTPUT = (
    f"{config['output_dir']}/asv_kraken2_taxonomy.tsv"
    if config["method"] == "kraken2"
    else f"{config['output_dir']}/asv_dada2_taxonomy.tsv"
)

# Metadata handling
METADATA = config.get("metadata", "none")
SAMPLES = []
if METADATA != "none" and os.path.exists(METADATA):
    meta_df = pd.read_csv(METADATA, sep="\t", comment="#")
    SAMPLES = list(meta_df["SampleID"])
else:
    # Fallback: list all FASTQ basenames
    SAMPLES = [f.split(config["read1_suffix"])[0] for f in os.listdir(config["input_dir"]) if f.endswith(config["read1_suffix"])]

# Utility functions for paired-end reads
def R1(wildcards): return f"{config['input_dir']}/{wildcards.sample}{config['read1_suffix']}"
def R2(wildcards): return f"{config['input_dir']}/{wildcards.sample}{config['read2_suffix']}"

# Workflow endpoints
if METADATA != "none":
    rule all:
        input:
            f"{config['output_dir']}/asv_table.tsv",
            f"{config['output_dir']}/asv_sequences.fasta",
            TAXONOMY_OUTPUT,
            f"{config['output_dir']}/asv_predicted_kos.tsv",
            f"{config['output_dir']}/diversity/alpha_diversity.tsv",
            f"{config['output_dir']}/diversity/beta_diversity.tsv",
            f"{config['output_dir']}/diversity/beta_pcoa.png",
            f"{config['output_dir']}/diversity/rarefaction_curve.png",
            f"{config['output_dir']}/taxonomy/taxonomy_barplot.png",
            f"{config['output_dir']}/ko/ko_pca.png"
else:
    rule all:
        input:
            f"{config['output_dir']}/asv_table.tsv",
            f"{config['output_dir']}/asv_sequences.fasta",
            TAXONOMY_OUTPUT,
            f"{config['output_dir']}/asv_predicted_kos.tsv"

##########################################
# Core rules
##########################################

rule run_dada2:
    input:
        fastq_dir=config["input_dir"]
    output:
        asv_table=f"{config['output_dir']}/asv_table.tsv",
        seqtab_rds=f"{config['output_dir']}/seqtab.rds",
        asv_fasta=f"{config['output_dir']}/asv_sequences.fasta"
    threads: config["threads"]
    script: "scripts/run_dada2.R"

rule kraken2_taxonomy:
    input:
        fasta=f"{config['output_dir']}/asv_sequences.fasta"
    output:
        taxonomy_report=f"{config['output_dir']}/kraken_output.tsv",
        taxonomy_table=f"{config['output_dir']}/asv_kraken2_taxonomy.tsv",
        raw_output=f"{config['output_dir']}/kraken_output_raw.tsv"
    threads: config["threads"]
    shell:
        """
        kraken2 --db {config[kraken_db]} \
                --threads {threads} \
                --use-names \
                --report {output.taxonomy_report} \
                --output {output.raw_output} \
                --confidence 0.1 \
                {input.fasta}
        python scripts/parse_kraken2_taxonomy.py {output.raw_output} {output.taxonomy_table} {config[kraken_db]}
        """

rule dada2_taxonomy:
    input:
        fasta=f"{config['output_dir']}/asv_sequences.fasta",
        seqtab=f"{config['output_dir']}/seqtab.rds"
    output:
        taxonomy=f"{config['output_dir']}/asv_dada2_taxonomy.tsv"
    params:
        taxa_db=config["taxa_db"]
    shell:
        """
        Rscript scripts/assign_dada2_taxonomy.R {input.fasta} {input.seqtab} {output.taxonomy} {params.taxa_db}
        """

rule infer_kos:
    input:
        asv_fasta=f"{config['output_dir']}/asv_sequences.fasta",
        ko_fasta=config["ko_db_fasta"],
        ko_npz=config["ko_npz"],
        ko_genomes=config["ko_genomes"],
        ko_list=config["ko_list"]
    output:
        tsv=f"{config['output_dir']}/asv_predicted_kos.tsv"
    threads: config["threads"]
    shell:
        """
        python scripts/infer_kos.py \
            --asv_fasta {input.asv_fasta} \
            --ko_fasta {input.ko_fasta} \
            --ko_npz {input.ko_npz} \
            --ko_genomes {input.ko_genomes} \
            --ko_list {input.ko_list} \
            --output {output.tsv} \
            --threads {threads} \
            --identity_threshold {config[identity_threshold]}
        """

##########################################
# Extra rules if metadata provided
##########################################

if METADATA != "none":

    rule diversity_metrics:
        input:
            asv_table=f"{config['output_dir']}/asv_table.tsv",
            taxonomy=TAXONOMY_OUTPUT,
            metadata=METADATA
        output:
            alpha=f"{config['output_dir']}/diversity/alpha_diversity.tsv",
            beta=f"{config['output_dir']}/diversity/beta_diversity.tsv",
            beta_pcoa=f"{config['output_dir']}/diversity/beta_pcoa.png"
        script: "scripts/diversity.py"

    rule rarefaction_curve:
        input:
            asv_table=f"{config['output_dir']}/asv_table.tsv",
            metadata=METADATA
        output:
            f"{config['output_dir']}/diversity/rarefaction_curve.png"
        script: "scripts/rarefaction.py"

    rule taxonomy_barplot:
        input:
            asv_table=f"{config['output_dir']}/asv_table.tsv",
            taxonomy=TAXONOMY_OUTPUT,
            metadata=METADATA
        output:
            f"{config['output_dir']}/taxonomy/taxonomy_barplot.png"
        script: "scripts/taxonomy_barplot.py"

    rule ko_pca:
        input:
            asv_table=f"{config['output_dir']}/asv_table.tsv",
            ko_table=f"{config['output_dir']}/asv_predicted_kos.tsv",
            metadata=METADATA
        output:
            f"{config['output_dir']}/ko/ko_pca.png"
        script: "scripts/ko_pca.py"
    
# Prioritise taxonomy rules
ruleorder: kraken2_taxonomy > dada2_taxonomy

