# Snakefile

configfile: "config.yaml"

rule all:
    input:
        "results/asv_table.tsv",
        "results/taxa_assignments.tsv",
        "results/seqtab.rds",
        "results/taxa.rds"

rule run_dada2:
    input:
        fastq_dir = config["input_dir"],
        taxa_db = config["taxa_db"]
    output:
        asv_table = "results/asv_table.tsv",
        taxa = "results/taxa_assignments.tsv",
        seqtab_rds = "results/seqtab.rds",
        taxa_rds = "results/taxa.rds"
    threads: config["threads"]
    script:
        "scripts/run_dada2.R"
