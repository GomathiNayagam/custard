import pandas as pd
import os
import sys
import re
import shlex
from Bio.Seq import Seq

configfile: "config.yaml"

# =============================================================================
# Step 0: Validate config at DAG-build time (catches errors at --dry-run)
# =============================================================================
sys.path.insert(0, os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
from validate_config import validate
validate(config)

# =============================================================================
# Global setup
# =============================================================================
if "method" not in config:
    config["method"] = "dada2"
if config["method"] not in ["kraken2", "dada2"]:
    raise ValueError("config['method'] must be 'kraken2' or 'dada2'")

TAXONOMY_OUTPUT = (
    f"{config['output_dir']}/asv_kraken2_taxonomy.tsv"
    if config["method"] == "kraken2"
    else f"{config['output_dir']}/asv_dada2_taxonomy.tsv"
)

METADATA = config.get("metadata", "none")
SAMPLES  = []
if METADATA != "none" and os.path.exists(METADATA):
    meta_df = pd.read_csv(METADATA, sep="\t", comment="#")
    SAMPLES = list(meta_df["SampleID"])
else:
    SAMPLES = [
        f.split(config["read1_suffix"])[0]
        for f in os.listdir(config["input_dir"])
        if f.endswith(config["read1_suffix"])
    ]

# Sample IDs become shell arguments, file paths, and Snakemake wildcards.
# Restrict them to a safe character set so a name from metadata.tsv can't break
# a shell command (word-splitting, globbing, injection) or a downstream path.
_SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9._-]+$")
_bad_samples = [s for s in SAMPLES if not _SAFE_SAMPLE_ID.match(str(s))]
if _bad_samples:
    raise ValueError(
        "Invalid sample ID(s): "
        + ", ".join(repr(s) for s in _bad_samples)
        + ". Sample IDs may contain only letters, digits, '.', '_' and '-' "
        "(no spaces or shell metacharacters)."
    )

# Sanity check thresholds (with defaults if not in config)
SC            = config.get("sanity_check", {})
WARN_LOSS_PCT = SC.get("warn_loss_pct", 80)
FAIL_LOSS_PCT = SC.get("fail_loss_pct", 90)

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

PRIMER_F    = config["cutadapt"]["primer_f"]
PRIMER_R    = config["cutadapt"]["primer_r"]
RC_PRIMER_F = reverse_complement(PRIMER_F)
RC_PRIMER_R = reverse_complement(PRIMER_R)

import gzip

def _iupac_to_regex(seq):
    iupac = {
        'A':'A','C':'C','G':'G','T':'T',
        'R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]',
        'K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]',
        'H':'[ACT]','V':'[ACG]','N':'[ACGT]'
    }
    return ''.join(iupac.get(b.upper(), b) for b in seq)

def reads_have_primers(input_dir, read1_suffix, primer_f, n_check=200):
    files = sorted([
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(read1_suffix)
    ])
    if not files:
        return True
    f = files[0]
    primer_re = re.compile(_iupac_to_regex(primer_f[:12]))
    hits, total = 0, 0
    opener = gzip.open if f.endswith('.gz') else open
    with opener(f, 'rt') as fh:
        for _ in range(n_check):
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline(); fh.readline()
            total += 1
            if primer_re.match(seq):
                hits += 1
    frac = hits / total if total > 0 else 0
    print(f"  Primer detection: {hits}/{total} reads start with primer ({frac:.0%})")
    return frac > 0.10

PRIMERS_PRESENT = reads_have_primers(
    config["input_dir"], config["read1_suffix"], PRIMER_F
)
if PRIMERS_PRESENT:
    print("Primers detected in input FASTQs -- cutadapt will trim them.")
else:
    print("WARNING: Primers NOT detected -- cutadapt RC trimming only.")

def R1(wildcards):
    return f"{config['input_dir']}/{wildcards.sample}{config['read1_suffix']}"
def R2(wildcards):
    return f"{config['input_dir']}/{wildcards.sample}{config['read2_suffix']}"

TRIMMED_DIR  = f"{config['output_dir']}/trimmed"
FILTERED_DIR = f"{config['output_dir']}/filtered"   # DADA2 filterAndTrim output dir

# =============================================================================
# Workflow endpoints
# =============================================================================
_core = [
    f"{config['output_dir']}/asv_table.tsv",
    f"{config['output_dir']}/asv_sequences.fasta",
    TAXONOMY_OUTPUT,
    f"{config['output_dir']}/asv_predicted_kos.tsv",
    f"{config['output_dir']}/read_tracking_qc.tsv",
    f"{config['output_dir']}/discarded/discarded_reads_taxonomy.tsv",
    f"{config['output_dir']}/asv_vsearch_taxonomy.tsv",
    f"{config['output_dir']}/asv_table_annotated.tsv",
    f"{config['output_dir']}/report.html",
    f"{config['output_dir']}/taxonomy/krona.html",
]
_diversity = [
    f"{config['output_dir']}/diversity/alpha_diversity.tsv",
    f"{config['output_dir']}/diversity/beta_diversity.tsv",
    f"{config['output_dir']}/diversity/beta_pcoa.png",
    f"{config['output_dir']}/diversity/rarefaction_curve.png",
    f"{config['output_dir']}/taxonomy/taxonomy_barplot.png",
    f"{config['output_dir']}/ko/ko_pca.png",
]

if METADATA != "none":
    rule all:
        input: _core + _diversity
else:
    rule all:
        input: _core

# =============================================================================
# Setup rule — one-time database download from Zenodo
# -----------------------------------------------------------------------------
# Downloads a single db.tar.gz from Zenodo containing all CUSTARD databases.
# The tarball is expected to extract to a  db/  subdirectory containing:
#
#   db/dada2_formatted_rpoc_amplicon.fasta
#   db/kraken2_amplicons.fasta
#   db/all_genomes_kos.npz
#   db/all_genomes_kos_genomes.txt
#   db/all_genomes_kos_kos.txt
#   db/kraken_amplicon/hash.k2d  (and other Kraken2 files)
#
# Usage (first-time only):
#   snakemake --cores 4 setup
#
# The sentinel file data/db/.setup_complete is written on success.
# Re-running setup will skip if the sentinel exists (Snakemake DAG).
# To force re-download: rm data/db/.setup_complete && snakemake --cores 4 setup
# =============================================================================
rule setup:
    output:
        sentinel = "data/db/.setup_complete"
    params:
        zenodo_url = config.get(
            "download_base_url",
            "https://zenodo.org/records/17065095/files/db.tar.gz?download=1"
        ),
        db_dir = "data/db",
    shell:
        """
        set -euo pipefail
        mkdir -p {params.db_dir}

        echo "╔══ CUSTARD Database Setup ══════════════════════════════════"
        echo "║  Downloading from Zenodo..."
        echo "║  URL: {params.zenodo_url}"
        echo "╚════════════════════════════════════════════════════════════"
        echo ""

        TARBALL="{params.db_dir}/db.tar.gz"

        # Download
        wget --progress=bar:force:noscroll \
             -O "$TARBALL" \
             "{params.zenodo_url}"

        echo ""
        echo "Download complete. Extracting..."

        # Extract — tarball is expected to produce a db/ subdirectory
        # Strip leading db/ path component and place files into {params.db_dir}
        tar -xzf "$TARBALL" \
            --strip-components=1 \
            -C "{params.db_dir}/"

        rm "$TARBALL"
        echo "Tarball removed."

        # Validate expected files are present
        MISSING=0
        for f in \
            rpocdb_v3_dada2_toGenus.fasta \
            rpocdb_v3_addSpecies.fasta \
            rpocdb_v3_vsearch.fasta \
            kraken2_amplicons.fasta \
            all_genomes_kos.npz \
            all_genomes_kos_genomes.txt \
            all_genomes_kos_kos.txt
        do
            if [ ! -f "{params.db_dir}/$f" ]; then
                echo "  ✗ Missing after extract: $f"
                MISSING=$((MISSING+1))
            else
                echo "  ✓ $f"
            fi
        done

        if [ "$MISSING" -gt 0 ]; then
            echo ""
            echo "ERROR: $MISSING expected files missing after extraction."
            exit 1
        fi

        # Build Kraken2 database from vsearch FASTA + NCBI taxonomy
        if [ ! -f "{params.db_dir}/kraken_amplicon/hash.k2d" ]; then
            echo ""
            echo "Building Kraken2 database (one-time, requires internet)..."
            python scripts/build_kraken2_db.py \
                --db-dir {params.db_dir} \
                --threads {config[threads]}
        else
            echo "  ✓ kraken_amplicon/hash.k2d (Kraken2 DB exists)"
        fi

        echo ""
        echo "╔══ Setup complete ══════════════════════════════════════════"
        echo "║  All databases ready in {params.db_dir}/"
        echo "╚════════════════════════════════════════════════════════════"
        date > {output.sentinel}
        """

# =============================================================================
# Rule 1: Cutadapt — primer trimming
# =============================================================================
rule cutadapt_trim:
    input:
        r1 = R1,
        r2 = R2
    output:
        r1  = f"{TRIMMED_DIR}/{{sample}}{config['read1_suffix']}",
        r2  = f"{TRIMMED_DIR}/{{sample}}{config['read2_suffix']}",
        log = f"{TRIMMED_DIR}/{{sample}}_cutadapt.log"
    params:
        primer_f        = lambda wc: PRIMER_F,
        primer_r        = lambda wc: PRIMER_R,
        rc_primer_f     = lambda wc: RC_PRIMER_F,
        rc_primer_r     = lambda wc: RC_PRIMER_R,
        error_rate      = lambda wc: config["cutadapt"]["error_rate"],
        min_length      = lambda wc: config["cutadapt"]["minimum_length"],
        discard_flag    = lambda wc: "--discard-untrimmed" if PRIMERS_PRESENT else "",
        fwd_primer_args = lambda wc: (
            f"-g {PRIMER_F} -G {PRIMER_R}" if PRIMERS_PRESENT else ""
        )
    threads: config["threads"]
    shell:
        """
        cutadapt \
            {params.fwd_primer_args} \
            -a {params.rc_primer_r} \
            -A {params.rc_primer_f} \
            {params.discard_flag} \
            --error-rate {params.error_rate} \
            --minimum-length {params.min_length} \
            --cores {threads} \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2} \
            > {output.log} 2>&1
        """

# =============================================================================
# Rule 2: DADA2 — ASV inference
# =============================================================================
rule run_dada2:
    input:
        trimmed_r1 = expand(
            f"{TRIMMED_DIR}/{{sample}}{config['read1_suffix']}", sample=SAMPLES
        ),
        trimmed_r2 = expand(
            f"{TRIMMED_DIR}/{{sample}}{config['read2_suffix']}", sample=SAMPLES
        )
    output:
        asv_table     = f"{config['output_dir']}/asv_table.tsv",
        seqtab_rds    = f"{config['output_dir']}/seqtab.rds",
        asv_fasta     = f"{config['output_dir']}/asv_sequences.fasta",
        read_tracking = f"{config['output_dir']}/read_tracking.tsv",
        qc_plots      = f"{config['output_dir']}/qc_plots.pdf"
    params:
        fastq_dir    = TRIMMED_DIR,
        filtered_dir = FILTERED_DIR,
        amplicon_len = config["dada2"]["amplicon_len"],
        primer_len_f = config["dada2"]["primer_len_f"],
        primer_len_r = config["dada2"]["primer_len_r"],
        maxEE_F      = config["dada2"]["maxEE_F"],
        maxEE_R      = config["dada2"]["maxEE_R"],
        pool         = config["dada2"]["pool"],
        min_overlap  = config["dada2"]["min_overlap"],
        read_length  = config["dada2"]["read_length"],
        truncLen_F   = config["dada2"].get("truncLen_F", None),
        truncLen_R   = config["dada2"].get("truncLen_R", None)
    threads: config["threads"]
    script: "scripts/run_dada2.R"

# =============================================================================
# Rule 3: Read tracking sanity check
# -----------------------------------------------------------------------------
# Runs immediately after DADA2. Checks for catastrophic or high read loss
# at each step. Writes a colour-coded retention table used by the HTML report.
# Exit code 1 = hard failure (all samples failing); 2 = warnings only.
# =============================================================================
rule check_read_tracking:
    input:
        tracking = f"{config['output_dir']}/read_tracking.tsv"
    output:
        qc_report = f"{config['output_dir']}/read_tracking_qc.tsv"
    params:
        warn_loss_pct = WARN_LOSS_PCT,
        fail_loss_pct = FAIL_LOSS_PCT,
    shell:
        """
        python scripts/check_read_tracking.py \
            --tracking       {input.tracking} \
            --output         {output.qc_report} \
            --warn_loss_pct  {params.warn_loss_pct} \
            --fail_loss_pct  {params.fail_loss_pct} \
        || EXIT=$?
        # Exit 2 = warnings only — allow pipeline to continue
        if [ "${{EXIT:-0}}" -eq 2 ]; then
            echo "Read tracking: warnings present but continuing."
            exit 0
        fi
        exit ${{EXIT:-0}}
        """

# =============================================================================
# Rule 4: Extract and classify filterAndTrim-discarded reads
# -----------------------------------------------------------------------------
# Diffs cutadapt-trimmed FASTQs against DADA2 filtered FASTQs to recover
# reads removed by filterAndTrim, then classifies them with Kraken2.
# Output tells you whether the discarded reads contain real signal (e.g.
# low-quality but correctly-targeted reads) or off-target contamination.
# =============================================================================
rule classify_discarded_reads:
    input:
        trimmed   = expand(
            f"{TRIMMED_DIR}/{{sample}}{config['read1_suffix']}", sample=SAMPLES
        ),
        # DADA2 filterAndTrim writes passing reads to filtered_dir
        # This rule depends on run_dada2 having completed
        asv_table = f"{config['output_dir']}/asv_table.tsv"
    output:
        fasta    = f"{config['output_dir']}/discarded/discarded_reads.fasta",
        taxonomy = f"{config['output_dir']}/discarded/discarded_reads_taxonomy.tsv",
        report   = f"{config['output_dir']}/discarded/kraken_discarded_report.tsv",
        raw      = f"{config['output_dir']}/discarded/kraken_discarded_raw.tsv"
    params:
        trimmed_dir  = TRIMMED_DIR,
        filtered_dir = FILTERED_DIR,
        samples      = " ".join(shlex.quote(s) for s in SAMPLES),
        read1_suffix = config["read1_suffix"],
        out_dir      = f"{config['output_dir']}/discarded"
    threads: config["threads"]
    shell:
        """
        mkdir -p {params.out_dir}

        # Step 1: Extract discarded reads into a merged FASTA
        python scripts/extract_discarded_reads.py \
            --trimmed_dir  {params.trimmed_dir} \
            --filtered_dir {params.filtered_dir} \
            --output_fasta {output.fasta} \
            --read1_suffix {params.read1_suffix} \
            --samples      {params.samples}

        # Step 2: Classify with Kraken2
        if [ -s {output.fasta} ]; then
            kraken2 --db {config[kraken_db]} \
                    --threads {threads} \
                    --use-names \
                    --report  {output.report} \
                    --output  {output.raw} \
                    --confidence 0.1 \
                    {output.fasta}
            python scripts/parse_kraken2_taxonomy.py \
                {output.raw} \
                {output.taxonomy} \
                {config[kraken_db]}
        else
            echo "No discarded reads found (all reads passed filterAndTrim)."
            touch {output.report} {output.raw}
            echo -e "ASV\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies" \
                > {output.taxonomy}
        fi
        """

# =============================================================================
# Rule 5: Kraken2 taxonomy (ASVs)
# =============================================================================
rule kraken2_taxonomy:
    input:
        fasta = f"{config['output_dir']}/asv_sequences.fasta"
    output:
        taxonomy_report = f"{config['output_dir']}/kraken_output.tsv",
        taxonomy_table  = f"{config['output_dir']}/asv_kraken2_taxonomy.tsv",
        raw_output      = f"{config['output_dir']}/kraken_output_raw.tsv"
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
        python scripts/parse_kraken2_taxonomy.py \
            {output.raw_output} \
            {output.taxonomy_table} \
            {config[kraken_db]}
        """

# =============================================================================
# Rule 6: DADA2 taxonomy (fallback)
# =============================================================================
rule dada2_taxonomy:
    input:
        fasta  = f"{config['output_dir']}/asv_sequences.fasta",
        seqtab = f"{config['output_dir']}/seqtab.rds"
    output:
        taxonomy = f"{config['output_dir']}/asv_dada2_taxonomy.tsv"
    params:
        taxa_db = config["taxa_db"],
        addsp_db = config.get("addspecies_db", ""),
        # min_boot lived in config but the rule never passed it. The script then
        # used its own hard-coded 50, so editing the config had no effect.
        min_boot = config.get("dada2", {}).get("min_boot", 50)
    shell:
        """
        Rscript scripts/assign_dada2_taxonomy.R \
            {input.fasta} {input.seqtab} {output.taxonomy} {params.taxa_db} \
            "{params.addsp_db}" {params.min_boot}
        """

# =============================================================================
# Rule 7: KO inference — weighted top-N vsearch + batch matrix indexing
# =============================================================================
rule infer_kos:
    input:
        asv_fasta  = f"{config['output_dir']}/asv_sequences.fasta",
        ko_fasta   = config["ko_db_fasta"],
        ko_npz     = config["ko_npz"],
        ko_genomes = config["ko_genomes"],
        ko_list    = config["ko_list"]
    output:
        tsv = f"{config['output_dir']}/asv_predicted_kos.tsv"
    threads: config["threads"]
    shell:
        """
        python scripts/infer_kos.py \
            --asv_fasta          {input.asv_fasta} \
            --ko_fasta           {input.ko_fasta} \
            --ko_npz             {input.ko_npz} \
            --ko_genomes         {input.ko_genomes} \
            --ko_list            {input.ko_list} \
            --output             {output.tsv} \
            --threads            {threads} \
            --identity_threshold {config[identity_threshold]} \
            --min_coverage       {config[min_coverage]} \
            --top_n_hits         {config[ko_inference][top_n_hits]}
        """

# =============================================================================
# Rule 8: Vsearch species-level taxonomy
# -----------------------------------------------------------------------------
# Runs vsearch against the DADA2 reference FASTA at high identity (97% default)
# with strict coverage (90% default) to achieve species-level resolution.
# Takes top-N hits, resolves consensus taxonomy by weighted majority vote.
# Output: asv_vsearch_taxonomy.tsv with vsearch_domain...vsearch_species columns.
# =============================================================================
rule vsearch_species_taxonomy:
    input:
        asv_fasta = f"{config['output_dir']}/asv_sequences.fasta",
        ref_fasta = config.get("vsearch_db", config["taxa_db"])
    output:
        tsv = f"{config['output_dir']}/asv_vsearch_taxonomy.tsv"
    threads: config["threads"]
    shell:
        """
        python scripts/vsearch_species.py \
            --asv_fasta    {input.asv_fasta} \
            --ref_fasta    {input.ref_fasta} \
            --output       {output.tsv} \
            --threads      {threads} \
            --identity     {config[vsearch_taxonomy][identity]} \
            --min_coverage {config[vsearch_taxonomy][min_coverage]} \
            --top_n        {config[vsearch_taxonomy][top_n]}
        """

# =============================================================================
# Rule 9: Annotated ASV table
# -----------------------------------------------------------------------------
# Merges all taxonomy sources into one wide table with column order:
#   ASV | domain…species (primary) | vsearch_domain…vsearch_species |
#   vsearch_identity | vsearch_top_hit | vsearch_n_hits |
#   <sample1> … <sampleN> | total
# Sorted by total reads descending. Original asv_table.tsv is unchanged.
# =============================================================================
rule annotate_asv_table:
    input:
        asv_table        = f"{config['output_dir']}/asv_table.tsv",
        primary_taxonomy = TAXONOMY_OUTPUT,
        vsearch_taxonomy = f"{config['output_dir']}/asv_vsearch_taxonomy.tsv"
    output:
        annotated = f"{config['output_dir']}/asv_table_annotated.tsv"
    params:
        method = config["method"]
    shell:
        """
        python scripts/annotate_asv_table.py \
            --asv_table        {input.asv_table} \
            --primary_taxonomy {input.primary_taxonomy} \
            --vsearch_taxonomy {input.vsearch_taxonomy} \
            --output           {output.annotated} \
            --method           {params.method}
        """

# =============================================================================
# Rule 10: HTML summary report — always generated
# =============================================================================
rule generate_report:
    input:
        asv_table        = f"{config['output_dir']}/asv_table.tsv",
        taxonomy         = TAXONOMY_OUTPUT,
        read_tracking    = f"{config['output_dir']}/read_tracking.tsv",
        read_tracking_qc = f"{config['output_dir']}/read_tracking_qc.tsv",
        asv_fasta        = f"{config['output_dir']}/asv_sequences.fasta",
        discarded_tax    = f"{config['output_dir']}/discarded/discarded_reads_taxonomy.tsv",
        ko_table         = f"{config['output_dir']}/asv_predicted_kos.tsv",
        qc_plots         = f"{config['output_dir']}/qc_plots.pdf",
        config_yaml      = "config.yaml"
    output:
        html = f"{config['output_dir']}/report.html"
    params:
        metadata = METADATA if METADATA != "none" else ""
    shell:
        """
        python scripts/generate_report.py \
            --asv_table          {input.asv_table} \
            --taxonomy           {input.taxonomy} \
            --read_tracking      {input.read_tracking} \
            --read_tracking_qc   {input.read_tracking_qc} \
            --asv_fasta          {input.asv_fasta} \
            --discarded_taxonomy {input.discarded_tax} \
            --ko_table           {input.ko_table} \
            --qc_plots_pdf       {input.qc_plots} \
            --config_yaml        {input.config_yaml} \
            --output             {output.html} \
            $([ -n "{params.metadata}" ] && echo "--metadata {params.metadata}" || true)
        """

# =============================================================================
# Rules 11–14: Diversity and visualisation (metadata required)
# =============================================================================
if METADATA != "none":

    rule diversity_metrics:
        input:
            asv_table = f"{config['output_dir']}/asv_table.tsv",
            taxonomy  = TAXONOMY_OUTPUT,
            metadata  = METADATA
        output:
            alpha     = f"{config['output_dir']}/diversity/alpha_diversity.tsv",
            beta      = f"{config['output_dir']}/diversity/beta_diversity.tsv",
            beta_pcoa = f"{config['output_dir']}/diversity/beta_pcoa.png"
        script: "scripts/diversity.py"

    rule rarefaction_curve:
        input:
            asv_table = f"{config['output_dir']}/asv_table.tsv",
            metadata  = METADATA
        output:
            f"{config['output_dir']}/diversity/rarefaction_curve.png"
        script: "scripts/rarefaction.py"

    rule taxonomy_barplot:
        input:
            asv_table = f"{config['output_dir']}/asv_table.tsv",
            taxonomy  = TAXONOMY_OUTPUT,
            metadata  = METADATA
        output:
            f"{config['output_dir']}/taxonomy/taxonomy_barplot.png"
        script: "scripts/taxonomy_barplot.py"

    rule ko_pca:
        input:
            asv_table = f"{config['output_dir']}/asv_table.tsv",
            ko_table  = f"{config['output_dir']}/asv_predicted_kos.tsv",
            metadata  = METADATA
        output:
            f"{config['output_dir']}/ko/ko_pca.png"
        script: "scripts/ko_pca.py"

# =============================================================================
# Krona interactive taxonomy chart
# =============================================================================
rule krona_chart:
    input:
        asv_table = f"{config['output_dir']}/asv_table.tsv",
        taxonomy  = TAXONOMY_OUTPUT
    output:
        krona = f"{config['output_dir']}/taxonomy/krona.html"
    shell:
        """
        python scripts/krona_chart.py {input.asv_table} {input.taxonomy} {output.krona}
        """

# =============================================================================
# Rule ordering
# =============================================================================
ruleorder: kraken2_taxonomy > dada2_taxonomy
