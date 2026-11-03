"""Auto-generate config.yaml from CLI arguments and detected file patterns."""

import os
import glob
import yaml


SUFFIX_PATTERNS = [
    ("_L001_R1_001.fastq.gz", "_L001_R2_001.fastq.gz"),
    ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
    ("_R1_001.fq.gz", "_R2_001.fq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
]


def detect_read_suffixes(input_dir):
    """Scan input directory and detect the R1/R2 suffix pattern."""
    files = os.listdir(input_dir)
    gz_files = [f for f in files if f.endswith(('.fastq.gz', '.fq.gz'))]

    if not gz_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir}")

    for r1_suf, r2_suf in SUFFIX_PATTERNS:
        r1_count = sum(1 for f in gz_files if f.endswith(r1_suf))
        r2_count = sum(1 for f in gz_files if f.endswith(r2_suf))
        if r1_count > 0 and r1_count == r2_count:
            return r1_suf, r2_suf

    raise ValueError(
        f"Could not detect paired-end read suffixes in {input_dir}.\n"
        f"Found files like: {gz_files[:3]}\n"
        f"Expected patterns: {[p[0] for p in SUFFIX_PATTERNS]}"
    )


def find_metadata(input_dir):
    """Look for metadata.tsv in common locations."""
    candidates = [
        os.path.join(input_dir, "metadata.tsv"),
        os.path.join(input_dir, "..", "metadata.tsv"),
        "metadata.tsv",
        os.path.join("config", "metadata.tsv"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return os.path.abspath(path)
    return None


def validate_metadata(path):
    """Check that metadata has Group column filled in. Returns (ok, message)."""
    if not os.path.isfile(path):
        return False, f"File not found: {path}"

    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]

    if len(lines) < 2:
        return False, "Metadata file is empty."

    header = lines[0].split('\t')
    if 'SampleID' not in header:
        return False, "Metadata file missing 'SampleID' column."
    if 'Group' not in header:
        return False, "Metadata file missing 'Group' column."

    group_idx = header.index('Group')
    empty_groups = []
    for line in lines[1:]:
        fields = line.split('\t')
        group = fields[group_idx].strip() if group_idx < len(fields) else ''
        if not group:
            sample = fields[0] if fields else '?'
            empty_groups.append(sample)

    if empty_groups:
        return False, (
            f"The 'Group' column is empty for: {', '.join(empty_groups)}\n"
            f"Please edit {path} and fill in the Group for each sample."
        )

    return True, "OK"


def generate_metadata_template(input_dir, r1_suffix, output_path="metadata.tsv"):
    """Generate a metadata template with SampleID pre-filled, Group left for user.

    Returns the path to the generated file, or None if file already exists.
    """
    if os.path.exists(output_path):
        return None

    samples = sorted([
        f.replace(r1_suffix, '')
        for f in os.listdir(input_dir)
        if f.endswith(r1_suffix)
    ])

    if not samples:
        return None

    with open(output_path, 'w') as f:
        f.write("SampleID\tGroup\n")
        for s in samples:
            f.write(f"{s}\t\n")

    return os.path.abspath(output_path)


def count_samples(input_dir, r1_suffix):
    """Count how many samples are in the input directory."""
    files = [f for f in os.listdir(input_dir) if f.endswith(r1_suffix)]
    return len(files)


def find_pipeline_dir():
    """Locate the CUSTARD pipeline directory (where Snakefile lives)."""
    candidates = [
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        os.getcwd(),
    ]
    for d in candidates:
        if os.path.isfile(os.path.join(d, "Snakefile")):
            return d
    raise FileNotFoundError(
        "Cannot find CUSTARD Snakefile. Are you running from the CUSTARD directory?"
    )


def build_config(input_dir, output_dir="results", metadata=None, threads=4,
                 method="dada2", pipeline_dir=None):
    """Build a config dictionary for the CUSTARD pipeline."""
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)

    if pipeline_dir is None:
        pipeline_dir = find_pipeline_dir()

    r1_suf, r2_suf = detect_read_suffixes(input_dir)
    n_samples = count_samples(input_dir, r1_suf)

    if metadata is None:
        metadata = find_metadata(input_dir)

    db_dir = os.path.join(pipeline_dir, "data", "db")

    config = {
        "input_dir": input_dir,
        "output_dir": output_dir,

        "download_base_url": "https://zenodo.org/records/17065095/files/db.tar.gz?download=1",

        "taxa_db": os.path.join(db_dir, "rpocdb_v3.1_dada2_toGenus.fasta"),
        "addspecies_db": os.path.join(db_dir, "rpocdb_v3.1_addSpecies.fasta"),
        "vsearch_db": os.path.join(db_dir, "rpocdb_v3.1_vsearch.fasta"),
        "kraken_db": os.path.join(db_dir, "kraken_amplicon"),
        "ko_db_fasta": os.path.join(db_dir, "kraken2_amplicons.fasta"),
        "ko_npz": os.path.join(db_dir, "all_genomes_kos.npz"),
        "ko_genomes": os.path.join(db_dir, "all_genomes_kos_genomes.txt"),
        "ko_list": os.path.join(db_dir, "all_genomes_kos_kos.txt"),

        "method": method,
        "threads": threads,
        "read1_suffix": r1_suf,
        "read2_suffix": r2_suf,
        "metadata": metadata if metadata else "none",

        "sanity_check": {
            "warn_loss_pct": 80,
            "fail_loss_pct": 90,
        },

        "identity_threshold": 0.4,
        "min_coverage": 0.5,
        "ko_inference": {"top_n_hits": 10},

        "vsearch_taxonomy": {
            "identity": 0.97,
            "min_coverage": 0.90,
            "top_n": 5,
        },

        "cutadapt": {
            "primer_f": "MAYGARAARMGNATGYTNCARGA",
            "primer_r": "GMCATYTGRTCNCCRTCRAA",
            "error_rate": 0.1,
            "minimum_length": 100,
        },

        "dada2": {
            "amplicon_len": 521,
            "primer_len_f": 23,
            "primer_len_r": 20,
            "read_length": 300,
            "truncLen_F": 260,
            "truncLen_R": 240,
            "maxEE_F": 2,
            "maxEE_R": 2,
            "pool": False,
            "min_overlap": 20,
            "min_boot": 50,
        },
    }

    return config, n_samples


def write_config(config, path="config.yaml"):
    """Write config dictionary to YAML file."""
    with open(path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
