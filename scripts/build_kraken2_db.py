#!/usr/bin/env python3
"""Build Kraken2 database from rpocDB v3.1 amplicons.

Runs during `custard setup` after the v3.1 sequence files are downloaded.
Requires: kraken2-build, internet access (for NCBI taxonomy download).
"""

import os
import sys
import csv
import subprocess
import tempfile


def build_kraken2_db(db_dir, meta_tsv=None, threads=4):
    """Build Kraken2 DB from rpocDB v3.1 vsearch FASTA.

    Args:
        db_dir: path to data/db/ where rpocdb_v3.1_vsearch.fasta lives
        meta_tsv: optional GTDB metadata TSV for taxid mapping
        threads: CPU threads for kraken2-build
    """
    vsearch_fasta = os.path.join(db_dir, "rpocdb_v3.1_vsearch.fasta")
    kraken_dir = os.path.join(db_dir, "kraken_amplicon")

    if os.path.exists(os.path.join(kraken_dir, "hash.k2d")):
        print("[kraken2] Database already exists — skipping build.")
        return True

    if not os.path.exists(vsearch_fasta):
        print(f"[kraken2] vsearch FASTA not found: {vsearch_fasta}", file=sys.stderr)
        return False

    os.makedirs(kraken_dir, exist_ok=True)

    # Step 1: Create Kraken2-formatted FASTA with taxids
    print("[kraken2] Mapping amplicon accessions to NCBI taxids...", flush=True)

    # Try to find GTDB metadata for taxid mapping
    if meta_tsv is None:
        candidates = [
            os.path.join(db_dir, "bac120_metadata_r232.tsv"),
            os.path.expanduser("~/rpocDB_database_development/bac120_metadata_r232.tsv"),
        ]
        for c in candidates:
            if os.path.exists(c):
                meta_tsv = c
                break

    if meta_tsv and os.path.exists(meta_tsv):
        acc_to_taxid = {}
        with open(meta_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                clean = row['accession'].replace('RS_', '').replace('GB_', '')
                taxid = row.get('ncbi_taxid', '')
                if taxid and taxid != 'none':
                    acc_to_taxid[clean] = taxid
        print(f"[kraken2] Loaded {len(acc_to_taxid):,} taxid mappings from metadata")
    else:
        # Fallback: use taxid 2 (Bacteria) for all sequences
        print("[kraken2] No metadata found — using taxid=2 (Bacteria) for all sequences")
        acc_to_taxid = None

    kraken_fasta = os.path.join(kraken_dir, "library.fasta")
    n_seqs = 0
    with open(vsearch_fasta) as fin, open(kraken_fasta, 'w') as fout:
        acc = None
        seq = []
        for line in fin:
            if line.startswith('>'):
                if acc and seq:
                    if acc_to_taxid is not None:
                        taxid = acc_to_taxid.get(acc, '2')
                    else:
                        taxid = '2'
                    fout.write(f">{acc}|kraken:taxid|{taxid}\n{''.join(seq)}\n")
                    n_seqs += 1
                header = line.strip()[1:]
                acc = header.split(';')[0]
                seq = []
            else:
                seq.append(line.strip())
        if acc and seq:
            if acc_to_taxid is not None:
                taxid = acc_to_taxid.get(acc, '2')
            else:
                taxid = '2'
            fout.write(f">{acc}|kraken:taxid|{taxid}\n{''.join(seq)}\n")
            n_seqs += 1

    print(f"[kraken2] Formatted {n_seqs:,} sequences", flush=True)

    # Step 2: Download taxonomy
    print("[kraken2] Downloading NCBI taxonomy...", flush=True)
    tax_dir = os.path.join(kraken_dir, "taxonomy")
    os.makedirs(tax_dir, exist_ok=True)

    tax_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tax_gz = os.path.join(tax_dir, "taxdump.tar.gz")

    result = subprocess.run(
        ["wget", "-q", tax_url, "-O", tax_gz],
        capture_output=True, text=True, timeout=300
    )
    if result.returncode != 0:
        print(f"[kraken2] Failed to download taxonomy: {result.stderr}", file=sys.stderr)
        return False

    subprocess.run(["tar", "-xzf", tax_gz, "-C", tax_dir], check=True)
    os.remove(tax_gz)

    # Step 3: Add library and build
    print("[kraken2] Adding sequences to library...", flush=True)
    subprocess.run(
        ["kraken2-build", "--add-to-library", kraken_fasta, "--db", kraken_dir],
        capture_output=True, text=True, check=True
    )

    print(f"[kraken2] Building database ({threads} threads)...", flush=True)
    result = subprocess.run(
        ["kraken2-build", "--build", "--db", kraken_dir, "--threads", str(threads)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"[kraken2] Build failed: {result.stderr}", file=sys.stderr)
        return False

    # Step 4: Clean up bulky files
    for d in ["taxonomy", "library"]:
        full = os.path.join(kraken_dir, d)
        if os.path.isdir(full):
            import shutil
            shutil.rmtree(full)
    for f in [kraken_fasta, os.path.join(kraken_dir, "unmapped.txt")]:
        if os.path.exists(f):
            os.remove(f)

    print(f"[kraken2] Database built: {kraken_dir}", flush=True)
    db_size = sum(
        os.path.getsize(os.path.join(kraken_dir, f))
        for f in os.listdir(kraken_dir) if f.endswith('.k2d')
    ) / 1e6
    print(f"[kraken2] Size: {db_size:.0f} MB", flush=True)
    return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Build Kraken2 DB from rpocDB v3")
    parser.add_argument("--db-dir", required=True, help="Path to data/db/")
    parser.add_argument("--metadata", default=None, help="GTDB metadata TSV")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    success = build_kraken2_db(args.db_dir, args.metadata, args.threads)
    sys.exit(0 if success else 1)
