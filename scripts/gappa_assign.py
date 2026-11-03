#!/usr/bin/env python3
"""
scripts/gappa_assign.py
=======================
Called by Snakemake rule gappa_assign.

Steps:
  1. Parse taxonomy from DADA2-formatted FASTA headers into GAPPA taxonomy file.
  2. Run gappa examine assign on the EPA-ng .jplace file.
  3. Normalise GAPPA's per_query output into a clean TSV with one row per ASV.

Input header format (from dada2_formatted_rpoc_amplicon.fasta):
  >GCF_016756535.1;taxonomy=d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;
   o__Rhizobiales;f__Rhizobiaceae;g__Mesorhizobium;s__Mesorhizobium sp016756535

GAPPA taxonomy file format (two columns, tab-separated, no header):
  GCF_016756535.1<TAB>d__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;
                      o__Rhizobiales;f__Rhizobiaceae;g__Mesorhizobium;
                      s__Mesorhizobium sp016756535

Output TSV columns:
  asv_id | placement_taxonomy | placement_lwr | placement_rank | placement_confident
"""

import os
import re
import sys
import json
import subprocess
import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake bindings
# ---------------------------------------------------------------------------
jplace_path   = snakemake.input["jplace"]
db_fasta_path = snakemake.input["db_fasta"]
out_tax_path  = snakemake.output["placement_tax"]
ref_tax_path  = snakemake.output["ref_tax_file"]
log_path      = snakemake.output["log"]
out_dir       = snakemake.params["out_dir"]
min_lwr       = float(snakemake.params["min_lwr"])
num_best_pl   = int(snakemake.params["num_best_pl"])
threads       = snakemake.threads

os.makedirs(out_dir, exist_ok=True)

def log(msg, fh=None):
    print(f"[gappa_assign] {msg}", flush=True)
    if fh:
        fh.write(f"[gappa_assign] {msg}\n")
        fh.flush()

# ---------------------------------------------------------------------------
# Step 1: Parse DADA2 FASTA headers → GAPPA taxonomy file
# ---------------------------------------------------------------------------
with open(log_path, "w") as logf:
    log("Parsing taxonomy from DADA2 FASTA headers...", logf)

    written = 0
    skipped = 0
    with open(db_fasta_path) as fh, open(ref_tax_path, "w") as out:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            # Match: accession;taxonomy=d__X;p__Y;...
            m = re.match(r'^([^;]+);taxonomy=(.+)$', header)
            if m:
                acc     = m.group(1).strip()
                lineage = m.group(2).strip()
                # GAPPA expects semicolon-delimited lineage: d__X;p__Y;c__Z;...
                # Already in this format from the DADA2 header — write directly
                out.write(f"{acc}\t{lineage}\n")
                written += 1
            else:
                skipped += 1

    log(f"  Taxonomy entries written: {written}", logf)
    if skipped:
        log(f"  Headers skipped (unexpected format): {skipped}", logf)
    if written == 0:
        log("ERROR: No taxonomy entries parsed. Check db_fasta header format.", logf)
        sys.exit(1)

    # ---------------------------------------------------------------------------
    # Step 2: Run GAPPA examine assign
    # ---------------------------------------------------------------------------
    log("Running GAPPA examine assign...", logf)

    gappa_cmd = [
        "gappa", "examine", "assign",
        "--jplace-path",       jplace_path,
        "--taxon-file",        ref_tax_path,
        "--best-n-taxa",       str(num_best_pl),
        "--min-lwr",           str(min_lwr),
        "--per-query-results",
        "--out-dir",           out_dir,
        "--threads",           str(threads),
        "--redo",
    ]

    log(f"  Command: {' '.join(gappa_cmd)}", logf)
    result = subprocess.run(
        gappa_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    logf.write(result.stdout)
    if result.returncode != 0:
        log(f"ERROR: GAPPA exited with code {result.returncode}", logf)
        sys.exit(1)

    log("GAPPA completed.", logf)

    # ---------------------------------------------------------------------------
    # Step 3: Find GAPPA's per-query output file
    # GAPPA writes: per_query_assign.tsv  (or assign_per_query.tsv in some versions)
    # ---------------------------------------------------------------------------
    candidates = [
        os.path.join(out_dir, "per_query_assign.tsv"),
        os.path.join(out_dir, "assign_per_query.tsv"),
    ]
    gappa_out = None
    for c in candidates:
        if os.path.exists(c):
            gappa_out = c
            break
    if gappa_out is None:
        # Fall back: find any TSV in out_dir that isn't our ref_tax_file
        tsv_files = [
            os.path.join(out_dir, f)
            for f in os.listdir(out_dir)
            if f.endswith(".tsv")
            and os.path.join(out_dir, f) not in (ref_tax_path, out_tax_path)
        ]
        if tsv_files:
            gappa_out = sorted(tsv_files)[-1]
            log(f"  Using fallback GAPPA output: {gappa_out}", logf)
        else:
            log("ERROR: Cannot find GAPPA per-query output file.", logf)
            sys.exit(1)

    # ---------------------------------------------------------------------------
    # Step 4: Normalise GAPPA output → clean TSV
    # ---------------------------------------------------------------------------
    # GAPPA per-query TSV typically has columns:
    #   name | LWR | fract | entropy | taxonomy
    # Column names vary by GAPPA version — we detect them dynamically.

    log(f"Normalising GAPPA output: {gappa_out}", logf)

    gappa_df = pd.read_csv(gappa_out, sep="\t", dtype=str)
    gappa_df.columns = [c.strip().lower() for c in gappa_df.columns]

    log(f"  GAPPA output columns: {list(gappa_df.columns)}", logf)

    # Identify key columns robustly
    name_col = next(
        (c for c in gappa_df.columns if c in ("name", "query", "sequence", "asv")),
        gappa_df.columns[0]
    )
    lwr_col = next(
        (c for c in gappa_df.columns if "lwr" in c or "weight" in c),
        None
    )
    tax_col = next(
        (c for c in gappa_df.columns
         if "tax" in c or "lineage" in c or "assign" in c),
        gappa_df.columns[-1]
    )

    log(f"  Using: name='{name_col}', lwr='{lwr_col}', taxonomy='{tax_col}'", logf)

    # Build normalised output
    rows = []
    for _, row in gappa_df.iterrows():
        asv_id   = str(row[name_col]).strip()
        taxonomy = str(row[tax_col]).strip() if pd.notna(row[tax_col]) else ""
        lwr      = float(row[lwr_col]) if lwr_col and pd.notna(row[lwr_col]) else 0.0
        confident = lwr >= min_lwr

        # Parse the deepest non-empty rank for the rank_assigned column
        # Lineage format: d__Bacteria;p__Pseudomonadota;...
        rank_assigned = "unclassified"
        if taxonomy and taxonomy not in ("", "nan", "unclassified"):
            parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
            if parts:
                rank_assigned = parts[-1]   # deepest assigned rank

        rows.append({
            "asv_id"              : asv_id,
            "placement_taxonomy"  : taxonomy if taxonomy != "nan" else "",
            "placement_lwr"       : round(lwr, 4),
            "placement_rank"      : rank_assigned,
            "placement_confident" : confident,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_tax_path, sep="\t", index=False)

    n_confident   = out_df["placement_confident"].sum()
    n_total       = len(out_df)
    n_unclassified = (out_df["placement_taxonomy"] == "").sum()

    log(f"  Total ASVs processed    : {n_total}", logf)
    log(f"  Confident (LWR>={min_lwr}) : {n_confident} ({n_confident/max(n_total,1):.1%})", logf)
    log(f"  Unclassified            : {n_unclassified}", logf)
    log(f"  Output -> {out_tax_path}", logf)
