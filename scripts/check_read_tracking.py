#!/usr/bin/env python3
"""
check_read_tracking.py — Post-DADA2 sanity check on read retention per sample.

Reads read_tracking.tsv produced by run_dada2.R and checks for:
  - Any sample losing more than warn_loss_pct % of reads at any single step
  - Any sample with < min_final_pct % reads surviving to final merged output
  - All samples failing (run-level failure — likely a config issue)

Exits with code 1 (fail) or 2 (warn-only) so Snakemake can be configured
to treat warnings as non-fatal.

read_tracking.tsv expected columns (DADA2 standard):
  sample | input | filtered | denoised_F | denoised_R | merged | nonchimeric
"""

import sys
import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tracking",       required=True,
                   help="Path to read_tracking.tsv")
    p.add_argument("--output",         required=True,
                   help="Path to write sanity check report TSV")
    p.add_argument("--warn_loss_pct",  type=float, default=80,
                   help="Warn if any step loses more than this %% of reads (default: 80)")
    p.add_argument("--fail_loss_pct",  type=float, default=90,
                   help="Fail if any step loses more than this %% of reads (default: 90)")
    p.add_argument("--min_final_pct",  type=float, default=10,
                   help="Fail if final retention is below this %% for ALL samples (default: 10)")
    return p.parse_args()


STEP_PAIRS = [
    ("input",      "filtered",    "filterAndTrim"),
    ("filtered",   "merged",      "denoise+merge"),
    ("merged",     "nonchimeric", "chimera_removal"),
]


def check(df, warn_loss, fail_loss, min_final):
    issues    = []   # (sample, step, pct_lost, severity)
    per_sample = []

    for _, row in df.iterrows():
        sample = row.get("sample", row.name)
        input_reads = row.get("input", 0)
        if input_reads == 0:
            issues.append((sample, "input", 100.0, "FAIL",
                           "Sample has 0 input reads"))
            continue

        for col_in, col_out, step_name in STEP_PAIRS:
            if col_in not in row or col_out not in row:
                continue
            n_in  = row[col_in]
            n_out = row[col_out]
            if n_in == 0:
                continue
            pct_lost = 100.0 * (1 - n_out / n_in)
            if pct_lost >= fail_loss:
                issues.append((sample, step_name, pct_lost, "FAIL",
                               f"Lost {pct_lost:.1f}% at {step_name} "
                               f"({n_in:,} -> {n_out:,} reads)"))
            elif pct_lost >= warn_loss:
                issues.append((sample, step_name, pct_lost, "WARN",
                               f"Lost {pct_lost:.1f}% at {step_name} "
                               f"({n_in:,} -> {n_out:,} reads)"))

        # Overall retention
        final     = row.get("nonchimeric", row.get("merged", 0))
        final_pct = 100.0 * final / input_reads
        per_sample.append({
            "sample":      sample,
            "input":       int(input_reads),
            "final":       int(final),
            "retention_%": round(final_pct, 1),
        })

    return issues, per_sample


def main():
    args = parse_args()

    try:
        df = pd.read_csv(args.tracking, sep="\t")
    except Exception as e:
        sys.exit(f"ERROR: Could not read read_tracking.tsv: {e}")

    # Normalise column names — DADA2 sometimes writes 'Sample' capitalised
    df.columns = [c.lower() for c in df.columns]
    if "sample" not in df.columns:
        df = df.reset_index()
        df = df.rename(columns={"index": "sample"})

    issues, per_sample = check(
        df,
        warn_loss  = args.warn_loss_pct,
        fail_loss  = args.fail_loss_pct,
        min_final  = args.min_final_pct,
    )

    # Write per-sample retention table
    report_df = pd.DataFrame(per_sample)
    report_df.to_csv(args.output, sep="\t", index=False)

    # Print issues
    has_fail = any(i[3] == "FAIL" for i in issues)
    has_warn = any(i[3] == "WARN" for i in issues)

    print("\n=== Read Tracking Sanity Check ===")
    if not issues:
        print("  All samples passed retention thresholds.")
    else:
        for sample, step, pct, severity, msg in issues:
            icon = "✗" if severity == "FAIL" else "⚠"
            print(f"  {icon} [{sample}] {msg}")

    # Summary table
    print(f"\n  {'Sample':<30} {'Input':>10} {'Final':>10} {'Retention':>12}")
    print(f"  {'-'*30} {'-'*10} {'-'*10} {'-'*12}")
    for r in per_sample:
        flag = ""
        if r["retention_%"] < args.min_final_pct:
            flag = " ✗"
        elif r["retention_%"] < (100 - args.warn_loss_pct):
            flag = " ⚠"
        print(f"  {r['sample']:<30} {r['input']:>10,} {r['final']:>10,} "
              f"{r['retention_%']:>11.1f}%{flag}")

    # Check global failure: ALL samples below min_final_pct
    final_pcts = [r["retention_%"] for r in per_sample]
    all_failing = all(p < args.min_final_pct for p in final_pcts) and len(final_pcts) > 0

    print("")
    if all_failing:
        print("FATAL: Every sample is below the minimum final retention threshold.")
        print("  Likely causes:")
        print("  1. Wrong truncLen — reads too short to merge after truncation.")
        print("     Check qc_plots.pdf and set truncLen_F/truncLen_R in config.yaml.")
        print("  2. Wrong primers — cutadapt is discarding all reads.")
        print("     Check trimmed/ FASTQ sizes and cutadapt log files.")
        print("  3. Wrong amplicon_len — insert length calculation is off.")
        print("  Aborting. Fix config and rerun.")
        sys.exit(1)
    elif has_fail:
        print("FAIL: One or more samples exceeded the loss threshold.")
        print("      Check qc_plots.pdf and the samples listed above.")
        print("      The pipeline will continue but results may be unreliable.")
        sys.exit(1)
    elif has_warn:
        print("WARN: One or more samples show high read loss at a step.")
        print("      Pipeline continues. Review qc_plots.pdf.")
        sys.exit(2)
    else:
        print("Read tracking: OK")


if __name__ == "__main__":
    main()
