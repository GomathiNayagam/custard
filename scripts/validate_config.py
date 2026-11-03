#!/usr/bin/env python3
"""
validate_config.py — Validate CUSTARD config.yaml at Snakemake DAG-build time.

Called from the top of the Snakefile before any rules are defined.
Exits with a clear error message if any problem is found — surfaces at
`snakemake --dry-run` rather than deep into an actual run.
"""

import os
import sys


def validate(config):
    errors   = []
    warnings = []

    # -------------------------------------------------------------------------
    # Required top-level keys
    # -------------------------------------------------------------------------
    required_keys = [
        "input_dir", "output_dir", "taxa_db", "kraken_db",
        "ko_db_fasta", "ko_npz", "ko_genomes", "ko_list",
        "threads", "read1_suffix", "read2_suffix",
        "identity_threshold", "min_coverage",
    ]
    for k in required_keys:
        if k not in config:
            errors.append(f"Missing required config key: '{k}'")

    if errors:
        _report(errors, warnings)   # can't continue without required keys

    # -------------------------------------------------------------------------
    # Path existence checks
    # -------------------------------------------------------------------------
    path_keys = [
        ("input_dir",   True,  "directory"),
        ("taxa_db",     False, "file"),
        ("ko_db_fasta", False, "file"),
        ("ko_npz",      False, "file"),
        ("ko_genomes",  False, "file"),
        ("ko_list",     False, "file"),
    ]
    for key, is_dir, kind in path_keys:
        p = config[key]
        if is_dir:
            if not os.path.isdir(p):
                errors.append(f"config['{key}'] {kind} not found: {p}")
        else:
            if not os.path.isfile(p):
                errors.append(f"config['{key}'] {kind} not found: {p}")

    # Kraken DB — check for hash.k2d sentinel of a complete build
    kraken_db = config["kraken_db"]
    if not os.path.isdir(kraken_db):
        errors.append(f"config['kraken_db'] directory not found: {kraken_db}")
    elif not os.path.exists(os.path.join(kraken_db, "hash.k2d")):
        errors.append(
            f"config['kraken_db'] exists but appears incomplete "
            f"(hash.k2d missing): {kraken_db}"
        )

    # Metadata (only checked when not "none")
    metadata = config.get("metadata", "none")
    if metadata != "none" and not os.path.isfile(metadata):
        errors.append(f"config['metadata'] file not found: {metadata}")

    # -------------------------------------------------------------------------
    # Numeric range checks
    # -------------------------------------------------------------------------
    for key, lo, hi in [
        ("identity_threshold", 0.0, 1.0),
        ("min_coverage",       0.0, 1.0),
        ("threads",            1,   512),
    ]:
        val = config.get(key)
        if val is not None:
            try:
                if not (lo <= float(val) <= hi):
                    errors.append(
                        f"config['{key}'] = {val} out of range [{lo}, {hi}]"
                    )
            except (TypeError, ValueError):
                errors.append(f"config['{key}'] must be numeric, got: {val!r}")

    # -------------------------------------------------------------------------
    # DADA2 sub-section
    # -------------------------------------------------------------------------
    dada2 = config.get("dada2", {})
    if not dada2:
        errors.append("Missing 'dada2:' section in config")
    else:
        for k in ["amplicon_len", "primer_len_f", "primer_len_r",
                  "read_length", "maxEE_F", "maxEE_R", "min_overlap"]:
            if k not in dada2:
                errors.append(f"Missing dada2.{k}")

        # pool value — must be unquoted in YAML
        pool = dada2.get("pool")
        if pool not in ("pseudo", True, False):
            errors.append(
                f"dada2.pool must be  pseudo  false  or  true  (unquoted). "
                f"Got: {pool!r}. Quoted 'false'/'true' causes an R error."
            )

        # truncLen — both or neither
        has_f = "truncLen_F" in dada2 and dada2["truncLen_F"] is not None
        has_r = "truncLen_R" in dada2 and dada2["truncLen_R"] is not None
        if has_f != has_r:
            errors.append(
                "dada2.truncLen_F and dada2.truncLen_R must both be set "
                "or both commented out — setting only one is not allowed."
            )
        if has_f and has_r:
            tF     = int(dada2["truncLen_F"])
            tR     = int(dada2["truncLen_R"])
            insert = (dada2.get("amplicon_len", 521)
                      - dada2.get("primer_len_f", 0)
                      - dada2.get("primer_len_r", 0))
            overlap = tF + tR - insert
            if overlap < 12:
                errors.append(
                    f"Merge overlap = {overlap} bp "
                    f"(truncLen_F={tF} + truncLen_R={tR} - insert={insert}). "
                    f"Minimum safe overlap is 12 bp — increase truncLen values."
                )
            elif overlap < 20:
                warnings.append(
                    f"Merge overlap = {overlap} bp (truncLen_F={tF}, "
                    f"truncLen_R={tR}). Recommend >= 20 bp for reliable merging."
                )

        if float(dada2.get("maxEE_R", 2)) > 5:
            warnings.append(
                f"dada2.maxEE_R = {dada2['maxEE_R']} is unusually permissive. "
                f"Values > 5 may pass too many low-quality reverse reads."
            )

    # -------------------------------------------------------------------------
    # Cutadapt sub-section
    # -------------------------------------------------------------------------
    cutadapt = config.get("cutadapt", {})
    if not cutadapt:
        errors.append("Missing 'cutadapt:' section in config")
    else:
        for k in ["primer_f", "primer_r", "error_rate", "minimum_length"]:
            if k not in cutadapt:
                errors.append(f"Missing cutadapt.{k}")
        if float(cutadapt.get("error_rate", 0.1)) > 0.3:
            warnings.append(
                f"cutadapt.error_rate = {cutadapt['error_rate']} is high. "
                f"Typical range: 0.1–0.15."
            )

    # -------------------------------------------------------------------------
    # KO inference sub-section
    # -------------------------------------------------------------------------
    ko_inf = config.get("ko_inference", {})
    if not ko_inf:
        errors.append(
            "Missing 'ko_inference:' section. "
            "Add:  ko_inference:\\n  top_n_hits: 10"
        )
    else:
        top_n = ko_inf.get("top_n_hits")
        if top_n is None:
            errors.append("Missing ko_inference.top_n_hits")
        elif int(top_n) < 1:
            errors.append(f"ko_inference.top_n_hits must be >= 1, got {top_n}")
        elif int(top_n) == 1:
            warnings.append(
                "ko_inference.top_n_hits = 1 (single-best-hit mode). "
                "Consider >= 5 for more robust KO predictions."
            )

    # -------------------------------------------------------------------------
    # Vsearch taxonomy sub-section
    # -------------------------------------------------------------------------
    vt = config.get("vsearch_taxonomy", {})
    if not vt:
        errors.append(
            "Missing 'vsearch_taxonomy:' section. "
            "Add:\n  vsearch_taxonomy:\n    identity: 0.97\n    min_coverage: 0.90\n    top_n: 5"
        )
    else:
        idt = vt.get("identity", 0.97)
        cov = vt.get("min_coverage", 0.90)
        try:
            if not (0.5 <= float(idt) <= 1.0):
                errors.append(f"vsearch_taxonomy.identity = {idt} out of range [0.5, 1.0]")
            if not (0.5 <= float(cov) <= 1.0):
                errors.append(f"vsearch_taxonomy.min_coverage = {cov} out of range [0.5, 1.0]")
        except (TypeError, ValueError):
            errors.append("vsearch_taxonomy.identity and min_coverage must be numeric")
        if float(idt) < 0.9:
            warnings.append(
                f"vsearch_taxonomy.identity = {idt} is below 90%. "
                f"Species-level calls typically require >= 97%."
            )

    # -------------------------------------------------------------------------
    # Sanity check thresholds (optional — defaults used if absent)
    # -------------------------------------------------------------------------
    sc = config.get("sanity_check", {})
    warn_pct = float(sc.get("warn_loss_pct", 80))
    fail_pct = float(sc.get("fail_loss_pct", 90))
    if warn_pct >= fail_pct:
        errors.append(
            f"sanity_check.warn_loss_pct ({warn_pct}) must be < "
            f"fail_loss_pct ({fail_pct})"
        )

    _report(errors, warnings)


def _report(errors, warnings):
    bar = "═" * 80
    print("", flush=True)
    if warnings:
        print(f"╔══ CUSTARD config warnings ══{'═'*51}", flush=True)
        for w in warnings:
            print(f"║  ⚠  {w}", flush=True)
        print(f"╚{bar}", flush=True)
    if errors:
        print(f"╔══ CUSTARD config ERRORS — fix before running ══{'═'*32}", flush=True)
        for e in errors:
            print(f"║  ✗  {e}", flush=True)
        print(f"╚{bar}", flush=True)
        print("Run `snakemake --dry-run` after fixing to verify.\n", flush=True)
        sys.exit(1)
    elif warnings:
        print("[validate_config] Warnings above — pipeline will still run.\n", flush=True)
    else:
        print("[validate_config] Config OK.\n", flush=True)


if __name__ == "__main__":
    import yaml
    if len(sys.argv) < 2:
        sys.exit("Usage: validate_config.py <config.yaml>")
    with open(sys.argv[1]) as f:
        cfg = yaml.safe_load(f)
    validate(cfg)
