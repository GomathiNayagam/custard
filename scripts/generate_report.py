#!/usr/bin/env python3
"""
generate_report.py — Self-contained HTML summary report for a CUSTARD run.

Always-on sections (no metadata required):
  1. Run summary
  2. Read tracking table (colour-coded retention)
  3. ASV length distribution
  4. Top-20 genera bar chart (all samples combined)
  5. Discarded reads taxonomy
  6. QC plots (embedded from qc_plots.pdf pages as images)
  7. Software versions
"""

import argparse
import base64
import os
import subprocess
import sys
from datetime import datetime
from io import BytesIO

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--asv_table",          required=True)
    p.add_argument("--taxonomy",           required=True)
    p.add_argument("--read_tracking",      required=True)
    p.add_argument("--read_tracking_qc",   required=True)
    p.add_argument("--asv_fasta",          required=True)
    p.add_argument("--discarded_taxonomy", default=None)
    p.add_argument("--ko_table",           default=None)
    p.add_argument("--qc_plots_pdf",       default=None,
                   help="Path to qc_plots.pdf from DADA2")
    p.add_argument("--config_yaml",        required=True)
    p.add_argument("--output",             required=True)
    p.add_argument("--metadata",           default=None)
    p.add_argument("--pipeline_version",   default="1.0.0")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def img_tag(b64, alt="plot", width="100%"):
    return (f'<img src="data:image/png;base64,{b64}" alt="{alt}" '
            f'style="max-width:{width};border-radius:6px;">')


def section(title, content, icon=""):
    return f"""
<section>
  <h2>{icon} {title}</h2>
  <div class="section-body">{content}</div>
</section>"""


def placeholder(msg):
    return f'<p class="placeholder">⬜ {msg}</p>'


def load_asv_table(path):
    """Load ASV table; return DataFrame or None. ASV IDs are always strings."""
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
        # Normalise ASV column name
        asv_col = next((c for c in df.columns
                        if c.upper() in ("ASV", "ASVID", "ASV_ID", "#OTU ID")),
                       df.columns[0])
        df = df.rename(columns={asv_col: "ASV"})
        return df
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------
def build_run_summary(args, asv_table, config_yaml):
    try:
        with open(config_yaml) as f:
            cfg_text = f.read()
    except Exception:
        cfg_text = "(config.yaml not readable)"

    if asv_table is not None:
        sample_cols = [c for c in asv_table.columns if c != "ASV"]
        n_samples   = len(sample_cols)
        n_asvs      = len(asv_table)   # one row = one ASV
    else:
        n_samples = n_asvs = "?"

    run_date = datetime.now().strftime("%Y-%m-%d %H:%M")
    html = f"""
<table class="summary-table">
  <tr><td>Date</td><td>{run_date}</td></tr>
  <tr><td>Pipeline version</td><td>{args.pipeline_version}</td></tr>
  <tr><td>Samples</td><td>{n_samples}</td></tr>
  <tr><td>Total ASVs</td><td>{n_asvs:,}</td></tr>
</table>
<details><summary>Config snapshot</summary>
<pre class="code-block">{cfg_text}</pre>
</details>"""
    return html


def build_read_tracking(tracking_path, qc_path):
    try:
        df = pd.read_csv(tracking_path, sep="\t")
        df.columns = [c.lower() for c in df.columns]
    except Exception as e:
        return f"<p>Could not load read_tracking.tsv: {e}</p>"

    if "sample" not in df.columns:
        df = df.reset_index().rename(columns={"index": "sample"})

    step_cols = [c for c in ["input", "filtered", "denoised_f", "denoised_r",
                              "merged", "nonchimeric"] if c in df.columns]

    rows_html = ""
    for _, row in df.iterrows():
        input_n = int(row.get("input", 1) or 1)
        cells   = f"<td>{row.get('sample', '?')}</td>"
        for col in step_cols:
            val = row.get(col, "")
            try:
                n   = int(val)
                pct = 100.0 * n / input_n
                cls = "good" if pct >= 70 else ("warn" if pct >= 40 else "bad")
                cells += (f'<td class="{cls}">{n:,}'
                          f'<br><small>{pct:.0f}%</small></td>')
            except Exception:
                cells += f"<td>{val}</td>"
        rows_html += f"<tr>{cells}</tr>"

    header = ("<tr><th>Sample</th>"
              + "".join(f"<th>{c}</th>" for c in step_cols)
              + "</tr>")
    table  = f'<table class="data-table">{header}{rows_html}</table>'

    qc_note = ""
    try:
        qc_df = pd.read_csv(qc_path, sep="\t")
        low   = qc_df[qc_df["retention_%"] < 40]
        if not low.empty:
            samples = ", ".join(low["sample"].astype(str).tolist())
            qc_note = (f'<p class="warn-box">⚠ Low retention (&lt;40%) in: '
                       f'<strong>{samples}</strong>. '
                       f'Check qc_plots.pdf and consider adjusting truncLen.</p>')
    except Exception:
        pass

    return qc_note + table


def build_asv_length_dist(asv_fasta_path):
    """Parse ASV FASTA and plot length distribution.
    Reads the full multi-line FASTA correctly — sequence lines only,
    no header lines, no blank lines.
    """
    try:
        lengths  = []
        seq_buf  = []
        with open(asv_fasta_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if seq_buf:
                        lengths.append(len("".join(seq_buf)))
                        seq_buf = []
                else:
                    seq_buf.append(line)
            if seq_buf:                          # last sequence
                lengths.append(len("".join(seq_buf)))

        if not lengths:
            return placeholder("No ASV sequences found in FASTA")

        unique, counts = np.unique(lengths, return_counts=True)
        modal_len      = unique[np.argmax(counts)]

        fig, ax = plt.subplots(figsize=(7, 3.5))
        ax.bar(unique, counts, color="#4C72B0", edgecolor="white", linewidth=0.5)
        ax.set_xlabel("ASV length (bp)")
        ax.set_ylabel("Number of ASVs")
        ax.set_title(
            f"ASV length distribution  (n={len(lengths):,} ASVs | "
            f"modal length = {modal_len} bp)"
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        return img_tag(fig_to_b64(fig), "ASV length distribution")
    except Exception as e:
        return placeholder(f"Could not generate length plot: {e}")


def build_taxonomy_barchart(taxonomy_path, asv_table, top_n=20):
    """Top-N genera bar chart, read-count weighted."""
    try:
        tax = pd.read_csv(taxonomy_path, sep="\t", dtype=str)
        tax.columns = [c.lower() for c in tax.columns]

        # Identify ASV column (first column or one named 'asv*')
        asv_col = next(
            (c for c in tax.columns if c.startswith("asv")),
            tax.columns[0]
        )
        # Identify genus column
        genus_col = next(
            (c for c in tax.columns if "genus" in c),
            None
        )
        if genus_col is None:
            return placeholder("No genus column found in taxonomy table")

        # Build per-ASV total read count from asv_table
        if asv_table is not None:
            sample_cols = [c for c in asv_table.columns if c != "ASV"]
            # Coerce count columns to numeric (they were loaded as str)
            counts_df   = asv_table[sample_cols].apply(pd.to_numeric, errors="coerce")
            totals      = counts_df.sum(axis=1)
            totals_df   = pd.DataFrame({"ASV": asv_table["ASV"].astype(str),
                                        "total": totals})
            # Merge on ASV id (both now str)
            tax[asv_col] = tax[asv_col].astype(str)
            merged = tax.merge(
                totals_df.rename(columns={"ASV": asv_col}),
                on=asv_col, how="left"
            )
        else:
            merged        = tax.copy()
            merged["total"] = 1

        merged["total"] = pd.to_numeric(merged["total"], errors="coerce").fillna(0)

        # Replace empty/unclassified genus with "Unclassified"
        merged[genus_col] = merged[genus_col].fillna("Unclassified").replace("", "Unclassified")

        genus_counts = (
            merged.groupby(genus_col)["total"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
        )

        if genus_counts.empty:
            return placeholder("No genus-level classifications found")

        fig, ax = plt.subplots(figsize=(9, 4.5))
        colors  = plt.cm.tab20.colors
        ax.barh(
            range(len(genus_counts)), genus_counts.values,
            color=[colors[i % 20] for i in range(len(genus_counts))]
        )
        ax.set_yticks(range(len(genus_counts)))
        ax.set_yticklabels(genus_counts.index, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Total reads")
        ax.set_title(f"Top {top_n} genera — all samples combined")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        return img_tag(fig_to_b64(fig), "Top genera")
    except Exception as e:
        return placeholder(f"Could not generate taxonomy chart: {e}")


def build_discarded_taxonomy(discarded_tax_path):
    if not discarded_tax_path or not os.path.exists(discarded_tax_path):
        return placeholder("Discarded reads taxonomy not available")
    try:
        df = pd.read_csv(discarded_tax_path, sep="\t", dtype=str)
        df.columns = [c.lower() for c in df.columns]
        genus_col = next((c for c in df.columns if "genus" in c), None)
        if genus_col is None:
            return placeholder("No genus column in discarded reads taxonomy")

        df[genus_col] = df[genus_col].fillna("Unclassified").replace("", "Unclassified")
        counts        = df[genus_col].value_counts().head(15)

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.barh(range(len(counts)), counts.values, color="#DD8452")
        ax.set_yticks(range(len(counts)))
        ax.set_yticklabels(counts.index, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Read count")
        ax.set_title("Top 15 genera in filterAndTrim-discarded reads")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()

        n_total      = len(df)
        n_classified = (df[genus_col] != "Unclassified").sum()
        note = (
            f"<p><strong>{n_total:,}</strong> discarded reads classified | "
            f"<strong>{n_classified:,}</strong> "
            f"({100*n_classified/max(n_total,1):.0f}%) assigned to genus level.</p>"
            f"<p class='note'>These reads failed DADA2 quality filtering. "
            f"If dominant taxa here match your expected community, consider "
            f"relaxing maxEE_R or adjusting truncLen.</p>"
        )
        return note + img_tag(fig_to_b64(fig), "Discarded reads taxonomy")
    except Exception as e:
        return placeholder(f"Could not generate discarded reads chart: {e}")


def build_qc_plots(qc_plots_pdf):
    """Convert qc_plots.pdf pages to PNG and embed inline."""
    if not qc_plots_pdf or not os.path.exists(qc_plots_pdf):
        return placeholder(
            "qc_plots.pdf not found. It is produced by DADA2 after a full run."
        )
    try:
        # Use pdftoppm (poppler) to convert PDF pages to PNG
        result = subprocess.run(
            ["pdftoppm", "-r", "150", "-png", qc_plots_pdf, "/tmp/custard_qc"],
            capture_output=True
        )
        import glob
        pages = sorted(glob.glob("/tmp/custard_qc-*.png"))
        if not pages:
            # Fallback: try convert (ImageMagick)
            subprocess.run(
                ["convert", "-density", "150", qc_plots_pdf,
                 "/tmp/custard_qc_im-%d.png"],
                capture_output=True
            )
            pages = sorted(glob.glob("/tmp/custard_qc_im-*.png"))

        if not pages:
            return placeholder(
                "Could not convert qc_plots.pdf to images. "
                "Install poppler-utils (pdftoppm) or ImageMagick."
            )

        imgs = []
        for page in pages:
            with open(page, "rb") as f:
                b64 = base64.b64encode(f.read()).decode()
            imgs.append(img_tag(b64, "QC plot", width="48%"))

        # Show 2 per row
        rows = []
        for i in range(0, len(imgs), 2):
            pair = imgs[i:i+2]
            rows.append(
                f'<div style="display:flex;gap:12px;margin-bottom:12px;">'
                + "".join(pair) + "</div>"
            )
        return "".join(rows)

    except FileNotFoundError:
        return placeholder(
            "qc_plots.pdf found but PDF-to-image conversion tools not available. "
            "Install poppler-utils: conda install -c conda-forge poppler"
        )
    except Exception as e:
        return placeholder(f"Could not render qc_plots.pdf: {e}")


def build_software_versions():
    tools = {
        "python":    ["python3",   "--version"],
        "snakemake": ["snakemake", "--version"],
        "cutadapt":  ["cutadapt",  "--version"],
        "kraken2":   ["kraken2",   "--version"],
        "vsearch":   ["vsearch",   "--version"],
    }
    rows = ""
    for name, cmd in tools.items():
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode().strip()
            ver = out.split("\n")[0]
        except Exception:
            ver = "not found in PATH"
        rows += f"<tr><td>{name}</td><td><code>{ver}</code></td></tr>"

    return f"""
<table class="summary-table">
  <tr><th>Tool</th><th>Version</th></tr>
  {rows}
</table>
<p class="note">R package versions are recorded in qc_plots.pdf session info.</p>"""


# ---------------------------------------------------------------------------
# CSS
# ---------------------------------------------------------------------------
CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       background: #f5f6fa; color: #222; line-height: 1.6; }
header { background: #1a1a2e; color: white; padding: 24px 40px; }
header h1 { font-size: 1.6rem; font-weight: 700; letter-spacing: 0.5px; }
header p  { font-size: 0.9rem; opacity: 0.7; margin-top: 4px; }
main { max-width: 1100px; margin: 32px auto; padding: 0 24px 48px; }
section { background: white; border-radius: 10px; padding: 28px 32px;
          margin-bottom: 24px; box-shadow: 0 1px 4px rgba(0,0,0,0.08); }
h2 { font-size: 1.15rem; font-weight: 600; color: #1a1a2e;
     border-bottom: 2px solid #f0f0f0; padding-bottom: 10px; margin-bottom: 18px; }
.section-body { font-size: 0.92rem; }
table.summary-table { border-collapse: collapse; width: 100%; }
table.summary-table td, table.summary-table th
  { padding: 8px 14px; border-bottom: 1px solid #eee; }
table.summary-table th { text-align: left; background: #f8f9fa; }
table.summary-table tr:last-child td { border-bottom: none; }
table.data-table { border-collapse: collapse; width: 100%; font-size: 0.84rem;
                   overflow-x: auto; display: block; }
table.data-table th { background: #1a1a2e; color: white;
                       padding: 8px 10px; text-align: left; white-space: nowrap; }
table.data-table td { padding: 6px 10px; border-bottom: 1px solid #f0f0f0;
                       text-align: right; }
table.data-table td:first-child { text-align: left; }
td.good { background: #d4edda; }
td.warn { background: #fff3cd; }
td.bad  { background: #f8d7da; font-weight: 600; }
.warn-box { background: #fff3cd; border-left: 4px solid #ffc107;
             padding: 10px 14px; border-radius: 4px; margin-bottom: 14px; }
.placeholder { color: #999; font-style: italic; padding: 16px 0; }
.note { color: #666; font-size: 0.85rem; margin-top: 8px; }
pre.code-block { background: #f8f9fa; border: 1px solid #e0e0e0;
                 border-radius: 6px; padding: 14px; font-size: 0.78rem;
                 overflow-x: auto; white-space: pre-wrap; margin-top: 10px; }
details summary { cursor: pointer; color: #4C72B0; font-weight: 500;
                  font-size: 0.9rem; margin-top: 10px; user-select: none; }
details summary:hover { text-decoration: underline; }
"""


def build_html(sections_html, run_date):
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CUSTARD Run Report</title>
  <style>{CSS}</style>
</head>
<body>
<header>
  <h1>🧫 CUSTARD Pipeline Report</h1>
  <p>rpoC amplicon sequencing · Bacterial community profiling · {run_date}</p>
</header>
<main>
{"".join(sections_html)}
</main>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args     = parse_args()
    run_date = datetime.now().strftime("%Y-%m-%d %H:%M")

    asv_table = load_asv_table(args.asv_table)

    sections = [
        section("Run Summary",
                build_run_summary(args, asv_table, args.config_yaml), "📋"),
        section("Read Tracking",
                build_read_tracking(args.read_tracking, args.read_tracking_qc), "📊"),
        section("ASV Length Distribution",
                build_asv_length_dist(args.asv_fasta), "📏"),
        section("Top 20 Genera — All Samples",
                build_taxonomy_barchart(args.taxonomy, asv_table), "🔬"),
        section("Discarded Reads — Taxonomy",
                build_discarded_taxonomy(args.discarded_taxonomy), "🗑️"),
        section("DADA2 QC Plots",
                build_qc_plots(args.qc_plots_pdf), "📈"),
        section("Software Versions",
                build_software_versions(), "⚙️"),
    ]

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        f.write(build_html(sections, run_date))

    print(f"[generate_report] Report written: {args.output}", flush=True)


if __name__ == "__main__":
    main()
