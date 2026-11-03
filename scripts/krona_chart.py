#!/usr/bin/env python3
"""Generate Krona interactive HTML sunburst from taxonomy data."""

import pandas as pd
import subprocess
import os
import sys
import tempfile

def main(asv_file, tax_file, output_html):
    os.makedirs(os.path.dirname(output_html), exist_ok=True)

    asv = pd.read_csv(asv_file, sep="\t", index_col=0)
    tax = pd.read_csv(tax_file, sep="\t")

    tax_ranks = [c for c in tax.columns if c in
                 ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    tax = tax.set_index('ASV_ID')

    krona_files = []
    sample_names = []

    tmpdir = tempfile.mkdtemp()

    for sample in asv.columns:
        sample_name = sample.replace('_F_filt.fastq.gz', '').replace('_L001', '')
        krona_path = os.path.join(tmpdir, f"{sample_name}.krona.txt")

        with open(krona_path, 'w') as f:
            for asv_id in asv.index:
                reads = asv.loc[asv_id, sample]
                if reads <= 0:
                    continue
                if asv_id in tax.index:
                    lineage = []
                    for rank in tax_ranks:
                        val = str(tax.loc[asv_id, rank]) if rank in tax.columns else 'Unclassified'
                        if val in ('nan', 'Unclassified', '', 'NA'):
                            val = 'Unclassified'
                        lineage.append(val)
                    f.write(f"{int(reads)}\t" + "\t".join(lineage) + "\n")
                else:
                    f.write(f"{int(reads)}\tUnclassified\n")

        krona_files.append(krona_path)
        sample_names.append(sample_name)

    cmd = ["ktImportText", "-o", output_html]
    for kf, sn in zip(krona_files, sample_names):
        cmd.append(f"{kf},{sn}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Cleanup
    for kf in krona_files:
        os.remove(kf)
    os.rmdir(tmpdir)

    if result.returncode != 0:
        print(f"[krona] Error: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    print(f"[krona] Generated: {output_html}")
    print(f"[krona] Samples: {', '.join(sample_names)}")

if __name__ == "__main__":
    # Called from shell: python krona_chart.py asv_table.tsv taxonomy.tsv output.html
    if len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        # Called from snakemake script directive
        main(snakemake.input["asv_table"], snakemake.input["taxonomy"], snakemake.output["krona"])
