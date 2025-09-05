#!/usr/bin/env python3

import sys
import pandas as pd
from ete3 import NCBITaxa

# Initialize NCBI taxonomy database
ncbi = NCBITaxa()

def get_taxonomy_lineage(taxid):
    """Retrieve full taxonomic lineage for a given taxid."""
    try:
        # Get lineage as list of taxids
        lineage = ncbi.get_lineage(taxid)
        if not lineage:
            return ["Unclassified"] + ["NA"] * 6

        # Map taxids to rank names
        rank_dict = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)

        # Initialize ranks
        ranks = {
            "superkingdom": "NA",
            "phylum": "NA",
            "class": "NA",
            "order": "NA",
            "family": "NA",
            "genus": "NA",
            "species": "NA"
        }

        # Fill in ranks
        for tid in lineage:
            rank = rank_dict.get(tid, "no rank")
            if rank in ranks and tid in names:
                ranks[rank] = names[tid].replace("uncultured ", "").split(" (taxid")[0].strip()

        return [
            ranks["superkingdom"],
            ranks["phylum"],
            ranks["class"],
            ranks["order"],
            ranks["family"],
            ranks["genus"],
            ranks["species"]
        ]
    except Exception:
        return ["Unclassified"] + ["NA"] * 6

def parse_kraken_output(raw_output, output_table, kraken_db):
    """Parse Kraken2 raw output and generate taxonomy table."""
    # Read raw Kraken2 output
    data = []
    with open(raw_output, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[0] == "C":  # Classified sequences
                asv_id = fields[1]
                tax_name = fields[2]
                # Extract taxid from name like "Bacteroidetes oral taxon 274 str. F0058 (taxid 575590)"
                taxid = tax_name.split("taxid ")[-1].rstrip(")") if "taxid" in tax_name else "0"
                try:
                    taxid = int(taxid)
                except ValueError:
                    taxid = 0
                lineage = get_taxonomy_lineage(taxid)
                data.append([asv_id] + lineage)
            else:
                # Unclassified sequences
                asv_id = fields[1]
                data.append([asv_id, "Unclassified"] + ["NA"] * 6)

    # Create DataFrame
    columns = ["ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    df = pd.DataFrame(data, columns=columns)

    # Save to TSV
    df.to_csv(output_table, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: parse_kraken2_taxonomy.py <kraken_raw_output> <output_taxonomy_table> <kraken_db>")
        sys.exit(1)

    raw_output = sys.argv[1]
    output_table = sys.argv[2]
    kraken_db = sys.argv[3]  # Not used directly but kept for compatibility
    parse_kraken_output(raw_output, output_table, kraken_db)
