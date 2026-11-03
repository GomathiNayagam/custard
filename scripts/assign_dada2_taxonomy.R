#!/usr/bin/env Rscript
# Two-step DADA2 taxonomy: assignTaxonomy (naive Bayes) + addSpecies (exact match)
#
# ---------------------------------------------------------------------------
# FIX 1 (critical) - wrong ASV IDs
# ---------------------------------------------------------------------------
# The previous version classified seqtab.rds. That table holds the PRE-chimera
# ASVs. The script then renumbered the rows ASV001..ASVn in seqtab column order.
# asv_table.tsv and asv_sequences.fasta hold the POST-chimera ASVs. They use
# their own numbering. The two numbering schemes agree only until the first
# chimera is removed.
# Measured effect on real data: 85.8% of 16S ASVs and 80.8% of rpoC ASVs
# received the taxonomy of a different ASV. 66.2% and 60.6% received a
# different genus. The fasta_file argument existed, but the script never used it.
# This version classifies the FASTA and keeps the FASTA ASV IDs.
#
# ---------------------------------------------------------------------------
# FIX 2 (critical) - addSpecies erased the naive-Bayes species column
# ---------------------------------------------------------------------------
# assignTaxonomy returns a "Species" column when the reference has 7 ranks.
# dada2::addSpecies appends a second column that is also named "Species".
# It then runs:
#     taxtab[!gen.match, "Species"] <- NA
# R matches a column name to the FIRST column with that name. This line
# therefore erased the naive-Bayes species wherever the exact-match genus
# disagreed.
# Measured effect on real data: 16S fell from 130 to 54 species calls. rpoC
# fell from 283 to 47. Both markers lost most of their valid calls.
# This version renames the naive-Bayes column BEFORE it calls addSpecies.
# The two species columns stay separate:
#     Species_naiveBayes - assignTaxonomy, bootstrap >= min_boot
#     Species_exact      - addSpecies, 100% match, genus-consistent
#
# Usage:
#   assign_dada2_taxonomy.R <fasta> <seqtab> <out> <taxa_db> [addspecies_db] [min_boot]

library(dada2)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 6) {
  stop("Usage: assign_dada2_taxonomy.R <fasta_file> <seqtab_file> <output_file> <taxa_db> [addspecies_db] [min_boot]")
}

fasta_file  <- args[1]
seqtab_file <- args[2]   # kept for interface compatibility; the script does not classify it
output_file <- args[3]
taxa_db     <- args[4]
addsp_db    <- if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL
min_boot    <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else 50

if (!file.exists(taxa_db))    stop(paste("Reference FASTA not found:", taxa_db))
if (!file.exists(fasta_file)) stop(paste("ASV FASTA not found:", fasta_file))

# --- Read the POST-chimera ASVs. These are the ASVs in asv_table.tsv. ---
asv     <- readDNAStringSet(fasta_file)
asv_ids <- sub(" .*", "", names(asv))      # ">ASV001 len=253" becomes "ASV001"
seqs    <- as.character(asv)
names(seqs) <- NULL
cat(sprintf("Loaded %d ASVs from %s\n", length(seqs), basename(fasta_file)))

if (anyDuplicated(asv_ids)) stop("The FASTA contains duplicate ASV IDs.")

cat(sprintf("Step 1: assignTaxonomy (%s), minBoot=%d...\n", basename(taxa_db), min_boot))
taxa <- assignTaxonomy(seqs, taxa_db, multithread = TRUE,
                       minBoot = min_boot, verbose = TRUE)

for (col in colnames(taxa)) {
  n <- sum(!is.na(taxa[, col]) & taxa[, col] != "")
  cat(sprintf("  %-10s: %d/%d (%.1f%%)\n", col, n, nrow(taxa), 100 * n / nrow(taxa)))
}

# Protect the naive-Bayes species column from addSpecies. See FIX 2 above.
if ("Species" %in% colnames(taxa)) {
  colnames(taxa)[colnames(taxa) == "Species"] <- "Species_naiveBayes"
}

if (!is.null(addsp_db) && file.exists(addsp_db)) {
  cat(sprintf("\nStep 2: addSpecies (%s)...\n", basename(addsp_db)))
  taxa <- addSpecies(taxa, addsp_db, verbose = TRUE)
  colnames(taxa)[colnames(taxa) == "Species"] <- "Species_exact"
  n_sp <- sum(!is.na(taxa[, "Species_exact"]) & taxa[, "Species_exact"] != "")
  cat(sprintf("  Species_exact assigned: %d/%d (%.1f%%)\n",
              n_sp, nrow(taxa), 100 * n_sp / nrow(taxa)))
} else {
  cat("\nNo addSpecies database provided. The script skips exact species matching.\n")
}

# assignTaxonomy keeps the sequences as rownames. The row order equals the input order.
if (nrow(taxa) != length(asv_ids)) {
  stop(sprintf("Row count mismatch: taxa=%d, ASVs=%d", nrow(taxa), length(asv_ids)))
}
if (!identical(rownames(taxa), seqs)) {
  stop("assignTaxonomy changed the row order. The script cannot map the ASV IDs safely.")
}

tax_df <- data.frame(ASV_ID = asv_ids, taxa, stringsAsFactors = FALSE,
                     check.names = FALSE)
tax_df$sequence <- seqs        # a downstream step can re-key by sequence if necessary

for (col in colnames(tax_df)) {
  if (col %in% c("ASV_ID", "sequence")) next
  tax_df[[col]][is.na(tax_df[[col]]) | tax_df[[col]] == ""] <- "Unclassified"
}

write.table(tax_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nSaved: %s  (%d ASVs. The IDs come from the FASTA.)\n",
            output_file, nrow(tax_df)))
