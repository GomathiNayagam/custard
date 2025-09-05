#!/usr/bin/env Rscript

# Load required libraries
library(dada2)
library(Biostrings)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: assign_dada2_taxonomy.R <fasta_file> <seqtab_file> <output_file> <taxa_db>")
}

fasta_file <- args[1]
seqtab_file <- args[2]
output_file <- args[3]
taxa_db <- args[4]

# Check if database path is valid
if (!file.exists(taxa_db)) {
  stop(paste("Reference FASTA file does not exist:", taxa_db))
}

# Read ASV sequences
seqs <- readDNAStringSet(fasta_file)
asv_ids <- names(seqs)

# Assign taxonomy
tax <- assignTaxonomy(seqs, taxa_db, multithread = TRUE, verbose = TRUE)

# Format taxonomy table
tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_df <- as.data.frame(matrix(NA, nrow = length(asv_ids), ncol = length(tax_levels) + 1))
colnames(tax_df) <- c("ASV_ID", tax_levels)
tax_df$ASV_ID <- asv_ids

for (i in 1:nrow(tax)) {
  for (j in seq_along(tax_levels)) {
    value <- tax[i, j]
    tax_df[i, tax_levels[j]] <- ifelse(is.na(value), "NA", value)
  }
}

# Fill unclassified kingdom entries
tax_df$Kingdom[is.na(tax_df$Kingdom) | tax_df$Kingdom == ""] <- "Unclassified"

# Write to file
write.table(tax_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

