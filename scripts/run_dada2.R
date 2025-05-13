#!/usr/bin/env Rscript

library(dada2)
library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- snakemake@input[["fastq_dir"]]
taxa_db <- snakemake@input[["taxa_db"]]
asv_table_out <- snakemake@output[["asv_table"]]
taxa_out <- snakemake@output[["taxa"]]
seqtab_rds <- snakemake@output[["seqtab_rds"]]
taxa_rds <- snakemake@output[["taxa_rds"]]

# 1. File handling
fns <- sort(list.files(input_dir, pattern = "_R[12].fastq.gz$", full.names = TRUE))
fnFs <- fns[grepl("_R1", fns)]
fnRs <- fns[grepl("_R2", fns)]

sample.names <- gsub("_R1.*", "", basename(fnFs))

# 2. Filtering
filtFs <- file.path("results", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("results", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
dir.create("results/filtered", showWarnings = FALSE)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240,200), 
                     maxN=0, maxEE=c(2,2), truncQ=2, 
                     compress=TRUE, multithread=TRUE)

# 3. Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# 4. Dereplication
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# 5. Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# 6. Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# 7. Make sequence table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, seqtab_rds)

# 8. Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

# 9. Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, taxa_db, multithread=TRUE)
saveRDS(taxa, taxa_rds)

# 10. Export ASV table
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0("ASV", seq_along(asv_seqs))

asv_tab <- t(seqtab.nochim)
colnames(asv_tab) <- sample.names
rownames(asv_tab) <- asv_headers

write.table(asv_tab, asv_table_out, sep="\t", quote=FALSE, col.names=NA)

# 11. Export taxonomy
taxa_df <- as.data.frame(taxa)
rownames(taxa_df) <- asv_headers
write.table(taxa_df, taxa_out, sep="\t", quote=FALSE, col.names=NA)

# 12. Save FASTA of ASVs for Step 2
asv_fasta <- Biostrings::DNAStringSet(asv_seqs)
names(asv_fasta) <- asv_headers
Biostrings::writeXStringSet(asv_fasta, "results/asv_sequences.fasta")
