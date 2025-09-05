#!/usr/bin/env Rscript

library(dada2)
library(tidyverse)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

input_dir <- snakemake@input[["fastq_dir"]]
asv_table_out <- snakemake@output[["asv_table"]]
seqtab_rds <- snakemake@output[["seqtab_rds"]]
asv_fasta_out <- snakemake@output[["asv_fasta"]]

# 1. File handling
fns <- sort(list.files(input_dir, pattern = "_R[12].fastq.gz$", full.names = TRUE))
fnFs <- fns[grepl("_R1", fns)]
fnRs <- fns[grepl("_R2", fns)]
sample.names <- gsub("_R1.*", "", basename(fnFs))

# 2. Filtering
filtFs <- file.path("results", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("results", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
dir.create("results/filtered", showWarnings = FALSE)

filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              truncLen=c(0, 0), maxN=0, maxEE=c(2,2), truncQ=2, 
              compress=TRUE, multithread=TRUE)

# 3–6. DADA2 pipeline
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# 7. Sequence table and chimera removal
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, seqtab_rds)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

# 8. Export ASV table
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0("ASV", seq_along(asv_seqs))

asv_tab <- t(seqtab.nochim)
colnames(asv_tab) <- sample.names
rownames(asv_tab) <- asv_headers
write.table(asv_tab, asv_table_out, sep="\t", quote=FALSE, col.names=NA)

# 9. Export ASVs as FASTA
asv_fasta <- Biostrings::DNAStringSet(asv_seqs)
names(asv_fasta) <- asv_headers
Biostrings::writeXStringSet(asv_fasta, asv_fasta_out)

