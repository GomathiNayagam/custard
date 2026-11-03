#!/usr/bin/env Rscript
# =============================================================================
# DADA2 rpoC Metabarcoding Pipeline
# Target gene  : rpoC (RNA polymerase beta subunit — bacterial phylogenetic marker)
# Primers      : F- MAYGARAARMGNATGYTNCARGA (23 bp, degenerate)
#                R- GMCATYTGRTCNCCRTCRAA   (20 bp, degenerate)
# Amplicon     : 521 bp with primers → 478 bp insert (post-primer)
# Sequencing   : 2×300 Illumina paired-end
# Overlap      : truncF(280) + truncR(260) - 478 = 62 bp
# NOTE         : Primers MUST be removed upstream by cutadapt before this script
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(ggplot2)
})

# =============================================================================
# 0. SNAKEMAKE I/O AND PARAMETERS
# =============================================================================

# --- Inputs / Outputs ---
# fastq_dir comes from params (not input) because Snakemake passes the
# trimmed directory path separately — input holds the per-sample file list
# used only for DAG dependency resolution.
input_dir     <- snakemake@params[["fastq_dir"]]
asv_table_out <- snakemake@output[["asv_table"]]
seqtab_rds    <- snakemake@output[["seqtab_rds"]]
asv_fasta_out <- snakemake@output[["asv_fasta"]]
track_out     <- snakemake@output[["read_tracking"]]
qc_plot_out   <- snakemake@output[["qc_plots"]]

# --- From config.yaml ---
TAXA_DB      <- snakemake@config[["taxa_db"]]
THREADS      <- snakemake@config[["threads"]]
OUTPUT_DIR   <- snakemake@config[["output_dir"]]
READ1_SUFFIX <- snakemake@config[["read1_suffix"]]
READ2_SUFFIX <- snakemake@config[["read2_suffix"]]

# --- From rule params (defined in Snakemake rule, sourced from config) ---
AMPLICON_WITH_PRIMERS <- snakemake@params[["amplicon_len"]]   # 521
F_LEN                 <- snakemake@params[["primer_len_f"]]   # 23
R_LEN                 <- snakemake@params[["primer_len_r"]]   # 20
MAX_EE_F              <- snakemake@params[["maxEE_F"]]        # 2
MAX_EE_R              <- snakemake@params[["maxEE_R"]]        # 2
# Pool argument must be R logical TRUE/FALSE or string "pseudo".
# YAML booleans (true/false) arrive as R logical; quoted strings ("false") arrive as
# character and are rejected by dada(). This block normalises all three input forms.
POOL_METHOD <- local({
  p <- snakemake@params[["pool"]]
  if (is.logical(p))                          return(p)           # TRUE / FALSE
  if (tolower(as.character(p)) == "pseudo")   return("pseudo")    # "pseudo"
  if (tolower(as.character(p)) == "true")     return(TRUE)        # "true" / "True"
  if (tolower(as.character(p)) == "false")    return(FALSE)       # "false" / "False"
  stop(sprintf("Invalid pool value '%s'. Use: true, false, or pseudo", p))
})
MIN_OVERLAP           <- snakemake@params[["min_overlap"]]    # 20 (or 30 for stricter merging)
READ_LENGTH           <- snakemake@params[["read_length"]]    # 300

# --- Optional manual truncLen overrides ---
# Set truncLen_F and truncLen_R in config/params to skip auto-calculation.
# Leave as NULL (omit from config) to use auto-calculated values.
# Useful when: auto-calc gives insufficient overlap, or after reviewing qc_plots.pdf.
#
# Robust null detection — Snakemake can pass missing values in several ways:
#   Python None  → R NULL     (is.null catches this)
#   Python None  → R NA       (is.na catches this)
#   Python None  → string "NULL" or "None" (string check catches this)
#   Python None  → integer NA_integer_     (is.na catches this)
# We also suppress the warning from as.integer(NA) since we check for NA explicitly.
parse_trunc_param <- function(key) {
  tryCatch({
    v <- snakemake@params[[key]]
    # NULL — param not present or Python None passed as R NULL
    if (is.null(v)) return(NULL)
    # NA — Python None passed as R NA
    if (length(v) == 0) return(NULL)
    if (all(is.na(v))) return(NULL)
    # String "NULL", "None", "NA", "" — Python None serialised as string
    if (is.character(v) && trimws(v) %in% c("NULL", "None", "NA", "")) return(NULL)
    # Valid integer — convert and check
    int_val <- suppressWarnings(as.integer(v))
    if (is.na(int_val)) return(NULL)
    if (int_val <= 0) {
      stop(sprintf("config dada2.%s = %d is not a positive integer.", key, int_val))
    }
    int_val
  }, error = function(e) {
    message(sprintf("  Warning: could not parse dada2.%s (%s) — using auto-calculation.",
                    key, conditionMessage(e)))
    NULL
  })
}

MANUAL_TRUNC_F <- parse_trunc_param("truncLen_F")
MANUAL_TRUNC_R <- parse_trunc_param("truncLen_R")

# Diagnostic — always printed so you can confirm what the R script actually received
message(sprintf(
  "truncLen params received: F=%s  R=%s  (NULL = auto-calculate)",
  if (is.null(MANUAL_TRUNC_F)) "NULL" else as.character(MANUAL_TRUNC_F),
  if (is.null(MANUAL_TRUNC_R)) "NULL" else as.character(MANUAL_TRUNC_R)
))

# --- Derived values ---
INSERT_LEN <- AMPLICON_WITH_PRIMERS - F_LEN - R_LEN   # 521 - 23 - 20 = 478
AMP_MIN    <- INSERT_LEN - 5                           # 473
AMP_MAX    <- INSERT_LEN + 5                           # 483

message("=== DADA2 rpoC Pipeline ===")
message(sprintf("Amplicon (with primers)  : %d bp", AMPLICON_WITH_PRIMERS))
message(sprintf("Insert (post-primer)     : %d bp", INSERT_LEN))
message(sprintf("ASV length filter range  : %d-%d bp", AMP_MIN, AMP_MAX))
message(sprintf("Max raw overlap (2x%d)  : %d bp", READ_LENGTH, READ_LENGTH * 2 - INSERT_LEN))

# =============================================================================
# 1. FILE HANDLING — detect trimmed FASTQs first (needed for read length measurement)
# =============================================================================
pattern_R1_pre <- gsub("\\.", "\\\\.", READ1_SUFFIX)
pattern_R2_pre <- gsub("\\.", "\\\\.", READ2_SUFFIX)
all_files_pre  <- sort(list.files(input_dir, full.names = TRUE))
fnFs_pre       <- all_files_pre[grepl(pattern_R1_pre, all_files_pre)]
fnRs_pre       <- all_files_pre[grepl(pattern_R2_pre, all_files_pre)]

if (length(fnFs_pre) == 0) {
  stop(paste0("No trimmed FASTQ files found in: ", input_dir,
              "\n  Check cutadapt_trim completed and read1_suffix matches filenames."))
}

# Check each file has reads — plotQualityProfile crashes with a cryptic BiocParallel
# error ("arguments imply differing number of rows: 0, 1") on any empty FASTQ.
# This happens when cutadapt discards all reads from a sample (no primer detected,
# or all reads below minimum_length). Identify and remove empty files early.
check_fastq_nonempty <- function(files) {
  # sapply instead of vapply — vapply requires FUN to return exactly length 1,
  # but readLines on a 0-byte file returns character(0) (length 0), not "".
  # We coerce the result explicitly: length-0 → FALSE, anything else → nchar > 0.
  sapply(files, function(f) {
    result <- tryCatch({
      con   <- gzfile(f, "r")
      first <- readLines(con, n = 1, warn = FALSE)
      close(con)
      # length(first)==0 means 0-byte or completely empty gz file
      if (length(first) == 0) return(FALSE)
      nchar(first) > 0
    }, error = function(e) {
      message(sprintf("  Could not open %s: %s", basename(f), conditionMessage(e)))
      FALSE
    })
    result
  }, USE.NAMES = FALSE)
}

message("Checking for empty FASTQ files...")
nonempty_f <- check_fastq_nonempty(fnFs_pre)
nonempty_r <- check_fastq_nonempty(fnRs_pre)
nonempty   <- nonempty_f & nonempty_r

if (any(!nonempty)) {
  empty_samples <- gsub(READ1_SUFFIX, "", basename(fnFs_pre[!nonempty]))
  message(sprintf(
    paste0("WARNING: %d sample(s) have empty trimmed FASTQs and will be excluded: %s\n",
           "  Likely cause: cutadapt found no primer matches (check %s/*_cutadapt.log)."),
    sum(!nonempty),
    paste(empty_samples, collapse = ", "),
    file.path(dirname(fnFs_pre[1]))
  ))
}

fnFs_pre <- fnFs_pre[nonempty]
fnRs_pre <- fnRs_pre[nonempty]

if (length(fnFs_pre) == 0) {
  stop("All samples are empty after cutadapt trimming. Check cutadapt logs.")
}
if (length(fnFs_pre) < 3) {
  warning(sprintf(
    "%d sample(s) remain after removing empty files. Error learning may be unreliable.",
    length(fnFs_pre)
  ))
}

# =============================================================================
# 2. MEASURE ACTUAL READ LENGTHS FROM TRIMMED FILES
# =============================================================================
# CRITICAL: Never assume reads are READ_LENGTH bp after cutadapt.
# Degenerate primer matching, --discard-untrimmed, and sequencer-side trimming
# all produce variable-length output. We sample up to 3 files to measure.

measure_read_lengths <- function(fastq_files, n_sample = 3, reads_per_file = 10000) {
  files       <- head(fastq_files, n_sample)
  all_lengths <- integer(0)
  for (f in files) {
    con  <- gzfile(f, "r")
    i    <- 0L
    lens <- integer(reads_per_file)
    while (i < reads_per_file) {
      header <- readLines(con, n = 1, warn = FALSE)
      if (length(header) == 0) break
      seq_line <- readLines(con, n = 1, warn = FALSE)
      readLines(con, n = 1, warn = FALSE)   # +
      readLines(con, n = 1, warn = FALSE)   # qual
      i <- i + 1L
      lens[i] <- nchar(seq_line)
    }
    close(con)
    all_lengths <- c(all_lengths, lens[seq_len(i)])
  }
  return(all_lengths)
}

message("Measuring actual read lengths from trimmed FASTQs...")
len_f <- measure_read_lengths(fnFs_pre)
len_r <- measure_read_lengths(fnRs_pre)

med_f  <- as.integer(median(len_f))
med_r  <- as.integer(median(len_r))
# Use 5th percentile as the floor, NOT the absolute minimum.
# Absolute min is not robust — a single adapter dimer or PhiX read
# (e.g., 105 bp) will collapse truncLen and destroy merging.
# The 5th percentile retains 95%% of reads while ignoring extreme outliers.
p05_f  <- as.integer(quantile(len_f, 0.05))
p05_r  <- as.integer(quantile(len_r, 0.05))

message(sprintf("F reads -- median: %d bp, 5th pct: %d bp, min: %d bp, max: %d bp",
                med_f, p05_f, min(len_f), max(len_f)))
message(sprintf("R reads -- median: %d bp, 5th pct: %d bp, min: %d bp, max: %d bp",
                med_r, p05_r, min(len_r), max(len_r)))

# Print full length distribution for diagnostics
message("F read length distribution (top 10 lengths):")
len_f_tab <- sort(table(len_f), decreasing = TRUE)
print(head(len_f_tab, 10))
message("R read length distribution (top 10 lengths):")
len_r_tab <- sort(table(len_r), decreasing = TRUE)
print(head(len_r_tab, 10))

if (med_f < READ_LENGTH * 0.8) {
  message(sprintf(paste0(
    "NOTE: Median F read length (%d bp) is shorter than READ_LENGTH=%d bp.\n",
    "  Reads appear to be already primer-trimmed upstream (expected).\n",
    "  truncLen will be based on OBSERVED lengths."),
    med_f, READ_LENGTH))
}

# =============================================================================
# 3. AUTO-CALCULATE TRUNCATION LENGTHS
# =============================================================================
# Uses 5th-percentile read length as the hard floor for truncLen.
# This is robust to outlier short reads (adapter dimers, PhiX fragments)
# which would otherwise collapse truncLen to an unusable value.
#
# Logic:
#   trunc = median - quality_buffer    (trim degraded tail)
#   trunc = min(trunc, p05)            (ensure 95%% of reads pass length check)
#   enforce truncF + truncR >= INSERT_LEN + MIN_OVERLAP
#
# Your data: F median=278, R median=281, p05_r will exclude the 105 bp outliers
# Expected result: truncF~258, truncR~241, overlap~21 bp (tight but valid)

auto_truncLen <- function(med_len_f, med_len_r, p05_len_f, p05_len_r,
                          insert_len, min_overlap,
                          q_drop_buffer_f = 20, q_drop_buffer_r = 40) {

  # Step 1: apply quality buffer to median
  trunc_f <- med_len_f - q_drop_buffer_f
  trunc_r <- med_len_r - q_drop_buffer_r

  # Step 2: cap at 5th percentile — ensures 95%% of reads pass filterAndTrim
  # (reads shorter than truncLen are silently discarded)
  trunc_f <- min(trunc_f, p05_len_f)
  trunc_r <- min(trunc_r, p05_len_r)

  # Step 3: enforce minimum overlap
  needed <- insert_len + min_overlap   # e.g. 478 + 20 = 498
  total  <- trunc_f + trunc_r

  if (total < needed) {
    deficit <- needed - total
    # Boost F first (higher quality end), then R
    boost_f <- min(deficit, p05_len_f - trunc_f)
    trunc_f <- trunc_f + boost_f
    remaining <- needed - (trunc_f + trunc_r)
    if (remaining > 0) {
      trunc_r <- trunc_r + min(remaining, p05_len_r - trunc_r)
    }
  }

  trunc_f <- max(trunc_f, 1L)
  trunc_r <- max(trunc_r, 1L)
  overlap <- trunc_f + trunc_r - insert_len

  message(sprintf("Auto truncLen -> F: %d, R: %d | Expected overlap: %d bp",
                  trunc_f, trunc_r, overlap))
  message(sprintf("  (based on: F median=%d p05=%d | R median=%d p05=%d)",
                  med_len_f, p05_len_f, med_len_r, p05_len_r))

  if (overlap < 0) {
    stop(sprintf(
      paste0("FATAL: Negative overlap (%d bp) — reads are too short to span the insert.\n",
             "  truncF=%d + truncR=%d < insert_len=%d\n",
             "  Your 2x300 reads (actual ~%d/%d bp after trimming) cannot cover a %d bp insert.\n",
             "  Check: (1) amplicon_len in config is correct, ",
             "(2) cutadapt log for over-trimming."),
      overlap, trunc_f, trunc_r, insert_len, med_len_f, med_len_r, insert_len))
  } else if (overlap < min_overlap) {
    warning(sprintf(
      paste0("Overlap (%d bp) is below MIN_OVERLAP=%d.\n",
             "  Merging may fail for some read pairs.\n",
             "  Consider reducing min_overlap to %d in config if most reads merge OK."),
      overlap, min_overlap, max(overlap - 2L, 10L)))
  }

  return(c(trunc_f, trunc_r))
}

# =============================================================================
# 3b. APPLY MANUAL OVERRIDE OR AUTO-CALCULATE
# =============================================================================
# Manual override completely bypasses auto_truncLen — auto is never called
# when both truncLen_F and truncLen_R are set in config. This means a low
# p05 read length (e.g. from a short-fragment population in R2) cannot cause
# a fatal error when you have already decided on truncLen values manually.
#
# Both must be set together — setting only one raises an error immediately.

if (!is.null(MANUAL_TRUNC_F) && !is.null(MANUAL_TRUNC_R)) {

  # --- Manual override path — auto_truncLen is NOT called ---
  message(sprintf("Manual truncLen override: F=%d, R=%d — skipping auto-calculation.",
                  MANUAL_TRUNC_F, MANUAL_TRUNC_R))

  # Warn if manual value exceeds p05 (will discard >5% of reads at filterAndTrim)
  if (MANUAL_TRUNC_F > p05_f) {
    warning(sprintf(
      paste0("truncLen_F=%d exceeds the 5th-percentile F read length (%d bp).\n",
             "  More than 5%% of F reads will be discarded by filterAndTrim.\n",
             "  Consider reducing truncLen_F to <= %d."),
      MANUAL_TRUNC_F, p05_f, p05_f))
  }
  if (MANUAL_TRUNC_R > p05_r) {
    warning(sprintf(
      paste0("truncLen_R=%d exceeds the 5th-percentile R read length (%d bp).\n",
             "  More than 5%% of R reads will be discarded by filterAndTrim.\n",
             "  Consider reducing truncLen_R to <= %d."),
      MANUAL_TRUNC_R, p05_r, p05_r))
  }

  # Warn if overlap is tight
  manual_overlap <- MANUAL_TRUNC_F + MANUAL_TRUNC_R - INSERT_LEN
  if (manual_overlap < MIN_OVERLAP) {
    warning(sprintf(
      paste0("truncLen_F=%d + truncLen_R=%d gives only %d bp overlap ",
             "(INSERT_LEN=%d, MIN_OVERLAP=%d).\n",
             "  Merging may fail. Increase truncLen values or reduce min_overlap."),
      MANUAL_TRUNC_F, MANUAL_TRUNC_R, manual_overlap, INSERT_LEN, MIN_OVERLAP))
  }

  message(sprintf("  Manual truncLen: F=%d, R=%d | Overlap=%d bp",
                  MANUAL_TRUNC_F, MANUAL_TRUNC_R, manual_overlap))

  TRUNC_F <- MANUAL_TRUNC_F
  TRUNC_R <- MANUAL_TRUNC_R

} else if (!is.null(MANUAL_TRUNC_F) && is.null(MANUAL_TRUNC_R)) {
  stop("truncLen_F is set but truncLen_R is missing. Set both or neither in config.")

} else if (is.null(MANUAL_TRUNC_F) && !is.null(MANUAL_TRUNC_R)) {
  stop("truncLen_R is set but truncLen_F is missing. Set both or neither in config.")

} else {
  # --- Auto-calculation path — only reached when both are NULL ---
  message("No manual truncLen — running auto-calculation...")
  trunc_lens <- auto_truncLen(
    med_len_f   = med_f,    med_len_r   = med_r,
    p05_len_f   = p05_f,    p05_len_r   = p05_r,
    insert_len  = INSERT_LEN,
    min_overlap = MIN_OVERLAP
  )
  TRUNC_F <- trunc_lens[1]
  TRUNC_R <- trunc_lens[2]
}

message(sprintf("Final truncLen: F=%d, R=%d | Overlap=%d bp | Mode=%s",
                TRUNC_F, TRUNC_R,
                TRUNC_F + TRUNC_R - INSERT_LEN,
                ifelse(!is.null(MANUAL_TRUNC_F), "MANUAL", "AUTO")))

# =============================================================================
# 4. FILE HANDLING — reuse pre-measured file lists, extract sample names
# =============================================================================
# fnFs_pre / fnRs_pre already validated above (length > 0, used for read measurement)
fnFs         <- fnFs_pre
fnRs         <- fnRs_pre
stopifnot("R1 and R2 file counts do not match" = length(fnFs) == length(fnRs))
sample.names <- gsub(READ1_SUFFIX, "", basename(fnFs))
message(sprintf("Found %d samples", length(sample.names)))

# =============================================================================
# 5. QUALITY PROFILE PLOTS
# Dashed lines show where truncation will be applied.
# If quality drops BEFORE the line → decrease truncLen in config.
# If quality is still good AT the line → you may increase truncLen.
# =============================================================================
message("Generating quality profile plots...")
# Use up to 6 non-empty files for QC plots
qc_idx <- seq_len(min(6, length(fnFs)))

# safe_qc_plot: wraps plotQualityProfile in tryCatch so a single
# problematic file (very few reads, truncated FASTQ) doesn't abort
# the entire pipeline — it just skips that file and logs a warning.
safe_qc_plot <- function(files, trunc_pos, color, label) {
  plots <- list()
  for (f in files) {
    p <- tryCatch(
      plotQualityProfile(f),
      error = function(e) {
        message(sprintf("  Skipping QC plot for %s: %s", basename(f), conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(p)) plots[[length(plots) + 1]] <- p
  }
  if (length(plots) == 0) {
    message(sprintf("WARNING: No valid %s quality profiles to plot.", label))
    return(invisible(NULL))
  }
  # Combine plots using patchwork if available, else plot individually
  combined <- tryCatch(
    {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(plots) +
          patchwork::plot_annotation(
            title = sprintf("%s reads — pre-filter | trunc%s = %d", label, label, trunc_pos)
          )
      } else {
        plots[[1]] +
          ggtitle(sprintf("%s reads — pre-filter | trunc%s = %d", label, label, trunc_pos))
      }
    },
    error = function(e) plots[[1]]
  )
  print(
    combined +
      geom_vline(xintercept = trunc_pos, color = color,
                 linetype = "dashed", linewidth = 0.8) +
      annotate("text", x = trunc_pos + 3, y = 10,
               label = paste0("trunc", label, " = ", trunc_pos),
               color = color, hjust = 0, size = 3.5)
  )
}

pdf(qc_plot_out, width = 11, height = 8)
  safe_qc_plot(fnFs[qc_idx], TRUNC_F, "red",  "F")
  safe_qc_plot(fnRs[qc_idx], TRUNC_R, "blue", "R")
dev.off()
message(sprintf("Quality profiles saved -> %s", qc_plot_out))

# =============================================================================
# 6. FILTER AND TRIM
# trimLeft = c(0,0): primers already removed upstream by cutadapt
# rm.phix  = TRUE : removes PhiX spike-in contamination
# =============================================================================
filt_dir <- file.path(OUTPUT_DIR, "filtered")
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))

message(sprintf("Filtering: truncLen=c(%d,%d), maxEE=c(%s,%s), threads=%d",
                TRUNC_F, TRUNC_R, MAX_EE_F, MAX_EE_R, THREADS))

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  trimLeft    = c(0, 0),
  truncLen    = c(TRUNC_F, TRUNC_R),
  maxN        = 0,
  maxEE       = c(MAX_EE_F, MAX_EE_R),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = THREADS
)

filter_summary <- data.frame(
  sample    = sample.names,
  reads_in  = out[, "reads.in"],
  reads_out = out[, "reads.out"],
  pct_pass  = round(100 * out[, "reads.out"] / pmax(out[, "reads.in"], 1), 1)
)
message("Filtering summary:")
print(filter_summary)

low_filter <- filter_summary$pct_pass < 50 & filter_summary$reads_in > 0
if (any(low_filter)) {
  warning(sprintf(
    "Low filtering yield (<50%%) in: %s\n  Consider increasing maxEE_R in config.",
    paste(sample.names[low_filter], collapse = ", ")
  ))
}

# Drop samples with no reads surviving
exists       <- file.exists(filtFs) & file.exists(filtRs)
if (any(!exists)) {
  message(sprintf("Dropping %d empty sample(s): %s",
                  sum(!exists), paste(sample.names[!exists], collapse = ", ")))
}
filtFs       <- filtFs[exists]
filtRs       <- filtRs[exists]
sample.names <- sample.names[exists]
out          <- out[exists, , drop = FALSE]

if (length(filtFs) == 0) stop("No samples remain after filtering — check truncLen and maxEE.")
if (length(filtFs) < 3)  warning("Fewer than 3 samples — error model may be unreliable.")

# =============================================================================
# 7. LEARN ERROR RATES
# =============================================================================
message("Learning error rates...")
errF <- learnErrors(filtFs, multithread = THREADS)
errR <- learnErrors(filtRs, multithread = THREADS)

pdf(sub("\\.pdf$", "_error_rates.pdf", qc_plot_out), width = 11, height = 8)
  print(plotErrors(errF, nominalQ = TRUE) + ggtitle("Error rates — Forward reads"))
  print(plotErrors(errR, nominalQ = TRUE) + ggtitle("Error rates — Reverse reads"))
dev.off()
# Good fit: black dots (observed) closely follow the fitted black line

# =============================================================================
# 8. DEREPLICATE
# =============================================================================
message("Dereplicating...")
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# =============================================================================
# 9. DENOISE
# pool = "pseudo" : good balance — catches rare taxa without full memory cost
# pool = TRUE     : use for <20 samples or very low-diversity communities
# pool = FALSE    : use for >100 samples to save memory
# =============================================================================
message(sprintf("Denoising (pool='%s')...", POOL_METHOD))
dadaFs <- dada(derepFs, err = errF, pool = POOL_METHOD, multithread = THREADS)
dadaRs <- dada(derepRs, err = errR, pool = POOL_METHOD, multithread = THREADS)

# =============================================================================
# 10. MERGE PAIRED READS
# minOverlap  = 20 : minimum reliable overlap
# maxMismatch = 0  : strict — any mismatch in overlap indicates error/chimera
# =============================================================================
message("Merging paired reads...")
mergers <- mergePairs(
  dadaFs, derepFs, dadaRs, derepRs,
  minOverlap  = MIN_OVERLAP,
  maxMismatch = 0,
  verbose     = TRUE
)

# =============================================================================
# 11. SEQUENCE TABLE + LENGTH FILTER
# =============================================================================
seqtab <- makeSequenceTable(mergers)

message(sprintf("ASV length distribution (should peak near %d bp):", INSERT_LEN))
print(table(nchar(getSequences(seqtab))))

seqtab_filt <- seqtab[, nchar(colnames(seqtab)) %in% AMP_MIN:AMP_MAX]
n_removed   <- ncol(seqtab) - ncol(seqtab_filt)
message(sprintf("Length filter: removed %d off-target ASVs (outside %d-%d bp)",
                n_removed, AMP_MIN, AMP_MAX))
seqtab <- seqtab_filt

if (ncol(seqtab) == 0) {
  stop(sprintf(
    "No ASVs remain after length filtering (%d-%d bp). Check amplicon_len in config.",
    AMP_MIN, AMP_MAX
  ))
}

saveRDS(seqtab, seqtab_rds)

# =============================================================================
# 10. CHIMERA REMOVAL
# =============================================================================
message("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = THREADS)

n_chim             <- ncol(seqtab) - ncol(seqtab.nochim)
pct_reads_retained <- round(100 * sum(seqtab.nochim) / sum(seqtab), 1)
message(sprintf("Chimera removal: %d chimeric ASVs removed, %.1f%% reads retained",
                n_chim, pct_reads_retained))

if (pct_reads_retained < 50) {
  warning("Less than 50% reads retained after chimera removal.\n",
          "  Possible causes: primers not fully removed by cutadapt, ",
          "wrong amplicon_len, or degraded samples.")
}

# =============================================================================
# 11. TAXONOMY — DADA2 assignTaxonomy with custom rpoC DB
# =============================================================================
message(sprintf("Assigning taxonomy using: %s", TAXA_DB))

taxa <- assignTaxonomy(
  seqtab.nochim,
  refFasta    = TAXA_DB,
  multithread = THREADS,
  minBoot     = 80,
  verbose     = TRUE
)

n_unclassified <- sum(is.na(taxa[, 1]))
message(sprintf("Taxonomy: %d / %d ASVs classified (%.1f%% unclassified at rank 1)",
                nrow(taxa) - n_unclassified, nrow(taxa),
                100 * n_unclassified / nrow(taxa)))

taxa_out <- file.path(OUTPUT_DIR, "taxonomy.tsv")
write.table(
  cbind(ASV = paste0("ASV", formatC(seq_len(nrow(taxa)),
                                    width = nchar(nrow(taxa)), flag = "0")),
        taxa),
  taxa_out, sep = "\t", quote = FALSE, row.names = FALSE
)
message(sprintf("Taxonomy table -> %s", taxa_out))

# =============================================================================
# 12. READ TRACKING
# =============================================================================
getN <- function(x) sum(getUniques(x))

track <- data.frame(
  sample    = sample.names,
  input     = out[, "reads.in"],
  filtered  = out[, "reads.out"],
  denoisedF = sapply(dadaFs,  getN),
  denoisedR = sapply(dadaRs,  getN),
  merged    = sapply(mergers, getN),
  nonchim   = rowSums(seqtab.nochim),
  stringsAsFactors = FALSE
)
track$pct_final <- round(100 * track$nonchim / pmax(track$input, 1), 1)

message("Read tracking:")
print(track)

low_yield <- track$pct_final < 10
if (any(low_yield)) {
  warning(sprintf("Very low yield (<10%%) in: %s",
                  paste(track$sample[low_yield], collapse = ", ")))
}

write.table(track, track_out, sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# 13. EXPORT ASV TABLE AND FASTA
# =============================================================================
asv_seqs_final <- colnames(seqtab.nochim)
n_asvs         <- length(asv_seqs_final)
asv_headers    <- paste0("ASV", formatC(seq_len(n_asvs),
                                        width = nchar(n_asvs), flag = "0"))

asv_tab           <- t(seqtab.nochim)
colnames(asv_tab) <- sample.names
rownames(asv_tab) <- asv_headers
write.table(asv_tab, asv_table_out, sep = "\t", quote = FALSE, col.names = NA)

asv_fasta        <- DNAStringSet(asv_seqs_final)
names(asv_fasta) <- paste0(asv_headers, " len=", nchar(asv_seqs_final))
writeXStringSet(asv_fasta, asv_fasta_out)

message(sprintf("ASV table -> %s  (%d ASVs x %d samples)",
                asv_table_out, n_asvs, length(sample.names)))
message(sprintf("FASTA     -> %s", asv_fasta_out))

# =============================================================================
# 14. FINAL SUMMARY
# =============================================================================
message("\n========== PIPELINE COMPLETE ==========")
message(sprintf("  Samples processed   : %d", length(sample.names)))
message(sprintf("  Final ASV count     : %d", n_asvs))
message(sprintf("  Median reads/sample : %d", as.integer(median(track$nonchim))))
message(sprintf("  Median %% retained  : %.1f%%", median(track$pct_final)))
message(sprintf("  TruncLen used       : F=%d, R=%d (overlap = %d bp)",
                TRUNC_F, TRUNC_R, TRUNC_F + TRUNC_R - INSERT_LEN))
message(sprintf("  Taxonomy DB         : %s", basename(TAXA_DB)))
message("=========================================")
message("NEXT STEPS:")
message("  1. Review QC plots — adjust truncLen params if quality drop misaligns")
message("  2. Check read tracking — investigate samples with <10% final retention")
message("  3. If method=kraken2 in config, pass asv_fasta to Kraken2 rule")
message("  4. Pass asv_fasta to KO inference rule for functional profiling")
