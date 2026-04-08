#!/usr/bin/env Rscript
# =============================================================================
# Format BIG40 FA/MD/ISOVF/OD GWAS summary stats for Multivariable MR
# =============================================================================
# Run AFTER download_mvmr_sumstats.sh
# Produces exposure-formatted files for use with mv_multiple()
# Requires: data.table, TwoSampleMR
# =============================================================================

library(data.table)
library(TwoSampleMR)

# --- Configuration ---
SUMSTATS_DIR <- "mvmr_gwas_sumstats"
OUTPUT_DIR   <- "mvmr_mr_ready"
SAMPLE_SIZE  <- 33224

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# --- Mapping: metric, method, tract, IDP number ---
# Generated from BIG40 IDP table structure
tract_info <- read.csv("mvmr_additional_metrics_mapping.csv", stringsAsFactors = FALSE)

message("=== Formatting ", nrow(tract_info), " GWAS files for MVMR ===\n")

for (i in 1:nrow(tract_info)) {
  metric  <- tract_info$metric[i]
  method  <- tract_info$method[i]
  tract   <- tract_info$tract[i]
  idp_num <- tract_info$big40_idp_number[i]

  # Find the downloaded file
  fname <- paste0(metric, "_", method, "_", tract, "_IDP", idp_num, ".txt.gz")
  fpath <- file.path(SUMSTATS_DIR, fname)

  if (!file.exists(fpath)) {
    message("  WARNING: File not found: ", fpath)
    next
  }

  message("Processing: ", metric, " ", tract, " (", method, ", IDP ", idp_num, ")")

  # Read data — BIG40 stats33k format: chr rsid pos a1 a2 beta se pval(-log10)
  dat <- fread(fpath, header = TRUE)
  colnames(dat) <- c("chr", "SNP", "pos", "other_allele", "effect_allele",
                      "beta", "se", "neglog10p")

  # Convert -log10(p) to p-value
  dat[, pval := 10^(-neglog10p)]

  # Add metadata
  dat[, samplesize := SAMPLE_SIZE]

  phenotype_name <- paste0(metric, "_", tract)
  dat[, Phenotype := phenotype_name]

  # Basic QC
  dat <- dat[se < 1 & se > 0]

  # Convert to plain data.frame (TwoSampleMR requirement)
  dat <- as.data.frame(dat)

  # Format for TwoSampleMR as exposure data
  exposure_dat <- format_data(
    dat,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    samplesize_col = "samplesize",
    phenotype_col = "Phenotype",
    chr_col = "chr",
    pos_col = "pos"
  )

  # Save full formatted data
  outfile <- file.path(OUTPUT_DIR,
                       paste0(metric, "_", method, "_", tract, "_mr_ready.tsv.gz"))
  fwrite(exposure_dat, outfile, sep = "\t")
  message("  Saved: ", outfile, " (", nrow(exposure_dat), " variants)")

  rm(dat, exposure_dat)
  gc()
}

message("\n=== DONE ===")
message("MR-ready files saved to: ", OUTPUT_DIR)
message("\nThese files are used alongside the ICVF files in icvf_mr_ready/")
message("for multivariable MR with mv_multiple() in TwoSampleMR")
