#!/usr/bin/env Rscript
# =============================================================================
# HLA Exclusion Sensitivity Analysis — Forward MR only
# =============================================================================
# Reruns forward MR after removing instruments in the HLA region (chr6:25-34Mb)
# Run from the same directory as run_bidirectional_mr.R
# =============================================================================

library(data.table)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(parallel)

# --- Command-line arguments ---
# Usage: Rscript script.R [--output-dir /path/to/output]
args <- commandArgs(trailingOnly = TRUE)
.output_dir_override <- NULL
.glioma_dir_override <- NULL
.ncores_override <- NULL
if (length(args) >= 2) {
  for (i in seq_len(length(args) - 1)) {
    if (args[i] == "--output-dir") .output_dir_override <- args[i + 1]
    if (args[i] == "--glioma-dir") .glioma_dir_override <- args[i + 1]
    if (args[i] == "--ncores") .ncores_override <- as.integer(args[i + 1])
  }
}


# =============================================================================
# CONFIGURATION (same as main script)
# =============================================================================

ICVF_DIR   <- "icvf_mr_ready"
GLIOMA_DIR <- if (!is.null(.glioma_dir_override)) .glioma_dir_override else "../20260326-GWAS_summary_stats/20260330a-results"
GLIOMA_SUBTYPES <- c("all_glioma", "idhwt", "idhmt", "idhmt_intact", "idhmt_codel")
GLIOMA_FILES <- c(
  "all_glioma/final/all_glioma_meta_summary_stats.tsv.gz",
  "idhwt/final/IDHwt_meta_summary_stats.tsv.gz",
  "idhmt/final/IDHmut_meta_summary_stats.tsv.gz",
  "idhmt_intact/final/IDHmut_1p19q_intact_meta_summary_stats.tsv.gz",
  "idhmt_codel/final/IDHmut_1p19q_codel_meta_summary_stats.tsv.gz"
)

LD_REF_PATH <- "/francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/EUR"
PLINK_BIN <- tryCatch(genetics.binaRies::get_plink_binary(),
  error = function(e) { p <- Sys.which("plink"); if (nchar(p) > 0) p else Sys.which("plink1.9") })

OUTPUT_DIR <- file.path(if (!is.null(.output_dir_override)) .output_dir_override else "mr_results", "hla_excluded")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

MIN_INSTRUMENTS <- 3
CLUMP_R2  <- 0.001
CLUMP_KB  <- 10000
# Memory-aware core limit: each mclapply fork copies the parent process memory
.mem_gb <- tryCatch({ mi <- readLines("/proc/meminfo", n = 1); as.numeric(sub(".*:\\s+(\\d+).*", "\\1", mi)) / 1024 / 1024 }, error = function(e) NA)
.mem_cores <- if (!is.na(.mem_gb)) max(4, floor(.mem_gb / 15)) else 16
N_CORES <- if (!is.null(.ncores_override)) .ncores_override else min(detectCores() - 1, .mem_cores, 60)

# HLA region definition (GRCh37)
HLA_CHR   <- 6
HLA_START <- 25000000
HLA_END   <- 34000000

# =============================================================================
# Local clump function (same as main script)
# =============================================================================

local_clump <- function(dat, snp_col = "SNP", pval_col = "pval.exposure") {
  bim <- fread(paste0(LD_REF_PATH, ".bim"), header = FALSE,
               col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  bim$chrpos <- paste0(bim$chr, ":", bim$pos)
  dat_df <- as.data.frame(dat)
  original_snps <- dat_df[[snp_col]]
  is_rsid <- grepl("^rs", original_snps[1])
  if (!is_rsid) {
    chrpos <- sub("^([^:]+:[^:]+).*", "\\1", original_snps)
    lookup <- setNames(bim$rsid, bim$chrpos)
    mapped_rsids <- lookup[chrpos]
    has_rsid <- !is.na(mapped_rsids)
    if (sum(has_rsid) < MIN_INSTRUMENTS) return(dat_df[FALSE, ])
    dat_df <- dat_df[has_rsid, ]
    dat_df$original_snp <- dat_df[[snp_col]]
    dat_df[[snp_col]] <- mapped_rsids[has_rsid]
  }
  tmpfile <- tempfile(fileext = ".tsv")
  fwrite(data.frame(SNP = dat_df[[snp_col]], P = dat_df[[pval_col]]), tmpfile, sep = "\t")
  tmpout <- tempfile()
  cmd <- sprintf("%s --bfile %s --clump %s --clump-p1 1 --clump-r2 %f --clump-kb %d --out %s --silent",
                 PLINK_BIN, LD_REF_PATH, tmpfile, CLUMP_R2, CLUMP_KB, tmpout)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  clumped_file <- paste0(tmpout, ".clumped")
  if (!file.exists(clumped_file)) { unlink(c(tmpfile, paste0(tmpout, "*"))); return(dat_df[FALSE, ]) }
  clumped <- tryCatch(fread(clumped_file), error = function(e) NULL)
  unlink(c(tmpfile, paste0(tmpout, "*")))
  if (is.null(clumped) || nrow(clumped) == 0) return(dat_df[FALSE, ])
  result <- dat_df[dat_df[[snp_col]] %in% clumped$SNP, ]
  if (!is_rsid && "original_snp" %in% names(result)) {
    result[[snp_col]] <- result$original_snp; result$original_snp <- NULL
  }
  return(result)
}

# =============================================================================
# HLA filter function
# =============================================================================

exclude_hla <- function(dat) {
  # Works with TwoSampleMR formatted data (chr.exposure / pos.exposure columns)
  chr_col <- intersect(c("chr.exposure", "chr.outcome", "chr"), names(dat))
  pos_col <- intersect(c("pos.exposure", "pos.outcome", "pos"), names(dat))

  if (length(chr_col) == 0 || length(pos_col) == 0) {
    # Try to parse from SNP column (chr:pos format)
    parts <- strsplit(dat$SNP, ":")
    dat$.chr <- as.integer(sapply(parts, `[`, 1))
    dat$.pos <- as.integer(sapply(parts, `[`, 2))
    chr_col <- ".chr"; pos_col <- ".pos"
  } else {
    chr_col <- chr_col[1]; pos_col <- pos_col[1]
  }

  in_hla <- !is.na(dat[[chr_col]]) & !is.na(dat[[pos_col]]) &
            dat[[chr_col]] == HLA_CHR &
            dat[[pos_col]] >= HLA_START &
            dat[[pos_col]] <= HLA_END

  n_removed <- sum(in_hla)
  if (n_removed > 0) message("    Removed ", n_removed, " HLA instruments (chr6:25-34Mb)")

  result <- dat[!in_hla, ]
  result$.chr <- NULL; result$.pos <- NULL
  return(result)
}

# =============================================================================
# Run forward MR with HLA exclusion
# =============================================================================

message("=== HLA Exclusion Sensitivity: Forward MR ===\n")

# Load glioma outcomes
glioma_outcomes <- list()
for (i in seq_along(GLIOMA_SUBTYPES)) {
  fpath <- file.path(GLIOMA_DIR, GLIOMA_FILES[i])
  if (!file.exists(fpath)) next
  dat <- fread(fpath, header = TRUE); dat <- as.data.frame(dat)
  dat$N <- dat$N_CASES + dat$N_CONTROLS
  dat$SNP_chrpos <- paste0(as.integer(dat$CHR), ":", dat$BP)
  out <- format_data(dat, type = "outcome",
    snp_col = "SNP_chrpos", beta_col = "BETA", se_col = "SE", pval_col = "P",
    effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "A1_FREQ",
    samplesize_col = "N", chr_col = "CHR", pos_col = "BP")
  out$outcome <- GLIOMA_SUBTYPES[i]
  glioma_outcomes[[GLIOMA_SUBTYPES[i]]] <- out
}

# Load, clump, and filter ICVF instruments
instrument_files <- list.files(ICVF_DIR, pattern = "_instruments_", full.names = TRUE)
clumped_exposures <- list()

for (inst_file in instrument_files) {
  pgs_id <- sub("_instruments.*", "", basename(inst_file))
  exposure_dat <- fread(inst_file, header = TRUE); exposure_dat <- as.data.frame(exposure_dat)
  if (nrow(exposure_dat) < MIN_INSTRUMENTS) next

  # Clump first
  clumped <- local_clump(exposure_dat, snp_col = "SNP", pval_col = "pval.exposure")
  if (nrow(clumped) < MIN_INSTRUMENTS) next

  # Convert to chr:pos
  if ("chr.exposure" %in% names(clumped) & "pos.exposure" %in% names(clumped)) {
    clumped$SNP <- paste0(as.integer(clumped$chr.exposure), ":", clumped$pos.exposure)
  }

  # EXCLUDE HLA
  message("  ", pgs_id, ": ", nrow(clumped), " instruments before HLA filter")
  clumped <- exclude_hla(clumped)
  message("    ", nrow(clumped), " instruments after HLA filter")

  if (nrow(clumped) >= MIN_INSTRUMENTS) {
    clumped_exposures[[inst_file]] <- clumped
  } else {
    message("    Skipping: too few instruments after HLA exclusion")
  }
}

# Build tasks and run in parallel
forward_tasks <- expand.grid(
  inst_file = names(clumped_exposures),
  subtype = names(glioma_outcomes),
  stringsAsFactors = FALSE
)
message("\n  Running ", nrow(forward_tasks), " forward MR tests (HLA excluded)...\n")

run_single_mr <- function(exposure_dat, outcome_dat, pgs_id, subtype) {
  result <- list(mr = NULL, het = NULL, pleio = NULL, steiger = NULL,
                 pgs_id = pgs_id, subtype = subtype, error = NULL)
  harmonised <- tryCatch(harmonise_data(exposure_dat, outcome_dat, action = 2), error = function(e) NULL)
  if (is.null(harmonised) || sum(harmonised$mr_keep) < MIN_INSTRUMENTS) {
    result$error <- "Too few SNPs"; return(result)
  }
  result$mr <- tryCatch(mr(harmonised, method_list = c("mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_weighted_mode")), error = function(e) NULL)
  result$het <- tryCatch(mr_heterogeneity(harmonised), error = function(e) NULL)
  result$pleio <- tryCatch(mr_pleiotropy_test(harmonised), error = function(e) NULL)
  result$steiger <- tryCatch(directionality_test(harmonised), error = function(e) NULL)
  return(result)
}

forward_raw <- if (nrow(forward_tasks) > 0) mclapply(seq_len(nrow(forward_tasks)), function(idx) {
  inst_file <- forward_tasks$inst_file[idx]
  subtype   <- forward_tasks$subtype[idx]
  pgs_id    <- sub("_instruments.*", "", basename(inst_file))
  run_single_mr(clumped_exposures[[inst_file]], glioma_outcomes[[subtype]], pgs_id, subtype)
}, mc.cores = N_CORES) else list()

# =============================================================================
# Collect and save
# =============================================================================

message("\n=== Saving HLA-excluded results ===\n")

mr_list <- list(); het_list <- list(); pleio_list <- list(); steiger_list <- list()
for (res in forward_raw) {
  if (is.null(res) || !is.list(res) || !is.null(res$error)) next
  if (!is.null(res$mr))      { tmp <- res$mr;      tmp$pgs_id <- res$pgs_id; tmp$subtype <- res$subtype; mr_list <- c(mr_list, list(tmp)) }
  if (!is.null(res$het))     { tmp <- res$het;     tmp$pgs_id <- res$pgs_id; tmp$subtype <- res$subtype; het_list <- c(het_list, list(tmp)) }
  if (!is.null(res$pleio))   { tmp <- res$pleio;   tmp$pgs_id <- res$pgs_id; tmp$subtype <- res$subtype; pleio_list <- c(pleio_list, list(tmp)) }
  if (!is.null(res$steiger)) { tmp <- res$steiger; tmp$pgs_id <- res$pgs_id; tmp$subtype <- res$subtype; steiger_list <- c(steiger_list, list(tmp)) }
}

if (length(mr_list) > 0) {
  fwd_df <- bind_rows(mr_list)
  fwrite(fwd_df, file.path(OUTPUT_DIR, "forward_mr_results_noHLA.tsv"), sep = "\t")
  if (length(het_list) > 0)     fwrite(bind_rows(het_list),     file.path(OUTPUT_DIR, "forward_heterogeneity_noHLA.tsv"), sep = "\t")
  if (length(pleio_list) > 0)   fwrite(bind_rows(pleio_list),   file.path(OUTPUT_DIR, "forward_pleiotropy_noHLA.tsv"), sep = "\t")
  if (length(steiger_list) > 0) fwrite(bind_rows(steiger_list), file.path(OUTPUT_DIR, "forward_steiger_noHLA.tsv"), sep = "\t")

  # Print comparison
  ivw <- fwd_df %>% filter(method == "Inverse variance weighted")
  sig <- ivw %>% filter(pval < 0.05)
  message("  Results: ", nrow(ivw), " tests, ", nrow(sig), " nominally significant")
  message("\n--- Significant IVW results (HLA excluded) ---")
  for (r in seq_len(nrow(sig))) {
    message(sprintf("  %s -> %s: b=%.4f, p=%.2e, OR=%.3f, nSNP=%d",
                    sig$pgs_id[r], sig$subtype[r], sig$b[r], sig$pval[r], exp(sig$b[r]), sig$nsnp[r]))
  }
} else {
  message("  No results produced")
}

message("\n--- Files created ---")
for (f in list.files(OUTPUT_DIR, full.names = TRUE)) message("  ", f)

message("\n--- Compare with main results ---")
message("  Main results: mr_results/forward_mr_results.tsv")
message("  HLA excluded: mr_results/hla_excluded/forward_mr_results_noHLA.tsv")
message("  If results are consistent, HLA pleiotropy is not driving the findings")
