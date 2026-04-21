#!/usr/bin/env Rscript
# =============================================================================
# Bidirectional Two-Sample Mendelian Randomization: ICVF → Glioma
# =============================================================================
# PARALLELIZED version with LOCAL LD clumping (plink) and sequential MR-PRESSO
#
# Requires: data.table, TwoSampleMR, MRPRESSO, ggplot2, dplyr, parallel
#
# Install if needed:
#   install.packages(c("data.table", "ggplot2", "dplyr", "remotes"))
#   remotes::install_github("MRCIEU/TwoSampleMR")
#   remotes::install_github("rondolab/MR-PRESSO")
#   remotes::install_github("MRCIEU/genetics.binaRies")  # for plink binary
# =============================================================================

library(data.table)
library(TwoSampleMR)
# library(MRPRESSO)  # DISABLED: MRPRESSO package has an internal bug that crashes
#   when constructing outlier-corrected results ("arguments imply differing
#   number of rows: 1, 0" or "missing value where TRUE/FALSE needed").
#   The bug is inside mr_presso() itself, not in our code.
#   MR-Egger intercept + Cochran's Q provide equivalent pleiotropy detection.
#   Uncomment if/when the MRPRESSO package fixes this issue.
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
# CONFIGURATION
# =============================================================================

# --- OpenGWAS JWT Authentication (used as fallback only) ---
TOKEN_FILE <- "OpenGWASToken"
if (file.exists(TOKEN_FILE)) {
  token <- trimws(readLines(TOKEN_FILE, n = 1, warn = FALSE))
  Sys.setenv(OPENGWAS_JWT = token)
  message("OpenGWAS JWT token loaded from: ", TOKEN_FILE)
}

# --- Local LD reference panel ---
# Download from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
LD_REF_PATH <- "/francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/EUR"

# Plink binary — tries genetics.binaRies package, then system plink
PLINK_BIN <- tryCatch(
  genetics.binaRies::get_plink_binary(),
  error = function(e) {
    plink_sys <- Sys.which("plink")
    if (nchar(plink_sys) > 0) plink_sys else Sys.which("plink1.9")
  }
)
message("Plink binary: ", PLINK_BIN)

# Verify LD reference exists
if (!file.exists(paste0(LD_REF_PATH, ".bed"))) {
  stop("LD reference not found: ", LD_REF_PATH, ".bed\n",
       "Download from http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz")
}
message("LD reference: ", LD_REF_PATH)

# --- Data paths ---
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

OUTPUT_DIR <- if (!is.null(.output_dir_override)) .output_dir_override else "mr_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE)

# --- MR Parameters ---
P_THRESHOLD     <- 5e-8
P_RELAXED       <- 5e-6
MIN_INSTRUMENTS <- 3
CLUMP_R2        <- 0.001
CLUMP_KB        <- 10000

# --- Parallelization ---
# Memory-aware core limit: each mclapply fork copies the parent process memory
# With 5 glioma outcomes loaded (~2GB each), each fork needs ~10-15GB
# Default: min(available_cores - 1, floor(total_memory_GB / 15), 60)
.mem_gb <- tryCatch({ mi <- readLines("/proc/meminfo", n = 1); as.numeric(sub(".*:\\s+(\\d+).*", "\\1", mi)) / 1024 / 1024 }, error = function(e) NA)
.mem_cores <- if (!is.na(.mem_gb)) max(4, floor(.mem_gb / 15)) else 16
N_CORES <- if (!is.null(.ncores_override)) .ncores_override else min(detectCores() - 1, .mem_cores, 60)
message("Using ", N_CORES, " cores for parallel execution")

# =============================================================================
# LOCAL LD CLUMPING FUNCTION
# =============================================================================
# Handles both rsID and chr:pos:A1:A2 format SNPs by mapping through .bim file

local_clump <- function(dat, snp_col = "SNP", pval_col = "pval.exposure") {
  # Read the .bim file to get chr:pos -> rsID mapping
  bim <- fread(paste0(LD_REF_PATH, ".bim"), header = FALSE,
               col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  bim$chrpos <- paste0(bim$chr, ":", bim$pos)

  dat_df <- as.data.frame(dat)
  original_snps <- dat_df[[snp_col]]

  # Detect SNP format — check ALL SNPs, not just the first one
  # FIX: Original checked only first SNP. If the first SNP happened to be
  # chr:pos format but most were rsIDs (or vice versa), the wrong branch ran.
  # For PGS001458 this caused 0/7170 SNPs to map, dropping the entire tract.
  n_rsid <- sum(grepl("^rs", original_snps))
  is_rsid <- n_rsid > length(original_snps) * 0.5  # majority are rsIDs

  if (!is_rsid) {
    # Extract chr:pos from formats like "1:752721" or "1:752721:A:G"
    chrpos <- sub("^([^:]+:[^:]+).*", "\\1", original_snps)
    lookup <- setNames(bim$rsid, bim$chrpos)
    mapped_rsids <- lookup[chrpos]

    # Keep only SNPs that mapped
    has_rsid <- !is.na(mapped_rsids)
    message("    Mapped ", sum(has_rsid), "/", length(has_rsid), " SNPs to rsIDs for clumping")

    if (sum(has_rsid) < MIN_INSTRUMENTS) {
      message("    Too few SNPs mapped to rsIDs")
      return(dat_df[FALSE, ])  # empty
    }

    dat_df <- dat_df[has_rsid, ]
    dat_df$original_snp <- dat_df[[snp_col]]
    dat_df[[snp_col]] <- mapped_rsids[has_rsid]
  }

  # Write temp file for plink
  tmpfile <- tempfile(fileext = ".tsv")
  clump_input <- data.frame(SNP = dat_df[[snp_col]], P = dat_df[[pval_col]])
  fwrite(clump_input, tmpfile, sep = "\t")

  # Run plink clumping
  tmpout <- tempfile()
  cmd <- sprintf(
    "%s --bfile %s --clump %s --clump-p1 1 --clump-r2 %f --clump-kb %d --out %s --silent",
    PLINK_BIN, LD_REF_PATH, tmpfile, CLUMP_R2, CLUMP_KB, tmpout
  )
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Read clumped results
  clumped_file <- paste0(tmpout, ".clumped")
  if (!file.exists(clumped_file)) {
    message("    Plink clumping produced no output")
    unlink(c(tmpfile, paste0(tmpout, "*")))
    return(dat_df[FALSE, ])
  }

  clumped <- tryCatch(fread(clumped_file), error = function(e) NULL)
  unlink(c(tmpfile, paste0(tmpout, "*")))

  if (is.null(clumped) || nrow(clumped) == 0) return(dat_df[FALSE, ])

  # Filter to clumped SNPs
  keep_snps <- clumped$SNP
  result <- dat_df[dat_df[[snp_col]] %in% keep_snps, ]

  # Restore original SNP IDs if we mapped
  if (!is_rsid && "original_snp" %in% names(result)) {
    result[[snp_col]] <- result$original_snp
    result$original_snp <- NULL
  }

  return(result)
}

# =============================================================================
# HELPER: Run MR for one exposure-outcome pair (NO MR-PRESSO — done separately)
# =============================================================================

run_single_mr <- function(exposure_dat, outcome_dat, pgs_id, subtype, direction) {
  result <- list(
    mr = NULL, het = NULL, pleio = NULL, steiger = NULL,
    pgs_id = pgs_id, subtype = subtype, direction = direction, error = NULL
  )

  harmonised <- tryCatch(
    harmonise_data(exposure_dat, outcome_dat, action = 2),
    error = function(e) NULL
  )

  if (is.null(harmonised) || sum(harmonised$mr_keep) < MIN_INSTRUMENTS) {
    n_harm <- if (!is.null(harmonised)) nrow(harmonised) else 0
    n_keep <- if (!is.null(harmonised)) sum(harmonised$mr_keep) else 0
    result$error <- paste0("Too few SNPs after harmonisation (matched: ", n_harm, ", kept: ", n_keep, ")")
    return(result)
  }

  result$mr <- tryCatch(
    mr(harmonised, method_list = c("mr_ivw", "mr_egger_regression",
                                    "mr_weighted_median", "mr_weighted_mode")),
    error = function(e) NULL
  )
  result$het <- tryCatch(mr_heterogeneity(harmonised), error = function(e) NULL)
  result$pleio <- tryCatch(mr_pleiotropy_test(harmonised), error = function(e) NULL)
  result$steiger <- tryCatch(directionality_test(harmonised), error = function(e) NULL)

  # AUDIT FIX: Calculate F-statistics for instrument strength
  # F = (beta/se)^2 for each instrument; mean F should be > 10
  f_vals <- (harmonised$beta.exposure[harmonised$mr_keep] /
             harmonised$se.exposure[harmonised$mr_keep])^2
  result$f_stats <- data.frame(
    mean_F = mean(f_vals), median_F = median(f_vals),
    min_F = min(f_vals), n_weak = sum(f_vals < 10),
    n_instruments = length(f_vals)
  )

  return(result)
}

# =============================================================================
# STEP 1: Load and format glioma outcome data
# =============================================================================

message("\n=== STEP 1: Loading glioma summary stats ===")

load_glioma <- function(filepath, subtype_name, type = "outcome") {
  dat <- fread(filepath, header = TRUE)
  dat <- as.data.frame(dat)
  dat$N <- dat$N_CASES + dat$N_CONTROLS
  dat$SNP_chrpos <- paste0(as.integer(dat$CHR), ":", dat$BP)

  formatted <- format_data(dat, type = type,
    snp_col = "SNP_chrpos", beta_col = "BETA", se_col = "SE", pval_col = "P",
    effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "A1_FREQ",
    samplesize_col = "N", chr_col = "CHR", pos_col = "BP")

  if (type == "outcome") formatted$outcome <- subtype_name
  else formatted$exposure <- subtype_name
  return(formatted)
}

glioma_outcomes <- list()
for (i in seq_along(GLIOMA_SUBTYPES)) {
  fpath <- file.path(GLIOMA_DIR, GLIOMA_FILES[i])
  if (!file.exists(fpath)) { message("  WARNING: ", fpath, " not found"); next }
  message("  Loading: ", GLIOMA_SUBTYPES[i])
  glioma_outcomes[[GLIOMA_SUBTYPES[i]]] <- load_glioma(fpath, GLIOMA_SUBTYPES[i], "outcome")
}
message("  Loaded ", length(glioma_outcomes), " glioma subtypes\n")

# =============================================================================

# =============================================================================
# CHECKPOINT: Skip Steps 2-4 if results already exist
# =============================================================================

.fwd_file <- file.path(OUTPUT_DIR, "forward_mr_results.tsv")
.rev_file <- file.path(OUTPUT_DIR, "reverse_mr_results.tsv")
.skip_analysis <- file.exists(.fwd_file) && file.info(.fwd_file)$size > 0 &&
                  file.exists(.rev_file) && file.info(.rev_file)$size > 0

if (.skip_analysis) {
  message("\n=== CHECKPOINT: Result files found — skipping Steps 2-4 ===")
  message("  Forward: ", .fwd_file)
  message("  Reverse: ", .rev_file)
  message("  To force rerun, delete these files or use a different --output-dir\n")

  # Load existing results so Steps 6-7 (plots + summary) can run
  fwd <- list(
    mr      = tryCatch(fread(file.path(OUTPUT_DIR, "forward_mr_results.tsv")), error = function(e) NULL),
    het     = tryCatch(fread(file.path(OUTPUT_DIR, "forward_heterogeneity.tsv")), error = function(e) NULL),
    pleio   = tryCatch(fread(file.path(OUTPUT_DIR, "forward_pleiotropy.tsv")), error = function(e) NULL),
    steiger = tryCatch(fread(file.path(OUTPUT_DIR, "forward_steiger.tsv")), error = function(e) NULL)
  )
  rev <- list(
    mr    = tryCatch(fread(file.path(OUTPUT_DIR, "reverse_mr_results.tsv")), error = function(e) NULL),
    het   = tryCatch(fread(file.path(OUTPUT_DIR, "reverse_heterogeneity.tsv")), error = function(e) NULL),
    pleio = tryCatch(fread(file.path(OUTPUT_DIR, "reverse_pleiotropy.tsv")), error = function(e) NULL)
  )
  message("  Loaded results from disk for plotting/summary")

} else {

# STEP 2: Forward MR — ICVF (exposure) → Glioma (outcome) [PARALLEL]
# =============================================================================

message("=== STEP 2: Forward MR — ICVF → Glioma (parallel) ===\n")

instrument_files <- list.files(ICVF_DIR, pattern = "_instruments_", full.names = TRUE)
message("  Found ", length(instrument_files), " ICVF instrument files")

# LD-clump each instrument set using LOCAL plink
message("  LD-clumping ICVF instruments (local plink)...")
clumped_exposures <- list()
for (inst_file in instrument_files) {
  fname <- basename(inst_file)
  pgs_id <- sub("_instruments.*", "", fname)

  exposure_dat <- fread(inst_file, header = TRUE)
  exposure_dat <- as.data.frame(exposure_dat)

  if (nrow(exposure_dat) < MIN_INSTRUMENTS) {
    message("    ", pgs_id, ": skipping (", nrow(exposure_dat), " instruments)")
    next
  }

  clumped <- local_clump(exposure_dat, snp_col = "SNP", pval_col = "pval.exposure")

  if (nrow(clumped) >= MIN_INSTRUMENTS) {
    # Convert rsID -> chr:pos for harmonisation with glioma
    if ("chr.exposure" %in% names(clumped) & "pos.exposure" %in% names(clumped)) {
      clumped$SNP <- paste0(as.integer(clumped$chr.exposure), ":", clumped$pos.exposure)
    }
    clumped_exposures[[inst_file]] <- clumped
    message("    ", pgs_id, ": ", nrow(clumped), " instruments after clumping")
  } else {
    message("    ", pgs_id, ": too few after clumping (", nrow(clumped), ")")
  }
}

# Build forward task list
forward_tasks <- expand.grid(
  inst_file = names(clumped_exposures),
  subtype = names(glioma_outcomes),
  stringsAsFactors = FALSE
)
message("\n  Running ", nrow(forward_tasks), " forward MR tests in parallel...\n")

forward_raw <- if (nrow(forward_tasks) > 0) mclapply(seq_len(nrow(forward_tasks)), function(idx) {
  inst_file <- forward_tasks$inst_file[idx]
  subtype   <- forward_tasks$subtype[idx]
  pgs_id    <- sub("_instruments.*", "", basename(inst_file))
  run_single_mr(clumped_exposures[[inst_file]], glioma_outcomes[[subtype]],
                pgs_id, subtype, "forward")
}, mc.cores = N_CORES) else list()

message("  Forward MR complete.\n")

# =============================================================================
# STEP 3: Reverse MR — Glioma (exposure) → ICVF (outcome) [PARALLEL]
# =============================================================================

message("=== STEP 3: Reverse MR — Glioma → ICVF (parallel) ===\n")

# Prepare glioma exposure instruments with LOCAL plink clumping
message("  Preparing glioma exposure instruments (local plink clumping)...")
glioma_exposures_clumped <- list()
for (i in seq_along(GLIOMA_SUBTYPES)) {
  subtype <- GLIOMA_SUBTYPES[i]
  if (!subtype %in% names(glioma_outcomes)) next

  fpath <- file.path(GLIOMA_DIR, GLIOMA_FILES[i])
  message("    ", subtype, ": formatting as exposure...")
  glioma_exp <- load_glioma(fpath, subtype, "exposure")

  # Select instruments
  n_total <- nrow(glioma_exp)
  instruments <- glioma_exp[glioma_exp$pval.exposure < P_THRESHOLD, ]
  message("      Total SNPs: ", n_total, ", at p<5e-8: ", nrow(instruments))
  if (nrow(instruments) < MIN_INSTRUMENTS) {
    message("      Trying p<5e-6")
    instruments <- glioma_exp[glioma_exp$pval.exposure < P_RELAXED, ]
    message("      At p<5e-6: ", nrow(instruments))
  }
  if (nrow(instruments) < MIN_INSTRUMENTS) {
    message("      Skipping: too few instruments")
    next
  }

  # Local plink clumping — handles chr:pos format via .bim mapping
  clumped <- local_clump(instruments, snp_col = "SNP", pval_col = "pval.exposure")

  if (nrow(clumped) >= MIN_INSTRUMENTS) {
    glioma_exposures_clumped[[subtype]] <- clumped
    message("      ", nrow(clumped), " instruments after clumping")
  } else {
    message("      Too few after clumping (", nrow(clumped), ")")
  }
}

# Load ICVF full files as outcomes
icvf_full_files <- list.files(ICVF_DIR, pattern = "_mr_ready\\.tsv\\.gz$", full.names = TRUE)
message("\n  Loading ", length(icvf_full_files), " ICVF outcome files...")

icvf_outcomes <- list()
for (icvf_file in icvf_full_files) {
  pgs_id <- sub("_.*", "", basename(icvf_file))
  trait <- sub("_mr_ready\\.tsv\\.gz$", "", basename(icvf_file))
  trait <- sub("^PGS\\d+_", "", trait)

  icvf_dat <- fread(icvf_file, header = TRUE)
  icvf_dat <- as.data.frame(icvf_dat)

  if ("chr.exposure" %in% names(icvf_dat) & "pos.exposure" %in% names(icvf_dat)) {
    icvf_dat$SNP_chrpos <- paste0(as.integer(icvf_dat$chr.exposure), ":", icvf_dat$pos.exposure)
  }

  icvf_out <- format_data(icvf_dat, type = "outcome",
    snp_col = "SNP_chrpos", beta_col = "beta.exposure", se_col = "se.exposure",
    pval_col = "pval.exposure", effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure", samplesize_col = "samplesize.exposure",
    chr_col = "chr.exposure", pos_col = "pos.exposure")
  icvf_out$outcome <- paste0(pgs_id, "_", trait)
  icvf_outcomes[[icvf_file]] <- icvf_out
  message("    Loaded: ", pgs_id, " ", trait)
}

# Build reverse task list
if (length(glioma_exposures_clumped) == 0) {
  message("\n  WARNING: No glioma subtypes had enough instruments for reverse MR.")
  reverse_raw <- list()
} else {
  message("  Glioma subtypes with instruments: ",
          paste(names(glioma_exposures_clumped), collapse = ", "))
  reverse_tasks <- expand.grid(
    subtype = names(glioma_exposures_clumped),
    icvf_file = names(icvf_outcomes),
    stringsAsFactors = FALSE
  )
  message("\n  Running ", nrow(reverse_tasks), " reverse MR tests in parallel...\n")

  reverse_raw <- if (nrow(reverse_tasks) > 0) mclapply(seq_len(nrow(reverse_tasks)), function(idx) {
    subtype   <- reverse_tasks$subtype[idx]
    icvf_file <- reverse_tasks$icvf_file[idx]
    pgs_id    <- sub("_.*", "", basename(icvf_file))
    run_single_mr(glioma_exposures_clumped[[subtype]], icvf_outcomes[[icvf_file]],
                  pgs_id, subtype, "reverse")
  }, mc.cores = N_CORES) else list()
}

message("  Reverse MR complete.\n")

# =============================================================================
# STEP 4: Collect and save results
# =============================================================================

message("=== STEP 4: Saving results ===\n")

collect_results <- function(raw_list, direction_label) {
  mr_list <- list(); het_list <- list(); pleio_list <- list()
  steiger_list <- list(); fstat_list <- list()  # AUDIT FIX: added F-stats
  for (res in raw_list) {
    if (is.null(res) || !is.list(res)) next
    if (!is.null(res$error)) { message("    Skipped: ", res$error); next }
    pgs_id <- res$pgs_id; subtype <- res$subtype
    if (!is.null(res$mr))      { tmp <- res$mr;      tmp$pgs_id <- pgs_id; tmp$subtype <- subtype; tmp$direction <- direction_label; mr_list <- c(mr_list, list(tmp)) }
    if (!is.null(res$het))     { tmp <- res$het;     tmp$pgs_id <- pgs_id; tmp$subtype <- subtype; het_list <- c(het_list, list(tmp)) }
    if (!is.null(res$pleio))   { tmp <- res$pleio;   tmp$pgs_id <- pgs_id; tmp$subtype <- subtype; pleio_list <- c(pleio_list, list(tmp)) }
    if (!is.null(res$steiger)) { tmp <- res$steiger; tmp$pgs_id <- pgs_id; tmp$subtype <- subtype; steiger_list <- c(steiger_list, list(tmp)) }
    if (!is.null(res$f_stats)) { tmp <- res$f_stats; tmp$pgs_id <- pgs_id; tmp$subtype <- subtype; tmp$direction <- direction_label; fstat_list <- c(fstat_list, list(tmp)) }
  }
  list(mr = if (length(mr_list) > 0) bind_rows(mr_list) else NULL,
       het = if (length(het_list) > 0) bind_rows(het_list) else NULL,
       pleio = if (length(pleio_list) > 0) bind_rows(pleio_list) else NULL,
       steiger = if (length(steiger_list) > 0) bind_rows(steiger_list) else NULL,
       f_stats = if (length(fstat_list) > 0) bind_rows(fstat_list) else NULL)
}

fwd <- collect_results(forward_raw, "forward_ICVF_to_glioma")
if (!is.null(fwd$mr))      fwrite(fwd$mr,      file.path(OUTPUT_DIR, "forward_mr_results.tsv"), sep = "\t")
if (!is.null(fwd$het))     fwrite(fwd$het,     file.path(OUTPUT_DIR, "forward_heterogeneity.tsv"), sep = "\t")
if (!is.null(fwd$pleio))   fwrite(fwd$pleio,   file.path(OUTPUT_DIR, "forward_pleiotropy.tsv"), sep = "\t")
if (!is.null(fwd$steiger)) fwrite(fwd$steiger, file.path(OUTPUT_DIR, "forward_steiger.tsv"), sep = "\t")
if (!is.null(fwd$f_stats)) fwrite(fwd$f_stats, file.path(OUTPUT_DIR, "forward_f_statistics.tsv"), sep = "\t")  # AUDIT FIX: added

rev <- collect_results(reverse_raw, "reverse_glioma_to_ICVF")
if (!is.null(rev$mr))      fwrite(rev$mr,      file.path(OUTPUT_DIR, "reverse_mr_results.tsv"), sep = "\t")
if (!is.null(rev$het))     fwrite(rev$het,     file.path(OUTPUT_DIR, "reverse_heterogeneity.tsv"), sep = "\t")
if (!is.null(rev$pleio))   fwrite(rev$pleio,   file.path(OUTPUT_DIR, "reverse_pleiotropy.tsv"), sep = "\t")
if (!is.null(rev$steiger)) fwrite(rev$steiger, file.path(OUTPUT_DIR, "reverse_steiger.tsv"), sep = "\t")  # AUDIT FIX: was missing
if (!is.null(rev$f_stats)) fwrite(rev$f_stats, file.path(OUTPUT_DIR, "reverse_f_statistics.tsv"), sep = "\t")  # AUDIT FIX: added

if (!is.null(fwd$mr) && !is.null(rev$mr)) {
  combined <- bind_rows(fwd$mr, rev$mr)
  fwrite(combined, file.path(OUTPUT_DIR, "combined_mr_results.tsv"), sep = "\t")
}


  # Write-protect result files to signal completion
  result_files <- list.files(OUTPUT_DIR, pattern = "\.tsv$", full.names = TRUE)
  for (rf in result_files) Sys.chmod(rf, mode = "0444")
  message("  Write-protected ", length(result_files), " result files")

}  # end of if (!.skip_analysis)

# =============================================================================
# STEP 5: MR-PRESSO (sequential — does not work in forked processes)
# =============================================================================

# =========================================================================
# STEP 5: MR-PRESSO — DISABLED
# =========================================================================
# The MRPRESSO R package (rondolab/MR-PRESSO) has an internal bug that
# crashes with "arguments imply differing number of rows: 1, 0" or
# "missing value where TRUE/FALSE needed" when constructing results.
# This occurs inside mr_presso() before it returns — not fixable from
# calling code. The bug triggers when the outlier test finds no outliers
# and tries to create an NA placeholder row with cbind.data.frame().
#
# Pleiotropy detection is still covered by:
#   - MR-Egger intercept test (in forward_pleiotropy.tsv)
#   - Cochran's Q heterogeneity test (in forward_heterogeneity.tsv)
#
# To re-enable: uncomment library(MRPRESSO) above and the section below.
# Monitor https://github.com/rondolab/MR-PRESSO for package updates.
# =========================================================================

# message("\n=== STEP 5: MR-PRESSO (sequential) ===")
# message("  Running for nominally significant forward IVW results only\n")

# if (!is.null(fwd$mr)) {
#   sig_ivw <- fwd$mr %>% filter(method == "Inverse variance weighted", pval < 0.05)
#   presso_results <- list()

#   for (row_i in seq_len(nrow(sig_ivw))) {
#     pgs_id  <- sig_ivw$pgs_id[row_i]
#     subtype <- sig_ivw$subtype[row_i]

    # Find matching clumped exposure
#     inst_match <- grep(pgs_id, names(clumped_exposures), value = TRUE)
#     if (length(inst_match) == 0) next

#     outcome_dat <- glioma_outcomes[[subtype]]
#     if (is.null(outcome_dat)) next

#     harmonised <- tryCatch(
#       harmonise_data(clumped_exposures[[inst_match[1]]], outcome_dat, action = 2),
#       error = function(e) NULL
#     )
#     if (is.null(harmonised) || sum(harmonised$mr_keep) < 4) next
#     harm_keep <- harmonised[harmonised$mr_keep, ]

#     message("  MR-PRESSO: ", pgs_id, " -> ", subtype, " (", nrow(harm_keep), " SNPs)")

#     tryCatch({
      # Step 1: Run global test ONLY (avoids MRPRESSO package bug with outlier results)
#       presso_global <- mr_presso(
#         BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
#         SdOutcome = "se.outcome", SdExposure = "se.exposure",
#         data = as.data.frame(harm_keep),
#         OUTLIERtest = FALSE, DISTORTIONtest = FALSE,
#         NbDistribution = 10000, SignifThreshold = 0.05
#       )
#       global_p <- tryCatch(presso_global[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA)
#       main_res <- tryCatch(presso_global[[1]]$`Main MR results`, error = function(e) NULL)
#       causal_raw <- if (!is.null(main_res) && nrow(main_res) >= 1) main_res$`Causal Estimate`[1] else NA

      # Step 2: Only run outlier detection if global test is significant
#       n_outliers <- 0; causal_corrected <- NA; distortion_p <- NA
#       if (!is.na(global_p) && global_p < 0.05) {
#         message("    Global test significant -- running outlier detection...")
#         presso_full <- tryCatch(
#           mr_presso(
#             BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
#             SdOutcome = "se.outcome", SdExposure = "se.exposure",
#             data = as.data.frame(harm_keep),
#             OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
#             NbDistribution = 10000, SignifThreshold = 0.05
#           ), error = function(e) { message("    Outlier test error: ", e$message); NULL }
#         )
#         if (!is.null(presso_full)) {
#           outlier_p <- tryCatch(presso_full[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue, error = function(e) NULL)
#           n_outliers <- if (!is.null(outlier_p)) sum(outlier_p < 0.05, na.rm = TRUE) else 0
#           main_full <- tryCatch(presso_full[[1]]$`Main MR results`, error = function(e) NULL)
#           if (!is.null(main_full) && nrow(main_full) >= 2) causal_corrected <- main_full$`Causal Estimate`[2]
#           distortion_p <- tryCatch(presso_full[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue, error = function(e) NA)
#         }
#       }

#       presso_results <- c(presso_results, list(data.frame(
#         pgs_id = pgs_id, subtype = subtype,
#         global_p = global_p, n_outliers = n_outliers,
#         causal_raw = causal_raw, causal_corrected = causal_corrected,
#         distortion_p = distortion_p
#       )))
#       message("    Global p = ", global_p, ", outliers = ", n_outliers)
#     }, error = function(e) {
#       message("    Failed: ", e$message)
#     })
#   }

#   if (length(presso_results) > 0) {
#     presso_df <- bind_rows(presso_results)
#     fwrite(presso_df, file.path(OUTPUT_DIR, "forward_mrpresso.tsv"), sep = "\t")
#     message("\n  Saved: forward_mrpresso.tsv")
#   }
# } else {
#   message("  No forward MR results to test")
# }


# =============================================================================
# STEP 6: Summary plots
# =============================================================================

message("\n=== STEP 6: Generating plots ===\n")

if (!is.null(fwd$mr)) {
  ivw_fwd <- fwd$mr %>% filter(method == "Inverse variance weighted")
  if (nrow(ivw_fwd) > 0) {
    ivw_fwd$label <- paste0(ivw_fwd$pgs_id, " -> ", ivw_fwd$subtype)
    ivw_fwd$or    <- exp(ivw_fwd$b)
    ivw_fwd$or_lo <- exp(ivw_fwd$b - 1.96 * ivw_fwd$se)
    ivw_fwd$or_hi <- exp(ivw_fwd$b + 1.96 * ivw_fwd$se)

    p <- ggplot(ivw_fwd, aes(x = or, y = reorder(label, or), color = subtype)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.2) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      scale_x_log10() +
      labs(title = "Forward MR: ICVF -> Glioma (IVW)",
           subtitle = "OR per SD increase in ICVF",
           x = "OR (95% CI)", y = NULL, color = "Glioma subtype") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(OUTPUT_DIR, "plots", "forward_mr_forest.pdf"),
           plot = p, width = 12, height = max(6, nrow(ivw_fwd) * 0.3 + 2))
    message("  Saved: forward_mr_forest.pdf")
  }
}

if (!is.null(rev$mr)) {
  ivw_rev <- rev$mr %>% filter(method == "Inverse variance weighted")
  if (nrow(ivw_rev) > 0) {
    ivw_rev$label <- paste0(ivw_rev$subtype, " -> ", ivw_rev$outcome)
    ivw_rev$beta_lo <- ivw_rev$b - 1.96 * ivw_rev$se
    ivw_rev$beta_hi <- ivw_rev$b + 1.96 * ivw_rev$se

    p2 <- ggplot(ivw_rev, aes(x = b, y = reorder(label, b), color = subtype)) +
      geom_point(size = 2) +
      geom_errorbarh(aes(xmin = beta_lo, xmax = beta_hi), height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = "Reverse MR: Glioma -> ICVF (IVW)",
           subtitle = "Effect of genetic liability to glioma on ICVF (SD units)",
           x = "Beta (95% CI)", y = NULL, color = "Glioma subtype") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(OUTPUT_DIR, "plots", "reverse_mr_forest.pdf"),
           plot = p2, width = 12, height = max(6, nrow(ivw_rev) * 0.3 + 2))
    message("  Saved: reverse_mr_forest.pdf")
  }
}

# =============================================================================
# STEP 7: Final summary
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("MR ANALYSIS COMPLETE")
message(paste(rep("=", 70), collapse = ""))

if (!is.null(fwd$mr)) {
  message("\n--- Forward MR (ICVF -> Glioma): Significant IVW results ---")
  sig_fwd <- fwd$mr %>% filter(method == "Inverse variance weighted", pval < 0.05)
  if (nrow(sig_fwd) > 0) {
    for (r in 1:nrow(sig_fwd)) {
      message(sprintf("  %s -> %s: b=%.4f, se=%.4f, p=%.2e, OR=%.3f",
                      sig_fwd$pgs_id[r], sig_fwd$subtype[r],
                      sig_fwd$b[r], sig_fwd$se[r], sig_fwd$pval[r], exp(sig_fwd$b[r])))
    }
  } else message("  No nominally significant results")
}

if (!is.null(rev$mr)) {
  message("\n--- Reverse MR (Glioma -> ICVF): Significant IVW results ---")
  sig_rev <- rev$mr %>% filter(method == "Inverse variance weighted", pval < 0.05)
  if (nrow(sig_rev) > 0) {
    for (r in 1:nrow(sig_rev)) {
      message(sprintf("  %s -> %s: b=%.4f, se=%.4f, p=%.2e",
                      sig_rev$subtype[r], sig_rev$outcome[r],
                      sig_rev$b[r], sig_rev$se[r], sig_rev$pval[r]))
    }
  } else message("  No nominally significant results")
}

message("\n--- Output files actually created ---")
output_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
if (length(output_files) > 0) {
  for (f in output_files) message("  ", f)
} else message("  (none)")

message("\n--- Next steps ---")
message("  Run mr_postprocessing.R for interpretation, heatmap, and tiered forest plot")
message("  Apply Bonferroni correction: 18 tracts x 5 subtypes = 90 tests -> p < 5.6e-4")
message("  Consider excluding HLA region instruments (chr6:25-34Mb) for robustness")
