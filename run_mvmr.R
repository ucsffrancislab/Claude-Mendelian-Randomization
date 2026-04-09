#!/usr/bin/env Rscript
# =============================================================================
# Multivariable MR: ICVF + FA + MD + ISOVF + OD → Glioma
# =============================================================================
# Tests whether ICVF's protective effect on glioma is independent of other
# diffusion metrics (FA, MD, ISOVF, OD) from the same white matter tracts.
#
# Requires: data.table, TwoSampleMR, ggplot2, dplyr, parallel
# =============================================================================

library(data.table)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(parallel)

# --- Command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
.output_dir_override <- NULL
if (length(args) >= 2) {
  for (i in seq_len(length(args) - 1)) {
    if (args[i] == "--output-dir") .output_dir_override <- args[i + 1]
  }
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# Data directories
ICVF_DIR <- "icvf_mr_ready"
MVMR_DIR <- "mvmr_mr_ready"
GLIOMA_DIR <- "../20260326-GWAS_summary_stats/20260330a-results"
MAPPING_FILE <- "mvmr_additional_metrics_mapping.csv"

# LD reference
LD_REF_PATH <- "/francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/EUR"
PLINK_BIN <- tryCatch(genetics.binaRies::get_plink_binary(),
  error = function(e) { p <- Sys.which("plink"); if (nchar(p) > 0) p else Sys.which("plink1.9") })

# Glioma subtypes
GLIOMA_SUBTYPES <- c("all_glioma", "idhwt", "idhmt", "idhmt_intact", "idhmt_codel")
GLIOMA_FILES <- c(
  "all_glioma/final/all_glioma_meta_summary_stats.tsv.gz",
  "idhwt/final/IDHwt_meta_summary_stats.tsv.gz",
  "idhmt/final/IDHmut_meta_summary_stats.tsv.gz",
  "idhmt_intact/final/IDHmut_1p19q_intact_meta_summary_stats.tsv.gz",
  "idhmt_codel/final/IDHmut_1p19q_codel_meta_summary_stats.tsv.gz"
)

# Output
OUTPUT_DIR <- if (!is.null(.output_dir_override)) .output_dir_override else "mvmr_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE)

# Parameters
P_THRESHOLD <- 5e-8
MIN_INSTRUMENTS <- 5  # MVMR needs more instruments than univariable
CLUMP_R2 <- 0.001
CLUMP_KB <- 10000
SAMPLE_SIZE <- 33224

# Parallelization
N_CORES <- min(detectCores() - 1, 60)
message("Using ", N_CORES, " cores")
message("Plink: ", PLINK_BIN)
message("LD ref: ", LD_REF_PATH)

# Metrics to include alongside ICVF
METRICS <- c("ICVF", "FA", "MD", "ISOVF", "OD")

# =============================================================================
# LOCAL LD CLUMPING (same as bidirectional MR script)
# =============================================================================

local_clump <- function(dat, snp_col = "SNP", pval_col = "pval.exposure") {
  bim <- fread(paste0(LD_REF_PATH, ".bim"), header = FALSE,
               col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
  bim$chrpos <- paste0(bim$chr, ":", bim$pos)

  dat_df <- as.data.frame(dat)
  original_snps <- dat_df[[snp_col]]

  n_rsid <- sum(grepl("^rs", original_snps))
  is_rsid <- n_rsid > length(original_snps) * 0.5

  if (!is_rsid) {
    chrpos <- sub("^([^:]+:[^:]+).*", "\\1", original_snps)
    lookup <- setNames(bim$rsid, bim$chrpos)
    mapped_rsids <- lookup[chrpos]
    has_rsid <- !is.na(mapped_rsids)
    message("    Mapped ", sum(has_rsid), "/", length(has_rsid), " SNPs to rsIDs for clumping")
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
    result[[snp_col]] <- result$original_snp
    result$original_snp <- NULL
  }
  return(result)
}

# =============================================================================
# STEP 1: Build tract-to-file mapping
# =============================================================================


# Tract name lookup: ICVF filename format → mapping CSV format
tract_name_lookup <- c(
  "body_corpus_callosum" = "Body_of_corpus_callosum",
  "cerebral_peduncle_(R)" = "Cerebral_peduncle_R",
  "cingulum_cingulate_gyrus_(L)" = "Cingulum_cingulate_gyrus_L",
  "cingulum_cingulate_gyrus_(R)" = "Cingulum_cingulate_gyrus_R",
  "cingulum_hippocampus_(L)" = "Cingulum_hippocampus_L",
  "cingulum_hippocampus_(R)" = "Cingulum_hippocampus_R",
  "genu_corpus_callosum" = "Genu_of_corpus_callosum",
  "middle_cerebellar_peduncle" = "Middle_cerebellar_peduncle",
  "posterior_limb_int_capsule_(L)" = "Posterior_limb_of_internal_capsule_L",
  "retrolenticular_int_capsule_(L)" = "Retrolenticular_part_of_internal_capsule_L",
  "retrolenticular_int_capsule_(R)" = "Retrolenticular_part_of_internal_capsule_R",
  "sagittal_stratum_(L)" = "Sagittal_stratum_L",
  "sagittal_stratum_(R)" = "Sagittal_stratum_R",
  "superior_corona_radiata_(L)" = "Superior_corona_radiata_L",
  "superior_corona_radiata_(R)" = "Superior_corona_radiata_R",
  "acoustic_radiation_(R)" = "ar_r",
  "forceps_major" = "fma",
  "parahippocampal_cingulum_(R)" = "cgh_r"
)

message("\n=== STEP 1: Building tract-to-file mapping ===\n")

# Load the additional metrics mapping
mvmr_map <- read.csv(MAPPING_FILE, stringsAsFactors = FALSE)

# Build a lookup: for each tract, which files contain which metrics
# ICVF files are in ICVF_DIR, other metrics in MVMR_DIR
tract_files <- list()

# Get unique tracts from the ICVF mapping
icvf_files <- list.files(ICVF_DIR, pattern = "_mr_ready\\.tsv\\.gz$", full.names = TRUE)

for (icvf_file in icvf_files) {
  # Parse PGS ID and tract name from filename
  pgs_id <- sub("_.*", "", basename(icvf_file))
  tract <- sub("_mr_ready\\.tsv\\.gz$", "", basename(icvf_file))
  tract <- sub("^PGS\\d+_", "", tract)

  # Map ICVF filename tract name to the mapping CSV tract name
  tract_clean <- sub("^ICVF_", "", tract)
  
  if (tract_clean %in% names(tract_name_lookup)) {
    tract_mapped <- tract_name_lookup[[tract_clean]]
  } else {
    tract_mapped <- tract_clean
  }

  map_matches <- mvmr_map[mvmr_map$tract == tract_mapped, ]

  if (nrow(map_matches) == 0) {
    message("  WARNING: No mapping found for tract: ", tract_clean, " (tried: ", tract_mapped, ")")
    next
  }

  # Build file paths for each metric
  files <- list(ICVF = icvf_file)
  method <- map_matches$method[1]

  for (metric in c("FA", "MD", "ISOVF", "OD")) {
    metric_row <- map_matches[map_matches$metric == metric, ]
    if (nrow(metric_row) > 0) {
      fname <- paste0(metric, "_", method, "_", metric_row$tract[1], "_mr_ready.tsv.gz")
      fpath <- file.path(MVMR_DIR, fname)
      if (file.exists(fpath)) {
        files[[metric]] <- fpath
      } else {
        message("  WARNING: Missing ", metric, " file: ", fpath)
      }
    }
  }

  tract_files[[tract_clean]] <- list(
    pgs_id = pgs_id,
    method = method,
    files = files
  )

  metrics_found <- paste(names(files), collapse = ", ")
  message("  ", tract_clean, " (", method, "): ", metrics_found)
}

message("\n  Tracts with all 5 metrics: ",
        sum(sapply(tract_files, function(x) length(x$files) == 5)),
        "/", length(tract_files))

# =============================================================================
# STEP 2: Load glioma outcome data
# =============================================================================

message("\n=== STEP 2: Loading glioma outcome data ===\n")

glioma_outcomes <- list()
for (i in seq_along(GLIOMA_SUBTYPES)) {
  fpath <- file.path(GLIOMA_DIR, GLIOMA_FILES[i])
  if (!file.exists(fpath)) next
  message("  Loading: ", GLIOMA_SUBTYPES[i])
  dat <- fread(fpath, header = TRUE)
  dat <- as.data.frame(dat)
  dat$N <- dat$N_CASES + dat$N_CONTROLS
  dat$SNP_chrpos <- paste0(as.integer(dat$CHR), ":", dat$BP)
  out <- format_data(dat, type = "outcome",
    snp_col = "SNP_chrpos", beta_col = "BETA", se_col = "SE", pval_col = "P",
    effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "A1_FREQ",
    samplesize_col = "N", chr_col = "CHR", pos_col = "BP")
  out$outcome <- GLIOMA_SUBTYPES[i]
  glioma_outcomes[[GLIOMA_SUBTYPES[i]]] <- out
}
message("  Loaded ", length(glioma_outcomes), " glioma subtypes\n")

# =============================================================================
# STEP 3: Run MVMR for each tract × subtype
# =============================================================================

message("=== STEP 3: Running MVMR ===\n")

run_mvmr_for_tract <- function(tract_name, tract_info, glioma_outcomes) {
  results <- list()
  files <- tract_info$files
  pgs_id <- tract_info$pgs_id

  # Only proceed if we have all 5 metrics
  if (length(files) < 5) {
    message("  ", tract_name, ": skipping (only ", length(files), "/5 metrics)")
    return(NULL)
  }

  message("  Processing: ", tract_name)

  # --- Load all exposure data and find instruments ---
  # For MVMR, instruments are SNPs associated with ANY of the exposures

  all_exposures <- list()
  all_instruments <- data.frame()

  for (metric in names(files)) {
    dat <- fread(files[[metric]], header = TRUE)
    dat <- as.data.frame(dat)

    # Create chr:pos SNP ID with integer chromosome
    if ("chr.exposure" %in% names(dat) & "pos.exposure" %in% names(dat)) {
      dat$SNP <- paste0(as.integer(dat$chr.exposure), ":", dat$pos.exposure)
    }

    dat$exposure <- paste0(metric, "_", tract_name)
    all_exposures[[metric]] <- dat

    # Collect genome-wide significant SNPs from this metric
    sig <- dat[dat$pval.exposure < P_THRESHOLD, ]
    if (nrow(sig) > 0) {
      all_instruments <- rbind(all_instruments,
        data.frame(SNP = sig$SNP, pval = sig$pval.exposure, source_metric = metric,
                   stringsAsFactors = FALSE))
    }
  }

  # Deduplicate — keep the most significant p-value per SNP
  if (nrow(all_instruments) == 0) {
    message("    No instruments at p<", P_THRESHOLD)
    return(NULL)
  }

  all_instruments <- all_instruments %>%
    group_by(SNP) %>%
    summarise(pval = min(pval), n_metrics = n(), .groups = "drop") %>%
    as.data.frame()

  message("    Combined instruments: ", nrow(all_instruments),
          " unique SNPs from ", length(unique(all_instruments$n_metrics)), " metric(s)")

  # --- LD clump the combined instrument set ---
  all_instruments$pval.exposure <- all_instruments$pval
  clumped <- local_clump(all_instruments, snp_col = "SNP", pval_col = "pval.exposure")

  if (nrow(clumped) < MIN_INSTRUMENTS) {
    message("    Too few instruments after clumping: ", nrow(clumped))
    return(NULL)
  }
  message("    After clumping: ", nrow(clumped), " instruments")

  clumped_snps <- clumped$SNP

  # --- For each glioma subtype, run MVMR ---
  for (subtype in names(glioma_outcomes)) {
    message("    vs ", subtype)

    outcome_dat <- glioma_outcomes[[subtype]]

    # Build the multi-exposure dataset
    # We need: for each clumped SNP, get beta/se from ALL 5 exposures + the outcome
    # TwoSampleMR's mv_multiple() expects a specific format

    # Approach: use mv_harmonise_data() which takes a list of exposure data
    # Each exposure needs the same set of SNPs

    exposure_list <- list()
    for (metric in METRICS) {
      exp_dat <- all_exposures[[metric]]
      # Filter to clumped SNPs
      exp_subset <- exp_dat[exp_dat$SNP %in% clumped_snps, ]

      if (nrow(exp_subset) == 0) {
        message("      No overlap for ", metric)
        next
      }

      # Ensure proper format
      exp_subset$id.exposure <- paste0(metric, "_", tract_name)
      exp_subset$exposure <- paste0(metric, "_", tract_name)
      exposure_list[[metric]] <- exp_subset
    }

    if (length(exposure_list) < 2) {
      message("      Need at least 2 exposures with data")
      next
    }

    # Combine all exposures
    mv_exposures <- do.call(rbind, exposure_list)

    # Run mv_harmonise_data
    mvdat <- tryCatch(
      mv_harmonise_data(mv_exposures, outcome_dat),
      error = function(e) { message("      Harmonisation failed: ", e$message); NULL }
    )

    if (is.null(mvdat)) next

    n_snps <- nrow(mvdat$exposure_beta)
    n_exposures <- ncol(mvdat$exposure_beta)
    message("      Harmonised: ", n_snps, " SNPs × ", n_exposures, " exposures")

    if (n_snps < MIN_INSTRUMENTS) {
      message("      Too few SNPs after harmonisation")
      next
    }

    # Run MVMR
    mvmr_res <- tryCatch(
      mv_multiple(mvdat),
      error = function(e) { message("      MVMR failed: ", e$message); NULL }
    )

    if (!is.null(mvmr_res) && !is.null(mvmr_res$result)) {
      res <- mvmr_res$result
      res$tract <- tract_name
      res$pgs_id <- pgs_id
      res$subtype <- subtype
      res$n_snps <- n_snps
      results <- c(results, list(res))

      # Print ICVF result
      icvf_row <- res[grepl("ICVF", res$exposure), ]
      if (nrow(icvf_row) > 0) {
        message(sprintf("      ICVF direct effect: b=%.4f, se=%.4f, p=%.2e",
                        icvf_row$b[1], icvf_row$se[1], icvf_row$pval[1]))
      }
    }
  }

  if (length(results) > 0) return(bind_rows(results))
  return(NULL)
}

# --- Run MVMR for each tract (sequential — each tract loads large files) ---

all_mvmr_results <- list()

for (tract_name in names(tract_files)) {
  res <- tryCatch(
    run_mvmr_for_tract(tract_name, tract_files[[tract_name]], glioma_outcomes),
    error = function(e) { message("  ERROR for ", tract_name, ": ", e$message); NULL }
  )
  if (!is.null(res)) {
    all_mvmr_results <- c(all_mvmr_results, list(res))
  }
}

# =============================================================================
# STEP 4: Save results
# =============================================================================

message("\n=== STEP 4: Saving results ===\n")

if (length(all_mvmr_results) > 0) {
  mvmr_df <- bind_rows(all_mvmr_results)
  fwrite(mvmr_df, file.path(OUTPUT_DIR, "mvmr_results.tsv"), sep = "\t")
  message("  Saved: mvmr_results.tsv (", nrow(mvmr_df), " rows)")

  # Extract ICVF-specific results
  icvf_results <- mvmr_df[grepl("ICVF", mvmr_df$exposure), ]
  fwrite(icvf_results, file.path(OUTPUT_DIR, "mvmr_icvf_direct_effects.tsv"), sep = "\t")
  message("  Saved: mvmr_icvf_direct_effects.tsv (", nrow(icvf_results), " rows)")
} else {
  message("  No MVMR results produced")
}

# =============================================================================
# STEP 5: Comparison plot — univariable vs multivariable ICVF effects
# =============================================================================

message("\n=== STEP 5: Generating plots ===\n")

if (length(all_mvmr_results) > 0 && nrow(icvf_results) > 0) {

  # --- Forest plot of ICVF direct effects ---
  icvf_results$label <- paste0(icvf_results$tract, " -> ", icvf_results$subtype)
  icvf_results$or <- exp(icvf_results$b)
  icvf_results$or_lo <- exp(icvf_results$b - 1.96 * icvf_results$se)
  icvf_results$or_hi <- exp(icvf_results$b + 1.96 * icvf_results$se)
  icvf_results$sig <- ifelse(icvf_results$pval < 0.05, "p<0.05", "NS")

  p_forest <- ggplot(icvf_results, aes(x = or, y = reorder(label, or), color = sig)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.25) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("p<0.05" = "#d62728", "NS" = "grey60")) +
    scale_x_log10() +
    labs(title = "MVMR: ICVF Direct Effect on Glioma",
         subtitle = "Adjusted for FA, MD, ISOVF, OD from same tract",
         x = "OR (95% CI, log scale)", y = NULL, color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "plots", "mvmr_icvf_forest.pdf"),
         plot = p_forest, width = 12, height = max(6, nrow(icvf_results) * 0.3 + 2))
  message("  Saved: mvmr_icvf_forest.pdf")

  # --- Heatmap of all metrics' direct effects ---
  mvmr_df$metric <- sub("_.*", "", mvmr_df$exposure)
  mvmr_df$or <- exp(mvmr_df$b)
  mvmr_df$star <- case_when(
    mvmr_df$pval < 0.05 / nrow(icvf_results) ~ "***",
    mvmr_df$pval < 0.01 ~ "**",
    mvmr_df$pval < 0.05 ~ "*",
    TRUE ~ ""
  )
  mvmr_df$label <- sprintf("%.2f%s", mvmr_df$or, mvmr_df$star)

  # Filter to all_glioma for a clean heatmap
  heatmap_data <- mvmr_df[mvmr_df$subtype == "all_glioma", ]

  if (nrow(heatmap_data) > 0) {
    p_heat <- ggplot(heatmap_data, aes(x = metric, y = reorder(tract, b), fill = b)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = label), size = 3) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "log(OR)") +
      labs(title = "MVMR: Direct Effects of Each Metric on All Glioma",
           subtitle = "OR shown. * p<0.05, ** p<0.01, *** Bonferroni",
           x = "Diffusion Metric", y = "White Matter Tract") +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank())

    ggsave(file.path(OUTPUT_DIR, "plots", "mvmr_all_metrics_heatmap.pdf"),
           plot = p_heat, width = 10, height = 10)
    ggsave(file.path(OUTPUT_DIR, "plots", "mvmr_all_metrics_heatmap.png"),
           plot = p_heat, width = 10, height = 10, dpi = 150)
    message("  Saved: mvmr_all_metrics_heatmap.pdf/png")
  }
}

# =============================================================================
# STEP 6: Summary
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("MULTIVARIABLE MR COMPLETE")
message(paste(rep("=", 70), collapse = ""))

if (length(all_mvmr_results) > 0) {
  message("\n--- ICVF Direct Effects (adjusted for FA, MD, ISOVF, OD) ---")
  sig_icvf <- icvf_results[icvf_results$pval < 0.05, ]
  if (nrow(sig_icvf) > 0) {
    for (r in seq_len(nrow(sig_icvf))) {
      message(sprintf("  %s -> %s: b=%.4f, se=%.4f, p=%.2e, OR=%.3f (nSNP=%d)",
                      sig_icvf$tract[r], sig_icvf$subtype[r],
                      sig_icvf$b[r], sig_icvf$se[r], sig_icvf$pval[r],
                      exp(sig_icvf$b[r]), sig_icvf$n_snps[r]))
    }
  } else {
    message("  No nominally significant ICVF direct effects")
  }

  message("\n--- Comparison with univariable MR ---")
  n_univar <- 52  # from our standard MR
  n_mvmr <- nrow(sig_icvf)
  message("  Univariable MR: 52/90 nominally significant")
  message("  MVMR (ICVF direct): ", n_mvmr, "/", nrow(icvf_results), " nominally significant")
  if (n_mvmr > 0) {
    message("  ICVF effect SURVIVES adjustment for other metrics")
    message("  -> Axonal density specifically, not general white matter health")
  } else {
    message("  ICVF effect ATTENUATED after adjustment")
    message("  -> Effect may be driven by shared genetic architecture with other metrics")
  }
}

message("\n--- Output files ---")
output_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
if (length(output_files) > 0) {
  for (f in output_files) message("  ", f)
} else {
  message("  (none)")
}
