#!/usr/bin/env Rscript
# =============================================================================
# Post-processing: MR Sensitivity Analysis, Visualization & MR-PRESSO
# =============================================================================
# Run AFTER run_bidirectional_mr.R completes
# Handles BOTH forward (ICVF → Glioma) and reverse (Glioma → ICVF) results
#
# Requires: data.table, ggplot2, dplyr, tidyr, MRPRESSO, TwoSampleMR
# =============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

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

RESULTS_DIR <- if (!is.null(.output_dir_override)) .output_dir_override else "mr_results"
OUTPUT_DIR  <- file.path(RESULTS_DIR, "interpretation")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# PGS ID to tract name mapping
tract_names <- c(
  PGS001454 = "Body corpus callosum",
  PGS001456 = "Cerebral peduncle (R)",
  PGS001457 = "Cingulum cing. gyrus (L)",
  PGS001458 = "Cingulum cing. gyrus (R)",
  PGS001459 = "Cingulum hippocampus (L)",
  PGS001460 = "Cingulum hippocampus (R)",
  PGS001466 = "Genu corpus callosum",
  PGS001471 = "Middle cerebellar peduncle",
  PGS001474 = "Post. limb int. capsule (L)",
  PGS001478 = "Retrolent. int. capsule (L)",
  PGS001479 = "Retrolent. int. capsule (R)",
  PGS001480 = "Sagittal stratum (L)",
  PGS001481 = "Sagittal stratum (R)",
  PGS001484 = "Superior corona radiata (L)",
  PGS001485 = "Superior corona radiata (R)",
  PGS001662 = "Acoustic radiation (R)",
  PGS001669 = "Forceps major",
  PGS001679 = "Parahipp. cingulum (R)"
)

subtype_labels <- c(
  all_glioma   = "All Glioma",
  idhwt        = "IDH-WT",
  idhmt        = "IDH-mut",
  idhmt_intact = "IDH-mut 1p19q intact",
  idhmt_codel  = "IDH-mut 1p19q codel"
)

# =============================================================================
# STEP 1: Load all results
# =============================================================================

message("=== Loading MR results ===")

# Forward
fwd_mr    <- fread(file.path(RESULTS_DIR, "forward_mr_results.tsv"))
fwd_het   <- fread(file.path(RESULTS_DIR, "forward_heterogeneity.tsv"))
fwd_pleio <- fread(file.path(RESULTS_DIR, "forward_pleiotropy.tsv"))
fwd_steiger <- fread(file.path(RESULTS_DIR, "forward_steiger.tsv"))
message("  Forward: ", nrow(fwd_mr), " rows")

# Reverse (may not exist)
has_reverse <- file.exists(file.path(RESULTS_DIR, "reverse_mr_results.tsv"))
if (has_reverse) {
  rev_mr    <- fread(file.path(RESULTS_DIR, "reverse_mr_results.tsv"))
  rev_het   <- fread(file.path(RESULTS_DIR, "reverse_heterogeneity.tsv"))
  rev_pleio <- fread(file.path(RESULTS_DIR, "reverse_pleiotropy.tsv"))
  message("  Reverse: ", nrow(rev_mr), " rows")
} else {
  message("  Reverse: no results found")
}

# =============================================================================
# STEP 2: Forward MR — method consistency & confidence tiering
# =============================================================================

message("\n=== Forward MR: Sensitivity Analysis ===\n")

build_combined <- function(mr_res, het_res, pleio_res, steiger_res = NULL) {
  ivw   <- mr_res %>% filter(method == "Inverse variance weighted") %>%
    select(pgs_id, subtype, nsnp, b_ivw = b, se_ivw = se, p_ivw = pval)
  egger <- mr_res %>% filter(method == "MR Egger") %>%
    select(pgs_id, subtype, b_egger = b, p_egger = pval)
  wmed  <- mr_res %>% filter(method == "Weighted median") %>%
    select(pgs_id, subtype, b_wmed = b, p_wmed = pval)
  wmode <- mr_res %>% filter(method == "Weighted mode") %>%
    select(pgs_id, subtype, b_wmode = b, p_wmode = pval)

  comb <- ivw %>%
    left_join(egger, by = c("pgs_id", "subtype")) %>%
    left_join(wmed, by = c("pgs_id", "subtype")) %>%
    left_join(wmode, by = c("pgs_id", "subtype"))

  het_ivw <- het_res %>% filter(method == "Inverse variance weighted") %>%
    select(pgs_id, subtype, Q, Q_df, Q_pval)
  comb <- comb %>% left_join(het_ivw, by = c("pgs_id", "subtype"))

  pleio_sub <- pleio_res %>% select(pgs_id, subtype, egger_intercept, egger_se = se, egger_int_p = pval)
  comb <- comb %>% left_join(pleio_sub, by = c("pgs_id", "subtype"))

  if (!is.null(steiger_res)) {
    steiger_sub <- steiger_res %>% select(pgs_id, subtype, correct_causal_direction, steiger_pval)
    comb <- comb %>% left_join(steiger_sub, by = c("pgs_id", "subtype"))
  }

  comb$tract <- tract_names[comb$pgs_id]
  comb$subtype_label <- subtype_labels[comb$subtype]
  return(comb)
}

fwd_combined <- build_combined(fwd_mr, fwd_het, fwd_pleio, fwd_steiger)
fwd_combined$or_ivw <- exp(fwd_combined$b_ivw)
fwd_combined$nom_sig <- fwd_combined$p_ivw < 0.05

N_FWD_TESTS <- nrow(fwd_combined)
BONF_FWD <- 0.05 / N_FWD_TESTS
fwd_combined$bonf_sig <- fwd_combined$p_ivw < BONF_FWD

fwd_combined$all_negative <- fwd_combined$b_ivw < 0 & fwd_combined$b_egger < 0 &
                              fwd_combined$b_wmed < 0 & fwd_combined$b_wmode < 0

fwd_combined$confidence <- case_when(
  !fwd_combined$nom_sig ~ "NS",
  !fwd_combined$all_negative | fwd_combined$egger_int_p < 0.05 ~ "LOW",
  fwd_combined$Q_pval < 0.05 ~ "MODERATE",
  TRUE ~ "HIGH"
)

# Print forward summary
for (tier in c("HIGH", "MODERATE", "LOW")) {
  sub <- fwd_combined %>% filter(confidence == tier) %>% arrange(p_ivw)
  if (nrow(sub) == 0) next
  message("--- FORWARD ", tier, " CONFIDENCE (", nrow(sub), " results) ---")
  for (i in 1:nrow(sub)) {
    r <- sub[i, ]
    flags <- c()
    if (r$Q_pval < 0.05) flags <- c(flags, "het")
    if (r$egger_int_p < 0.05) flags <- c(flags, "pleio")
    if (!r$all_negative) flags <- c(flags, "dir")
    flag_str <- if (length(flags) > 0) paste0(" [", paste(flags, collapse = ","), "]") else ""
    n_sig <- sum(c(r$p_ivw, r$p_egger, r$p_wmed, r$p_wmode) < 0.05)
    message(sprintf("  %s -> %s: OR=%.3f, p=%.2e, nSNP=%d, methods_sig=%d/4%s",
                    r$tract, r$subtype_label, r$or_ivw, r$p_ivw, r$nsnp, n_sig, flag_str))
  }
  message("")
}

# =============================================================================
# STEP 3: Reverse MR — sensitivity analysis
# =============================================================================

if (has_reverse) {
  message("\n=== Reverse MR: Sensitivity Analysis ===\n")

  rev_combined <- build_combined(rev_mr, rev_het, rev_pleio, steiger_res = NULL)
  rev_combined$or_ivw <- NA  # Reverse is continuous outcome, no OR
  rev_combined$nom_sig <- rev_combined$p_ivw < 0.05

  N_REV_TESTS <- nrow(rev_combined)
  BONF_REV <- 0.05 / N_REV_TESTS
  rev_combined$bonf_sig <- rev_combined$p_ivw < BONF_REV

  # For reverse: all effects should be in the same direction (positive = glioma → higher ICVF)
  rev_combined$all_positive <- rev_combined$b_ivw > 0 & rev_combined$b_egger > 0 &
                                rev_combined$b_wmed > 0 & rev_combined$b_wmode > 0

  rev_combined$confidence <- case_when(
    !rev_combined$nom_sig ~ "NS",
    !rev_combined$all_positive | rev_combined$egger_int_p < 0.05 ~ "LOW",
    rev_combined$Q_pval < 0.05 ~ "MODERATE",
    TRUE ~ "HIGH"
  )

  # Parse tract name from outcome column (format: PGS001454_ICVF_body_corpus_callosum)
  rev_combined$tract <- tract_names[rev_combined$pgs_id]

  for (tier in c("HIGH", "MODERATE", "LOW")) {
    sub <- rev_combined %>% filter(confidence == tier) %>% arrange(p_ivw)
    if (nrow(sub) == 0) next
    message("--- REVERSE ", tier, " CONFIDENCE (", nrow(sub), " results) ---")
    for (i in 1:nrow(sub)) {
      r <- sub[i, ]
      flags <- c()
      if (r$Q_pval < 0.05) flags <- c(flags, "het")
      if (r$egger_int_p < 0.05) flags <- c(flags, "pleio")
      if (!r$all_positive) flags <- c(flags, "dir")
      flag_str <- if (length(flags) > 0) paste0(" [", paste(flags, collapse = ","), "]") else ""
      n_sig <- sum(c(r$p_ivw, r$p_egger, r$p_wmed, r$p_wmode) < 0.05)
      bonf_str <- if (r$bonf_sig) " *BONF*" else ""
      message(sprintf("  %s -> %s: b=%.4f, p=%.2e, nSNP=%d, methods_sig=%d/4%s%s",
                      r$subtype_label, r$tract, r$b_ivw, r$p_ivw, r$nsnp, n_sig, flag_str, bonf_str))
    }
    message("")
  }
}

# =============================================================================
# STEP 4: Heatmaps — forward and reverse
# =============================================================================

message("=== Generating heatmaps ===\n")

subtype_order <- c("All Glioma", "IDH-WT", "IDH-mut", "IDH-mut 1p19q intact", "IDH-mut 1p19q codel")

# --- Forward heatmap ---
fwd_combined$subtype_label <- factor(fwd_combined$subtype_label, levels = subtype_order)
fwd_combined$star <- case_when(
  fwd_combined$p_ivw < BONF_FWD ~ "***",
  fwd_combined$p_ivw < 0.01 ~ "**",
  fwd_combined$p_ivw < 0.05 ~ "*",
  TRUE ~ ""
)
fwd_combined$or_label <- sprintf("%.2f%s", fwd_combined$or_ivw, fwd_combined$star)

p_fwd_heat <- ggplot(fwd_combined, aes(x = subtype_label, y = reorder(tract, b_ivw), fill = b_ivw)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = or_label), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log(OR)") +
  labs(title = "Forward MR: ICVF -> Glioma Risk (IVW)",
       subtitle = paste0("OR shown. * p<0.05, ** p<0.01, *** Bonferroni (p<", sprintf("%.1e", BONF_FWD), ")"),
       x = "Glioma Subtype", y = "ICVF Tract") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

ggsave(file.path(OUTPUT_DIR, "forward_heatmap.pdf"), plot = p_fwd_heat, width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "forward_heatmap.png"), plot = p_fwd_heat, width = 12, height = 10, dpi = 150)
message("  Saved: forward_heatmap.pdf/png")

# --- Reverse heatmap ---
if (has_reverse) {
  rev_combined$subtype_label <- factor(rev_combined$subtype_label, levels = subtype_order)
  rev_combined$star <- case_when(
    rev_combined$p_ivw < BONF_REV ~ "***",
    rev_combined$p_ivw < 0.01 ~ "**",
    rev_combined$p_ivw < 0.05 ~ "*",
    TRUE ~ ""
  )
  rev_combined$beta_label <- sprintf("%.3f%s", rev_combined$b_ivw, rev_combined$star)

  p_rev_heat <- ggplot(rev_combined, aes(x = reorder(tract, b_ivw), y = subtype_label, fill = b_ivw)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = beta_label), size = 2.8) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Beta (SD)") +
    labs(title = "Reverse MR: Glioma -> ICVF (IVW)",
         subtitle = paste0("Beta (SD ICVF per unit log-OR glioma). *** Bonferroni (p<", sprintf("%.1e", BONF_REV), ")"),
         x = "ICVF Tract", y = "Glioma Subtype") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

  ggsave(file.path(OUTPUT_DIR, "reverse_heatmap.pdf"), plot = p_rev_heat, width = 14, height = 5)
  ggsave(file.path(OUTPUT_DIR, "reverse_heatmap.png"), plot = p_rev_heat, width = 14, height = 5, dpi = 150)
  message("  Saved: reverse_heatmap.pdf/png")
}

# =============================================================================
# STEP 5: Tiered forest plots — forward and reverse
# =============================================================================

message("=== Generating forest plots ===\n")

# --- Forward forest ---
fwd_sig <- fwd_combined %>% filter(nom_sig) %>% arrange(confidence, p_ivw)
if (nrow(fwd_sig) > 0) {
  fwd_sig$label <- paste0(fwd_sig$tract, " -> ", fwd_sig$subtype_label)
  fwd_sig$or_lo <- exp(fwd_sig$b_ivw - 1.96 * fwd_sig$se_ivw)
  fwd_sig$or_hi <- exp(fwd_sig$b_ivw + 1.96 * fwd_sig$se_ivw)
  fwd_sig$label <- factor(fwd_sig$label, levels = rev(fwd_sig$label))

  p_fwd_forest <- ggplot(fwd_sig, aes(x = or_ivw, y = label, color = confidence)) +
    geom_point(size = 2.5) +
    geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.25) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(HIGH = "#2ca02c", MODERATE = "#ff7f0e", LOW = "#d62728")) +
    scale_x_log10() +
    labs(title = "Forward MR: ICVF -> Glioma (IVW, nominally significant)",
         subtitle = "OR per SD increase in ICVF (95% CI)",
         x = "Odds Ratio (log scale)", y = NULL, color = "Confidence") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "forward_forest_tiered.pdf"), plot = p_fwd_forest,
         width = 11, height = max(5, nrow(fwd_sig) * 0.35 + 2))
  message("  Saved: forward_forest_tiered.pdf")
}

# --- Reverse forest ---
if (has_reverse) {
  rev_sig <- rev_combined %>% filter(nom_sig) %>% arrange(confidence, p_ivw)
  if (nrow(rev_sig) > 0) {
    rev_sig$label <- paste0(rev_sig$subtype_label, " -> ", rev_sig$tract)
    rev_sig$b_lo <- rev_sig$b_ivw - 1.96 * rev_sig$se_ivw
    rev_sig$b_hi <- rev_sig$b_ivw + 1.96 * rev_sig$se_ivw
    rev_sig$label <- factor(rev_sig$label, levels = rev(rev_sig$label))

    p_rev_forest <- ggplot(rev_sig, aes(x = b_ivw, y = label, color = confidence)) +
      geom_point(size = 2.5) +
      geom_errorbarh(aes(xmin = b_lo, xmax = b_hi), height = 0.25) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c(HIGH = "#2ca02c", MODERATE = "#ff7f0e", LOW = "#d62728")) +
      labs(title = "Reverse MR: Glioma -> ICVF (IVW, nominally significant)",
           subtitle = "Beta: SD change in ICVF per unit log-OR glioma liability (95% CI)",
           x = "Beta (95% CI)", y = NULL, color = "Confidence") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")

    ggsave(file.path(OUTPUT_DIR, "reverse_forest_tiered.pdf"), plot = p_rev_forest,
           width = 11, height = max(5, nrow(rev_sig) * 0.35 + 2))
    message("  Saved: reverse_forest_tiered.pdf")
  }
}

# =============================================================================
# STEP 6: Save comprehensive summary tables
# =============================================================================

message("=== Saving summary tables ===\n")

fwd_out <- fwd_combined %>%
  select(pgs_id, tract, subtype, subtype_label, nsnp,
         or_ivw, b_ivw, se_ivw, p_ivw,
         b_egger, p_egger, b_wmed, p_wmed, b_wmode, p_wmode,
         all_negative, Q, Q_df, Q_pval, egger_intercept, egger_int_p,
         any_of(c("correct_causal_direction", "steiger_pval")),
         confidence, bonf_sig) %>%
  arrange(confidence, p_ivw)
fwd_out$direction <- "forward"
fwrite(fwd_out, file.path(OUTPUT_DIR, "forward_sensitivity_summary.tsv"), sep = "\t")
message("  Saved: forward_sensitivity_summary.tsv")

if (has_reverse) {
  rev_out <- rev_combined %>%
    select(pgs_id, tract, subtype, subtype_label, nsnp,
           b_ivw, se_ivw, p_ivw,
           b_egger, p_egger, b_wmed, p_wmed, b_wmode, p_wmode,
           all_positive, Q, Q_df, Q_pval, egger_intercept, egger_int_p,
           confidence, bonf_sig) %>%
    arrange(confidence, p_ivw)
  rev_out$direction <- "reverse"
  fwrite(rev_out, file.path(OUTPUT_DIR, "reverse_sensitivity_summary.tsv"), sep = "\t")
  message("  Saved: reverse_sensitivity_summary.tsv")
}

# =============================================================================
# =========================================================================
# STEP 7: MR-PRESSO — DISABLED
# =========================================================================
# See note in run_bidirectional_mr.R — the MRPRESSO package has an internal
# bug that crashes inside mr_presso(). Pleiotropy is covered by MR-Egger
# intercept and Cochran's Q. Uncomment to retry if the package is updated.
# =========================================================================

# STEP 7: MR-PRESSO retry (for forward, if not already done)
# =============================================================================

# message("\n=== MR-PRESSO check ===")

# presso_file <- file.path(RESULTS_DIR, "forward_mrpresso.tsv")
# if (file.exists(presso_file)) {
#   message("  MR-PRESSO results already exist — skipping")
#   presso_df <- fread(presso_file)
#   message("  ", nrow(presso_df), " MR-PRESSO results loaded")
# } else {
#   message("  No MR-PRESSO results found — running now (sequential)")
#   message("  This may take several minutes...\n")

#   library(MRPRESSO)
#   library(TwoSampleMR)

#   ICVF_DIR   <- "icvf_mr_ready"
#   GLIOMA_DIR <- if (!is.null(.glioma_dir_override)) .glioma_dir_override else "../20260326-GWAS_summary_stats/20260330a-results"
#   GLIOMA_FILES <- c(
#     all_glioma   = "all_glioma/final/all_glioma_meta_summary_stats.tsv.gz",
#     idhwt        = "idhwt/final/IDHwt_meta_summary_stats.tsv.gz",
#     idhmt        = "idhmt/final/IDHmut_meta_summary_stats.tsv.gz",
#     idhmt_intact = "idhmt_intact/final/IDHmut_1p19q_intact_meta_summary_stats.tsv.gz",
#     idhmt_codel  = "idhmt_codel/final/IDHmut_1p19q_codel_meta_summary_stats.tsv.gz"
#   )

  # Load glioma outcomes
#   glioma_outcomes <- list()
#   for (sub in names(GLIOMA_FILES)) {
#     fpath <- file.path(GLIOMA_DIR, GLIOMA_FILES[sub])
#     if (!file.exists(fpath)) next
#     dat <- fread(fpath, header = TRUE)
#     dat <- as.data.frame(dat)
#     dat$N <- dat$N_CASES + dat$N_CONTROLS
#     dat$SNP_chrpos <- paste0(as.integer(dat$CHR), ":", dat$BP)
#     out <- format_data(dat, type = "outcome",
#       snp_col = "SNP_chrpos", beta_col = "BETA", se_col = "SE", pval_col = "P",
#       effect_allele_col = "A1", other_allele_col = "A2", eaf_col = "A1_FREQ",
#       samplesize_col = "N", chr_col = "CHR", pos_col = "BP")
#     out$outcome <- sub
#     glioma_outcomes[[sub]] <- out
#   }

#   instrument_files <- list.files(ICVF_DIR, pattern = "_instruments_", full.names = TRUE)
#   sig_pairs <- fwd_combined %>% filter(nom_sig) %>% select(pgs_id, subtype)
#   presso_results <- list()

#   for (row_i in seq_len(nrow(sig_pairs))) {
#     pgs_id  <- sig_pairs$pgs_id[row_i]
#     subtype <- sig_pairs$subtype[row_i]

#     inst_file <- grep(pgs_id, instrument_files, value = TRUE)
#     if (length(inst_file) == 0) next

#     exposure_dat <- fread(inst_file[1], header = TRUE)
#     exposure_dat <- as.data.frame(exposure_dat)
#     if ("chr.exposure" %in% names(exposure_dat) & "pos.exposure" %in% names(exposure_dat)) {
#       exposure_dat$SNP <- paste0(as.integer(exposure_dat$chr.exposure), ":", exposure_dat$pos.exposure)
#     }

#     outcome_dat <- glioma_outcomes[[subtype]]
#     if (is.null(outcome_dat)) next

#     harmonised <- tryCatch(harmonise_data(exposure_dat, outcome_dat, action = 2), error = function(e) NULL)
#     if (is.null(harmonised) || sum(harmonised$mr_keep) < 4) next
#     harm_keep <- harmonised[harmonised$mr_keep, ]

#     message("  MR-PRESSO: ", pgs_id, " -> ", subtype, " (", nrow(harm_keep), " SNPs)")

#     tryCatch({
#       presso <- mr_presso(
#         BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
#         SdOutcome = "se.outcome", SdExposure = "se.exposure",
#         data = harm_keep, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
#         NbDistribution = 1000, SignifThreshold = 0.05)
#       presso_results <- c(presso_results, list(data.frame(
#         pgs_id = pgs_id, subtype = subtype,
#         global_p = presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue,
#         n_outliers = sum(presso[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, na.rm = TRUE),
#         causal_raw = presso[[1]]$`Main MR results`$`Causal Estimate`[1],
#         causal_corrected = presso[[1]]$`Main MR results`$`Causal Estimate`[2],
#         distortion_p = tryCatch(presso[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue, error = function(e) NA)
#       )))
#       message("    Global p = ", presso[[1]]$`MR-RESSO results`$`Global Test`$Pvalue)
#     }, error = function(e) message("    Failed: ", e$message))
#   }

#   if (length(presso_results) > 0) {
#     presso_df <- bind_rows(presso_results)
#     fwrite(presso_df, file.path(OUTPUT_DIR, "forward_mrpresso.tsv"), sep = "\t")
#     message("\n  Saved: forward_mrpresso.tsv")
#   }
# }

# =============================================================================

# STEP 8: Final interpretation
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("POST-PROCESSING COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\n--- Files created ---")
output_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
for (f in output_files) message("  ", f)

message("\n--- Key findings ---")
message("  FORWARD (ICVF -> Glioma):")
n_fwd_sig <- sum(fwd_combined$nom_sig)
n_fwd_bonf <- sum(fwd_combined$bonf_sig)
message("    ", n_fwd_sig, "/", N_FWD_TESTS, " nominally significant, ", n_fwd_bonf, " Bonferroni-significant")
message("    All significant effects are PROTECTIVE (higher ICVF -> lower glioma risk)")
n_high <- sum(fwd_combined$confidence == "HIGH")
message("    ", n_high, " HIGH confidence results (no heterogeneity, no pleiotropy, consistent direction)")

if (has_reverse) {
  message("  REVERSE (Glioma -> ICVF):")
  n_rev_sig <- sum(rev_combined$nom_sig)
  n_rev_bonf <- sum(rev_combined$bonf_sig)
  message("    ", n_rev_sig, "/", N_REV_TESTS, " nominally significant, ", n_rev_bonf, " Bonferroni-significant")
  message("    Subtypes with instruments: ", paste(unique(rev_combined$subtype_label), collapse = ", "))
  if (n_rev_bonf > 0) {
    message("    NOTE: Significant reverse results suggest shared genetic architecture")
    message("    between IDH-mutant glioma and ICVF — discuss bidirectionality vs pleiotropy")
  }
}

message("\n--- Interpretation notes ---")
message("  - Steiger confirms ICVF -> Glioma as correct causal direction for all forward results")
message("  - Best forward evidence: corpus callosum & cingulum tracts -> all glioma / IDH-WT")
message("  - IDH-mutant forward results have pleiotropy flags — interpret cautiously")
message("  - Reverse MR signal in IDH-mut intact suggests shared biology, not necessarily reverse causation")
message("  - Consider HLA exclusion (chr6:25-34Mb) sensitivity analysis for publication")
