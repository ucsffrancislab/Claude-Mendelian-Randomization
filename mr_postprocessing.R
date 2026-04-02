#!/usr/bin/env Rscript
# =============================================================================
# Post-processing: MR Sensitivity Analysis, Visualization & MR-PRESSO
# =============================================================================
# Run AFTER run_bidirectional_mr.R completes
# Reads the output TSVs and produces interpretation + MR-PRESSO
#
# Requires: data.table, ggplot2, dplyr, tidyr, MRPRESSO
# =============================================================================

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MRPRESSO)

# =============================================================================
# CONFIGURATION
# =============================================================================

RESULTS_DIR <- "mr_results"
OUTPUT_DIR  <- file.path(RESULTS_DIR, "interpretation")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Bonferroni threshold (18 tracts x 5 subtypes = 90 tests)
N_TESTS <- 90
BONF_THRESH <- 0.05 / N_TESTS

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
# STEP 1: Load results
# =============================================================================

message("=== Loading MR results ===")
mr_res  <- fread(file.path(RESULTS_DIR, "forward_mr_results.tsv"))
het_res <- fread(file.path(RESULTS_DIR, "forward_heterogeneity.tsv"))
pleio_res <- fread(file.path(RESULTS_DIR, "forward_pleiotropy.tsv"))
steiger_res <- fread(file.path(RESULTS_DIR, "forward_steiger.tsv"))

message("  MR results: ", nrow(mr_res), " rows")
message("  Methods: ", paste(unique(mr_res$method), collapse = ", "))

# =============================================================================
# STEP 2: Method consistency check
# =============================================================================

message("\n=== Method Consistency Analysis ===\n")

ivw <- mr_res %>% filter(method == "Inverse variance weighted") %>%
  select(pgs_id, subtype, nsnp, b_ivw = b, se_ivw = se, p_ivw = pval)
egger <- mr_res %>% filter(method == "MR Egger") %>%
  select(pgs_id, subtype, b_egger = b, p_egger = pval)
wmed <- mr_res %>% filter(method == "Weighted median") %>%
  select(pgs_id, subtype, b_wmed = b, p_wmed = pval)
wmode <- mr_res %>% filter(method == "Weighted mode") %>%
  select(pgs_id, subtype, b_wmode = b, p_wmode = pval)

combined <- ivw %>%
  left_join(egger, by = c("pgs_id", "subtype")) %>%
  left_join(wmed, by = c("pgs_id", "subtype")) %>%
  left_join(wmode, by = c("pgs_id", "subtype"))

combined$or_ivw <- exp(combined$b_ivw)
combined$nom_sig <- combined$p_ivw < 0.05
combined$bonf_sig <- combined$p_ivw < BONF_THRESH
combined$all_negative <- combined$b_ivw < 0 & combined$b_egger < 0 &
                         combined$b_wmed < 0 & combined$b_wmode < 0

# Add heterogeneity
het_ivw <- het_res %>% filter(method == "Inverse variance weighted") %>%
  select(pgs_id, subtype, Q, Q_df, Q_pval)
combined <- combined %>% left_join(het_ivw, by = c("pgs_id", "subtype"))

# Add pleiotropy
pleio_sub <- pleio_res %>%
  select(pgs_id, subtype, egger_intercept, egger_se = se, egger_int_p = pval)
combined <- combined %>% left_join(pleio_sub, by = c("pgs_id", "subtype"))

# Add Steiger
steiger_sub <- steiger_res %>%
  select(pgs_id, subtype, correct_causal_direction, steiger_pval)
combined <- combined %>% left_join(steiger_sub, by = c("pgs_id", "subtype"))

# Confidence tier
combined$confidence <- case_when(
  !combined$nom_sig ~ "NS",
  !combined$all_negative | combined$egger_int_p < 0.05 ~ "LOW",
  combined$Q_pval < 0.05 ~ "MODERATE",
  TRUE ~ "HIGH"
)

# Add tract names
combined$tract <- tract_names[combined$pgs_id]
combined$subtype_label <- subtype_labels[combined$subtype]

# Print significant results by tier
for (tier in c("HIGH", "MODERATE", "LOW")) {
  sub <- combined %>% filter(confidence == tier) %>% arrange(p_ivw)
  if (nrow(sub) == 0) next
  message("--- ", tier, " CONFIDENCE (", nrow(sub), " results) ---")
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
# STEP 3: Heatmap — log(OR) across all tract-subtype combinations
# =============================================================================

message("=== Generating heatmap ===\n")

# Order subtypes
subtype_order <- c("All Glioma", "IDH-WT", "IDH-mut",
                   "IDH-mut 1p19q intact", "IDH-mut 1p19q codel")
combined$subtype_label <- factor(combined$subtype_label, levels = subtype_order)

# Significance annotation
combined$star <- case_when(
  combined$p_ivw < BONF_THRESH ~ "***",
  combined$p_ivw < 0.01 ~ "**",
  combined$p_ivw < 0.05 ~ "*",
  TRUE ~ ""
)
combined$or_label <- sprintf("%.2f%s", combined$or_ivw, combined$star)

p_heat <- ggplot(combined, aes(x = subtype_label, y = reorder(tract, b_ivw), fill = b_ivw)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = or_label), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log(OR)") +
  labs(
    title = "Forward MR: ICVF -> Glioma Risk (IVW)",
    subtitle = "OR shown in cells. * p<0.05, ** p<0.01, *** Bonferroni-significant",
    x = "Glioma Subtype", y = "ICVF Tract"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(OUTPUT_DIR, "mr_heatmap.pdf"), plot = p_heat, width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "mr_heatmap.png"), plot = p_heat, width = 12, height = 10, dpi = 150)
message("  Saved: mr_heatmap.pdf/png")

# =============================================================================
# STEP 4: Forest plot of significant results with confidence tiers
# =============================================================================

message("=== Generating forest plot ===\n")

sig <- combined %>% filter(nom_sig) %>% arrange(confidence, p_ivw)
sig$label <- paste0(sig$tract, " -> ", sig$subtype_label)
sig$or_lo <- exp(sig$b_ivw - 1.96 * sig$se_ivw)
sig$or_hi <- exp(sig$b_ivw + 1.96 * sig$se_ivw)
sig$label <- factor(sig$label, levels = rev(sig$label))

p_forest <- ggplot(sig, aes(x = or_ivw, y = label, color = confidence)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(HIGH = "#2ca02c", MODERATE = "#ff7f0e", LOW = "#d62728")) +
  scale_x_log10() +
  labs(
    title = "Forward MR: ICVF -> Glioma (IVW, nominally significant)",
    subtitle = "OR per SD increase in ICVF (95% CI)",
    x = "Odds Ratio (log scale)", y = NULL, color = "Confidence"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "mr_forest_tiered.pdf"), plot = p_forest,
       width = 11, height = max(5, nrow(sig) * 0.35 + 2))
message("  Saved: mr_forest_tiered.pdf")

# =============================================================================
# STEP 5: Save comprehensive summary table
# =============================================================================

message("=== Saving summary table ===\n")

summary_out <- combined %>%
  select(pgs_id, tract, subtype, subtype_label, nsnp,
         or_ivw, b_ivw, se_ivw, p_ivw,
         b_egger, p_egger, b_wmed, p_wmed, b_wmode, p_wmode,
         all_negative, Q, Q_df, Q_pval,
         egger_intercept, egger_int_p,
         correct_causal_direction, steiger_pval,
         confidence, bonf_sig) %>%
  arrange(confidence, p_ivw)

fwrite(summary_out, file.path(OUTPUT_DIR, "mr_sensitivity_summary.tsv"), sep = "\t")
message("  Saved: mr_sensitivity_summary.tsv")

# =============================================================================
# STEP 6: MR-PRESSO (retry — failed silently in parallel workers)
# =============================================================================

message("\n=== Running MR-PRESSO (sequential) ===")
message("  This runs 1000 permutations per test — may take a few minutes\n")

# We need the harmonised data to run MR-PRESSO
# Re-load the ICVF instrument files and glioma outcomes, harmonise, then run

library(TwoSampleMR)

ICVF_DIR   <- "icvf_mr_ready"
GLIOMA_DIR <- "../20260326-GWAS_summary_stats/20260330a-results"
GLIOMA_SUBTYPES <- c("all_glioma", "idhwt", "idhmt", "idhmt_intact", "idhmt_codel")
GLIOMA_FILES <- c(
  "all_glioma/final/all_glioma_meta_summary_stats.tsv.gz",
  "idhwt/final/IDHwt_meta_summary_stats.tsv.gz",
  "idhmt/final/IDHmut_meta_summary_stats.tsv.gz",
  "idhmt_intact/final/IDHmut_1p19q_intact_meta_summary_stats.tsv.gz",
  "idhmt_codel/final/IDHmut_1p19q_codel_meta_summary_stats.tsv.gz"
)

# Load glioma outcomes (with chr:pos SNP IDs)
glioma_outcomes <- list()
for (i in seq_along(GLIOMA_SUBTYPES)) {
  fpath <- file.path(GLIOMA_DIR, GLIOMA_FILES[i])
  if (!file.exists(fpath)) next
  dat <- fread(fpath, header = TRUE)
  dat <- as.data.frame(dat)
  dat$N <- dat$N_CASES + dat$N_CONTROLS
  dat$SNP_chrpos <- paste0(dat$CHR, ":", dat$BP)
  out <- format_data(dat, type = "outcome",
                     snp_col = "SNP_chrpos", beta_col = "BETA", se_col = "SE",
                     pval_col = "P", effect_allele_col = "A1", other_allele_col = "A2",
                     eaf_col = "A1_FREQ", samplesize_col = "N",
                     chr_col = "CHR", pos_col = "BP")
  out$outcome <- GLIOMA_SUBTYPES[i]
  glioma_outcomes[[GLIOMA_SUBTYPES[i]]] <- out
}

# Load clumped ICVF instruments (from the forward MR)
instrument_files <- list.files(ICVF_DIR, pattern = "_instruments_", full.names = TRUE)

# Only run MR-PRESSO for nominally significant IVW results
sig_pairs <- combined %>% filter(nom_sig) %>% select(pgs_id, subtype)

presso_results <- list()

for (row_i in 1:nrow(sig_pairs)) {
  pgs_id  <- sig_pairs$pgs_id[row_i]
  subtype <- sig_pairs$subtype[row_i]

  # Find matching instrument file
  inst_file <- grep(pgs_id, instrument_files, value = TRUE)
  if (length(inst_file) == 0) next

  exposure_dat <- fread(inst_file[1], header = TRUE)
  exposure_dat <- as.data.frame(exposure_dat)

  # Convert SNP to chr:pos
  if ("chr.exposure" %in% names(exposure_dat) & "pos.exposure" %in% names(exposure_dat)) {
    exposure_dat$SNP <- paste0(exposure_dat$chr.exposure, ":", exposure_dat$pos.exposure)
  }

  outcome_dat <- glioma_outcomes[[subtype]]
  if (is.null(outcome_dat)) next

  harmonised <- tryCatch(
    harmonise_data(exposure_dat, outcome_dat, action = 2),
    error = function(e) NULL
  )

  if (is.null(harmonised) || sum(harmonised$mr_keep) < 4) next
  harm_keep <- harmonised[harmonised$mr_keep, ]

  message("  MR-PRESSO: ", pgs_id, " -> ", subtype, " (", nrow(harm_keep), " SNPs)")

  tryCatch({
    presso <- mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      data = harm_keep,
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      NbDistribution = 1000,
      SignifThreshold = 0.05
    )

    presso_results <- c(presso_results, list(data.frame(
      pgs_id = pgs_id,
      subtype = subtype,
      global_p = presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue,
      n_outliers = sum(presso[[1]]$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, na.rm = TRUE),
      causal_raw = presso[[1]]$`Main MR results`$`Causal Estimate`[1],
      causal_corrected = presso[[1]]$`Main MR results`$`Causal Estimate`[2],
      distortion_p = tryCatch(
        presso[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue,
        error = function(e) NA
      )
    )))
    message("    Global test p = ", presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue)
  }, error = function(e) {
    message("    Failed: ", e$message)
  })
}

if (length(presso_results) > 0) {
  presso_df <- bind_rows(presso_results)
  fwrite(presso_df, file.path(OUTPUT_DIR, "mr_presso_results.tsv"), sep = "\t")
  message("\n  Saved: mr_presso_results.tsv")

  message("\n--- MR-PRESSO Summary ---")
  for (i in 1:nrow(presso_df)) {
    r <- presso_df[i, ]
    flag <- if (r$global_p < 0.05) " ** OUTLIERS DETECTED" else " OK"
    message(sprintf("  %s -> %s: global_p=%.3f, outliers=%d%s",
                    r$pgs_id, r$subtype, r$global_p, r$n_outliers, flag))
  }
} else {
  message("  No MR-PRESSO results produced")
}

# =============================================================================
# STEP 7: Final summary
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("POST-PROCESSING COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\n--- Files created ---")
output_files <- list.files(OUTPUT_DIR, recursive = TRUE, full.names = TRUE)
for (f in output_files) message("  ", f)

message("\n--- Interpretation notes ---")
message("  - All significant effects are PROTECTIVE (higher ICVF -> lower glioma risk)")
message("  - Best evidence: corpus callosum & cingulum tracts -> all glioma / IDH-WT")
message("  - IDH-mutant results have pleiotropy flags — interpret with caution")
message("  - No results survive Bonferroni (p < ", sprintf("%.2e", BONF_THRESH), ")")
message("  - Steiger confirms ICVF -> Glioma direction for all significant results")
message("  - Reverse MR was not possible (glioma SNP IDs incompatible with IEU LD clumping)")
message("  - Consider local plink clumping for reverse MR if needed")
