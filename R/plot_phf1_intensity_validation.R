#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_phf1_intensity_validation.R
#
# Figure panels: 1F - upstream step
# Writes the per-cell p95 table (plots/phf1_intensity_validation/) that
# plot_phf1_intensity_validation_neurons.R draws as Fig. 1F. Its own all-cell panels are
# not in the paper.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_phf1_intensity_validation.R
#
# Validation / acceptance figure for the extracted PHF1 mask intensity:
# does the continuous intensity recover the manual PHF1 +/- calls?
#
# Produces the standard plot TRIPLE under plots/phf1_intensity_validation/, once
# per intensity statistic (mean and 95th-percentile), via run_validation():
#   plot_<name>.pdf                    pooled (all samples) violin of per-sample-z intensity by
#                                      manual PHF1 status; dots = per-sample means (mean +/- SEM);
#                                      paired-Wilcoxon p (on the 9 donor means) shown above the bracket
#   source_data_<name>.tsv             long-format, one row per drawn cell (violin)
#   source_data_<name>_sample_means.tsv per-sample means (the dots / test input)
#   stats_<name>.txt                   AUC + per-sample Wilcoxon + paired test + sessionInfo()
#   stats_<name>_auc.tsv               structured AUC/medians table (supp-table ready)
#
# names: phf1_intensity_vs_manual      (mean intensity: phf1_intensity_mean / _z)
#        phf1_intensity_p95_vs_manual  (95th-pct intensity: phf1_intensity_p95 / _z)

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(tidyr); library(ggplot2)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

outdir <- "plots/phf1_intensity_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# rank-based ROC-AUC (no pROC dependency)
auc_fn <- function(score, pos) {
  ok <- !is.na(score)
  score <- score[ok]; pos <- as.logical(pos[ok])
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

seu <- readRDS("seu_PHF1.rds")
md <- seu@meta.data
if (!"cell_id" %in% colnames(md)) md$cell_id <- rownames(md)   # seu already has a cell_id col

# run the full validation triple for one intensity statistic.
#   raw_col : raw intensity column, used for AUC / Wilcoxon / NA filter (e.g. phf1_intensity_mean)
#   z_col   : per-sample robust-z column, plotted on the y-axis (e.g. phf1_intensity_mean_z)
#   name    : output filename stem
#   ylab    : y-axis label
# Internally the two columns are carried as int_raw / int_z and renamed back to their
# real names only when writing the cell-level source data, so outputs are unchanged.
run_validation <- function(raw_col, z_col, name, ylab) {
  df <- md %>%
    transmute(
      cell_id   = cell_id,
      sample_id = as.character(sample_id),
      Braak     = factor(Braak, levels = braak_levels),
      celltype  = celltype,
      PHF1      = as.logical(PHF1 %in% c(TRUE, "TRUE", "True", 1, "1")),
      int_raw   = .data[[raw_col]],
      int_z     = .data[[z_col]]
    ) %>%
    filter(!is.na(int_raw))

  ## --- stats: AUC + Wilcoxon overall and per sample ---
  overall_auc <- auc_fn(df$int_raw, df$PHF1)

  per_sample <- df %>%
    group_by(sample_id) %>%
    summarise(
      n            = n(),
      n_phf1_pos   = sum(PHF1),
      median_pos   = median(int_raw[PHF1], na.rm = TRUE),
      median_neg   = median(int_raw[!PHF1], na.rm = TRUE),
      auc          = auc_fn(int_raw, PHF1),
      wilcox_p     = tryCatch(
        wilcox.test(int_raw ~ PHF1)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    arrange(desc(auc))

  write.table(per_sample,
              file.path(outdir, sprintf("stats_%s_auc.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ## --- source data (long; one row per drawn cell) ---
  write.table(df %>% rename(!!raw_col := int_raw, !!z_col := int_z),
              file.path(outdir, sprintf("source_data_%s.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ## --- plot: pooled across samples; violin = all cells, dots = per-sample means ---
  pd <- df %>% mutate(grp = factor(ifelse(PHF1, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")))

  # per-sample means (the inferential unit: each donor gives a paired +/- pair)
  sm <- df %>%
    group_by(sample_id, PHF1) %>%
    summarise(sample_mean = mean(int_z, na.rm = TRUE), .groups = "drop") %>%
    mutate(grp = factor(ifelse(PHF1, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")))
  write.table(sm, file.path(outdir, sprintf("source_data_%s_sample_means.tsv", name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # paired Wilcoxon on the per-sample means (donor-level). Full result -> stats log;
  # the figure shows the p-value only (no test name / n, to stay clean).
  sm_wide <- sm %>% select(sample_id, grp, sample_mean) %>%
    tidyr::pivot_wider(names_from = grp, values_from = sample_mean)
  paired_w <- tryCatch(wilcox.test(sm_wide[["PHF1+"]], sm_wide[["PHF1-"]], paired = TRUE),
                       error = function(e) NULL)
  p_lab <- if (is.null(paired_w)) "" else paste0("p = ", formatC(paired_w$p.value, format = "g", digits = 2))

  # group mean +/- SEM across the per-sample means (donor-level uncertainty)
  grp_mean <- sm %>% group_by(grp) %>%
    summarise(m = mean(sample_mean), sem = sd(sample_mean) / sqrt(n()), .groups = "drop") %>%
    mutate(xc = as.integer(grp))
  y_lo <- min(df$int_z, na.rm = TRUE)
  y_hi <- max(quantile(df$int_z, 0.99, na.rm = TRUE), max(sm$sample_mean, na.rm = TRUE))
  y_br <- y_hi * 1.06; tick <- (y_hi - y_lo) * 0.02

  p <- ggplot(pd, aes(x = grp, y = int_z)) +
    geom_violin(aes(fill = grp), scale = "width", linewidth = 0.2, alpha = 0.6, colour = NA) +
    geom_errorbar(data = grp_mean, aes(x = grp, ymin = m - sem, ymax = m + sem),
                  inherit.aes = FALSE, width = 0.16, linewidth = 0.3) +
    geom_segment(data = grp_mean, aes(x = xc - 0.2, xend = xc + 0.2, y = m, yend = m),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_jitter(data = sm, aes(x = grp, y = sample_mean), width = 0.1, height = 0,
                shape = 21, size = 1.5, fill = "white", colour = "black", stroke = 0.3) +
    annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br, linewidth = 0.3) +
    annotate("segment", x = 1, xend = 1, y = y_br, yend = y_br - tick, linewidth = 0.3) +
    annotate("segment", x = 2, xend = 2, y = y_br, yend = y_br - tick, linewidth = 0.3) +
    annotate("text", x = 1.5, y = y_br, label = p_lab, vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = c("PHF1-" = "grey75", "PHF1+" = "#BD0026"), guide = "none") +
    coord_cartesian(ylim = c(y_lo, y_br * 1.12)) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 8) + fig_theme +          # match other figures' text sizes
    # short 2-group x labels -> keep horizontal; drop fig_theme's big left margin
    # (that is only needed for the 45deg celltype labels elsewhere)
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"))

  ggsave(file.path(outdir, sprintf("plot_%s.pdf", name)),
         p, width = 4, height = 6.14, units = "cm", device = "pdf")

  ## --- human-readable stats log ---
  sink(file.path(outdir, sprintf("stats_%s.txt", name)))
  cat("PHF1 mask-intensity validation against manual PHF1 +/- calls\n")
  cat("===========================================================\n\n")
  cat(sprintf("Cells with extracted intensity: %d\n", nrow(df)))
  cat(sprintf("Manual PHF1+ among them: %d\n\n", sum(df$PHF1)))
  cat(sprintf("Overall ROC-AUC (%s vs manual PHF1): %.4f\n", raw_col, overall_auc))
  cat("(>0.5 = higher intensity in manual-positives; ~0.5 = no separation / misregistration)\n\n")
  cat("Per-sample AUC / medians / Wilcoxon:\n")
  print(per_sample, row.names = FALSE)
  cat(sprintf("\nPooled figure statistic — paired Wilcoxon on per-sample mean z (PHF1+ vs PHF1-),\n  n = %d donors, p = %s\n",
              nrow(sm_wide), if (is.null(paired_w)) "NA" else signif(paired_w$p.value, 4)))
  cat("(The figure's dots are these per-sample means; the violins show the full cell-level\n")
  cat(" distribution. The paired donor-level test is the inferential one.)\n")
  cat("\nAcceptance: every sample should show median_pos >> median_neg and AUC well above 0.5.\n")
  cat("Samples failing this are flagged here and should be excluded / re-registered,\n")
  cat("not silently merged.\n\n")
  cat("------------------------------------------------------------\n")
  print(sessionInfo())
  sink()

  cat(sprintf("[%s] Overall AUC = %.4f. Wrote triple to %s/\n", name, overall_auc, outdir))
  invisible(overall_auc)
}

# mean intensity + 95th-percentile intensity (bright-tail statistic)
run_validation("phf1_intensity_mean", "phf1_intensity_mean_z",
               "phf1_intensity_vs_manual", "PHF1 intensity (per-sample z)")
run_validation("phf1_intensity_p95", "phf1_intensity_p95_z",
               "phf1_intensity_p95_vs_manual", "PHF1 intensity, 95th pct (per-sample z)")
