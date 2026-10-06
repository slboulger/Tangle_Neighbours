#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_phf1_intensity_validation_neurons.R
#
# Figure panels: 1F
# Reads the per-cell p95 table written by plot_phf1_intensity_validation.R.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_phf1_intensity_validation_neurons.R
#
# Figure 1F: does the 95th-percentile PHF1 mask intensity recover the manual PHF1 +/- calls
# in NEURONS?
#
# Identical to the p95 arm of R/plot_phf1_intensity_validation.R (same violin + donor-mean
# figure, same rank AUC, same paired Wilcoxon on donor means) except that the cell set is the
# ten-label neuron pool (the nine named subtypes + "Unassigned Neuron"; NEURON_POOL_10 in
# R/nn3_utils.r, i.e. the cells that act as distance anchors). The robust z is NOT recomputed:
# phf1_intensity_p95_z keeps its Methods definition (within donor, over all cells), so the
# y-axis is on the same scale as the all-cell version and only the cell set changes.
#
# Input : plots/phf1_intensity_validation/source_data_phf1_intensity_p95_vs_manual.tsv
#         (the per-cell table written by R/plot_phf1_intensity_validation.R; avoids reading
#         the 7.8 GB object, so run that script first)
# Output: plots/phf1_intensity_validation_neurons/
#   plot_phf1_intensity_p95_vs_manual_neurons.pdf
#   source_data_phf1_intensity_p95_vs_manual_neurons.tsv               one row per drawn cell
#   source_data_phf1_intensity_p95_vs_manual_neurons_sample_means.tsv  per-donor means (the dots)
#   stats_phf1_intensity_p95_vs_manual_neurons.txt                     stats log (+ sessionInfo)
#   stats_phf1_intensity_p95_vs_manual_neurons_auc.tsv                 per-donor AUC / medians
#   stats_phf1_intensity_p95_vs_manual_neurons_effectsize.tsv          headline numbers + effect size

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
})
# HPC path first (canonical); the local mount is accepted so the script also runs off-cluster.
roots <- c("<PROJECT_ROOT>/phf1_v2",
           "<PROJECT_ROOT>/phf1_v2")
setwd(roots[dir.exists(roots)][1])
source("R/palettes.R")

NEURON_POOL_10 <- c(neuron_order, "Unassigned Neuron")   # == NEURON_POOL_10 in R/nn3_utils.r

indir  <- "plots/phf1_intensity_validation"
outdir <- "plots/phf1_intensity_validation_neurons"
name   <- "phf1_intensity_p95_vs_manual_neurons"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# rank-based ROC-AUC (no pROC dependency) - identical to R/plot_phf1_intensity_validation.R
auc_fn <- function(score, pos) {
  ok <- !is.na(score)
  score <- score[ok]; pos <- as.logical(pos[ok])
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

all_cells <- read.delim(file.path(indir, "source_data_phf1_intensity_p95_vs_manual.tsv"),
                        colClasses = c(sample_id = "character"))
stopifnot(nrow(all_cells) == 166840, sum(all_cells$PHF1) == 601)

df <- all_cells %>%
  filter(celltype %in% NEURON_POOL_10, !is.na(phf1_intensity_p95)) %>%
  transmute(cell_id, sample_id, Braak = factor(Braak, levels = braak_levels), celltype,
            PHF1 = as.logical(PHF1), int_raw = phf1_intensity_p95, int_z = phf1_intensity_p95_z)
stopifnot(sum(df$PHF1) == 398)   # the 398 tangle-bearing neurons that act as distance anchors

## --- stats: AUC overall and per donor -------------------------------------------------------
overall_auc <- auc_fn(df$int_raw, df$PHF1)
per_sample <- df %>%
  group_by(sample_id) %>%
  summarise(n = n(), n_phf1_pos = sum(PHF1),
            median_pos = median(int_raw[PHF1]), median_neg = median(int_raw[!PHF1]),
            auc = auc_fn(int_raw, PHF1), .groups = "drop") %>%
  arrange(desc(auc))
write.table(per_sample, file.path(outdir, sprintf("stats_%s_auc.tsv", name)),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- source data (one row per drawn cell) ---------------------------------------------------
write.table(df %>% rename(phf1_intensity_p95 = int_raw, phf1_intensity_p95_z = int_z),
            file.path(outdir, sprintf("source_data_%s.tsv", name)),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- per-donor means and the paired test (the inferential unit is the donor) ----------------
pd <- df %>% mutate(grp = factor(ifelse(PHF1, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")))
sm <- df %>%
  group_by(sample_id, PHF1) %>%
  summarise(sample_mean = mean(int_z, na.rm = TRUE), .groups = "drop") %>%
  mutate(grp = factor(ifelse(PHF1, "PHF1+", "PHF1-"), levels = c("PHF1-", "PHF1+")))
write.table(sm, file.path(outdir, sprintf("source_data_%s_sample_means.tsv", name)),
            sep = "\t", quote = FALSE, row.names = FALSE)
sm_wide <- sm %>% select(sample_id, grp, sample_mean) %>%
  pivot_wider(names_from = grp, values_from = sample_mean)
stopifnot(nrow(sm_wide) == 9, !anyNA(sm_wide))
paired_w <- wilcox.test(sm_wide[["PHF1+"]], sm_wide[["PHF1-"]], paired = TRUE)
p_lab <- paste0("p = ", formatC(paired_w$p.value, format = "g", digits = 2))

# effect size: paired difference of donor means with t-based 95% CI, donor-paired dz
dd <- sm_wide[["PHF1+"]] - sm_wide[["PHF1-"]]
tt <- t.test(dd)
pos_z <- df$int_z[df$PHF1]
eff <- tibble(
  cell_set = "neurons (NEURON_POOL_10)",
  n_cells = nrow(df), n_tangle_bearing = sum(df$PHF1), n_donors = nrow(sm_wide),
  median_z_tangle_bearing = median(pos_z), q25_z_tangle_bearing = unname(quantile(pos_z, 0.25)),
  q75_z_tangle_bearing = unname(quantile(pos_z, 0.75)),
  auc_pooled = overall_auc, auc_donor_min = min(per_sample$auc), auc_donor_max = max(per_sample$auc),
  median_p95_pos_min = min(per_sample$median_pos), median_p95_pos_max = max(per_sample$median_pos),
  median_p95_neg_min = min(per_sample$median_neg), median_p95_neg_max = max(per_sample$median_neg),
  paired_diff_mean = mean(dd), paired_diff_ci_lo = tt$conf.int[1], paired_diff_ci_hi = tt$conf.int[2],
  dz = mean(dd) / sd(dd), n_donors_same_direction = sum(sign(dd) == sign(mean(dd))),
  wilcox_paired_V = unname(paired_w$statistic), wilcox_paired_p = paired_w$p.value)
write.table(eff, file.path(outdir, sprintf("stats_%s_effectsize.tsv", name)),
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- plot: identical encoding to the published panel ----------------------------------------
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
  labs(x = NULL, y = "PHF1 intensity, 95th pct (per-sample z)") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"))
ggsave(file.path(outdir, sprintf("plot_%s.pdf", name)),
       p, width = 4, height = 6.14, units = "cm", device = "pdf")

## --- stats log ----------------------------------------------------------------------------
sink(file.path(outdir, sprintf("stats_%s.txt", name)))
cat("PHF1 p95 mask-intensity validation against manual PHF1 +/- calls - NEURONS ONLY\n")
cat("================================================================================\n\n")
cat("Cell set: NEURON_POOL_10 (", paste(NEURON_POOL_10, collapse = ", "), ")\n", sep = "")
cat("Robust z: phf1_intensity_p95_z as published (within donor, over ALL cells) - not recomputed.\n\n")
cat(sprintf("Neurons with extracted intensity: %d | manual PHF1+ among them: %d | donors: %d\n\n",
            nrow(df), sum(df$PHF1), nrow(sm_wide)))
cat(sprintf("Median robust z of tangle-bearing neurons: %.2f (IQR %.2f-%.2f)\n",
            eff$median_z_tangle_bearing, eff$q25_z_tangle_bearing, eff$q75_z_tangle_bearing))
cat(sprintf("Pooled ROC-AUC (raw p95 vs manual call): %.4f | per-donor range %.3f-%.3f\n",
            overall_auc, eff$auc_donor_min, eff$auc_donor_max))
cat(sprintf("Per-donor median raw p95: tangle-bearing %.1f-%.1f, tangle-free %.1f-%.1f\n\n",
            eff$median_p95_pos_min, eff$median_p95_pos_max, eff$median_p95_neg_min, eff$median_p95_neg_max))
cat("Per-donor AUC / medians:\n")
print(as.data.frame(per_sample), row.names = FALSE)
cat(sprintf("\nFigure statistic - two-sided paired Wilcoxon signed-rank on per-donor mean z (PHF1+ vs PHF1-),\n  n = %d donors, V = %g, p = %s\n",
            nrow(sm_wide), eff$wilcox_paired_V, signif(paired_w$p.value, 4)))
cat(sprintf("Effect size: paired difference %.2f z (95%% CI %.2f to %.2f), dz = %.2f, %d/%d donors same direction\n",
            eff$paired_diff_mean, eff$paired_diff_ci_lo, eff$paired_diff_ci_hi, eff$dz,
            eff$n_donors_same_direction, eff$n_donors))
cat("(Dots are per-donor means; violins show all neurons; the paired donor-level test is the\n")
cat(" inferential one.)\n\n")
cat("------------------------------------------------------------\n")
print(sessionInfo())
sink()

cat(sprintf("[%s] neurons %d (PHF1+ %d) | median z %.2f [%.2f, %.2f] | AUC %.4f (donors %.3f-%.3f) | paired Wilcoxon p = %.4g\n",
            name, nrow(df), sum(df$PHF1), eff$median_z_tangle_bearing, eff$q25_z_tangle_bearing,
            eff$q75_z_tangle_bearing, overall_auc, eff$auc_donor_min, eff$auc_donor_max, paired_w$p.value))
