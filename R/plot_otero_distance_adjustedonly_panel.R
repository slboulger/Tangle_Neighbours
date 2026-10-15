#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_otero_distance_adjustedonly_panel.R
#
# Figure panels: 3B
# Reads the outputs of phf1_module_distance_intensity_adjust.R; no model is fitted here.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_otero_distance_adjustedonly_panel.R
#
# Otero-Garcia Exc1/2 UP module score over distance to the nearest PHF1+ neuron in tangle-free
# Exc-IT-L2-3-CBLN2-HOPX neurons: the intensity-adjusted series of
# plot_distance_intensity_adjust_combined_<ct>_adjustedonly.pdf, drawn on its own (Otero-Garcia
# signature only, no PHF1+ reference lines).
#
# Pure re-plot: reads the rolling-mean source data written by
# R/phf1_module_distance_intensity_adjust.R (intensity-ADJUSTED rolling mean, 75 um window,
# 95% CI = +/- 1.96 SEM). Same colour (darken(#0072B2, 0.45)), theme, legend and
# save_fixed_panel() geometry as that script. The stats log quotes its LMM row.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # fig_theme

CT         <- "Exc-IT-L2-3-CBLN2-HOPX"
SIG        <- "otero_up"
LABEL      <- "Otero-Garcia Exc1/2 UP"
SRC_KEY    <- paste0(LABEL, " (adjusted)")                 # the adjusted series is drawn
SIG_COLOUR <- "#0072B2"                                    # SIG_DEF$otero_up$colour in phf1_module_distance_intensity_adjust.R
X_CAP      <- 1000
PARENT_DIR <- "plots/phf1_distance_intensity_adjust_1000um"
OUTPUT_DIR <- "plots/otero_distance_adjustedonly_panel"
NAME       <- sprintf("otero_distance_adjustedonly_%s", CT)
XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# identical to the helpers in phf1_module_distance_intensity_adjust.R
darken <- function(col, f = 0.45) {
  v <- as.numeric(grDevices::col2rgb(col)) * f
  grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
}
save_fixed_panel <- function(p, path) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  ggsave(path, g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")
}

roll_file <- file.path(PARENT_DIR, sprintf("source_data_distance_intensity_adjust_combined_%s_rollmean.tsv", CT))
roll <- read.delim(roll_file, check.names = FALSE) %>%
  filter(key == SRC_KEY) %>%
  mutate(key = factor(LABEL, levels = LABEL))
if (nrow(roll) == 0) stop("No rows for key '", SRC_KEY, "' in ", roll_file)

key_cols <- setNames(darken(SIG_COLOUR, 0.45), LABEL)

# make_adjusted_only_plot() from phf1_module_distance_intensity_adjust.R, minus the two reference-band layers
p <- ggplot(roll) +
  geom_ribbon(aes(dist_to_phf1_um, ymin = roll_mean - 1.96 * sem,
                  ymax = roll_mean + 1.96 * sem, fill = key), alpha = 0.15, colour = NA) +
  geom_line(aes(dist_to_phf1_um, roll_mean, colour = key), linewidth = 0.7) +
  scale_colour_manual(values = key_cols, name = NULL, drop = FALSE, limits = LABEL) +
  scale_fill_manual(values = key_cols, guide = "none", drop = FALSE) +
  labs(x = XLAB, y = "Module score") +
  coord_cartesian(xlim = c(0, X_CAP)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.text = element_text(size = 6), legend.key.width = grid::unit(14, "pt"),
        legend.key.height = grid::unit(9, "pt"), legend.position = "right",
        legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 5))

save_fixed_panel(p, file.path(OUTPUT_DIR, sprintf("plot_%s.pdf", NAME)))

write.table(roll %>% mutate(key = as.character(key), source_key = SRC_KEY,
                            ci_lo = roll_mean - 1.96 * sem, ci_hi = roll_mean + 1.96 * sem) %>%
              select(celltype, key, source_key, window_um, dist_to_phf1_um, roll_mean, sem,
                     ci_lo, ci_hi, n_window),
            file.path(OUTPUT_DIR, sprintf("source_data_%s.tsv", NAME)),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- stats log: the LMM row for this curve from phf1_module_distance_intensity_adjust.R (no refit) ----
summ <- read.delim(file.path(PARENT_DIR, "stats_distance_intensity_adjust_summary.tsv"), check.names = FALSE) %>%
  filter(celltype == CT, signature == SIG)
if (nrow(summ) != 1) stop("Expected one summary row for ", CT, " / ", SIG)

sink(file.path(OUTPUT_DIR, sprintf("stats_%s.txt", NAME)))
cat("Otero-Garcia Exc1/2 UP module score over distance to PHF1+ neuron -", CT, "\n")
cat("Re-plot of the adjusted series in ", PARENT_DIR, "/plot_distance_intensity_adjust_combined_", CT,
    "_adjustedonly.pdf\n", sep = "")
cat("with the PHF1 signature and the dotted PHF1+ reference lines/bands removed. No refit.\n\n")
cat("Drawn: model-free rolling mean of the PHF1-intensity-ADJUSTED module score\n")
cat("  (score - b_intensity * (intensity - mean)), window", unique(roll$window_um), "um,",
    "ribbon = mean +/- 1.96 SEM, grid of", nrow(roll), "points over 0 -", X_CAP, "um.\n")
cat("  Cells per window:", min(roll$n_window), "-", max(roll$n_window), "\n")
cat("  Source:", roll_file, "\n\n")
cat("Models (R/phf1_module_distance_intensity_adjust.R, Set 1; distance = log(d)/sd(log(d)), dist_sd =",
    sprintf("%.4f", summ$dist_sd), "):\n")
cat("  unadjusted: score ~ dist_scaled + nUMI_log + percent_neg + Sex + Age_s + PMI_s + (1 | sample_id)\n")
cat("  adjusted  : score ~ dist_scaled + phf1_intensity_p95_z + nUMI_log + percent_neg + Sex + Age_s + PMI_s + (1 | sample_id)\n")
cat("  n cells =", summ$n_cells, "\n\n")
cat("=== Effect sizes (distance slope per SD log-distance, 95% CI) ===\n")
cat(sprintf("  unadjusted: beta = %.3g [%.3g, %.3g], p = %.2g\n",
            summ$dist_slope_unadj, summ$dist_slope_unadj_ci_lo, summ$dist_slope_unadj_ci_hi, summ$dist_p_unadj))
cat(sprintf("  adjusted  : beta = %.3g [%.3g, %.3g], p = %.2g   <- the drawn series\n",
            summ$dist_slope_adj, summ$dist_slope_adj_ci_lo, summ$dist_slope_adj_ci_hi, summ$dist_p_adj))
cat(sprintf("  per log-unit (adjusted): %.3g [%.3g, %.3g]\n",
            summ$dist_slope_adj_perlog, summ$dist_ci_lo_adj_perlog, summ$dist_ci_hi_adj_perlog))
cat(sprintf("  attenuation by intensity adjustment: %.1f%% (%s)\n", summ$pct_attenuation, summ$verdict))
cat("  Cell-level fit on ~2.5e4 cells from 9 donors.\n\n")
sink(); capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, sprintf("stats_%s.txt", NAME)), append = TRUE)
cat("Wrote", OUTPUT_DIR, "\n")
