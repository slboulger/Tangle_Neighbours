#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# replot_intensity_adjust_forest_otero.R
#
# Figure panels: S1A
# Redraws from the source data of phf1_module_distance_intensity_adjust.R; no refit.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# replot_intensity_adjust_forest_otero.R
#
# The Otero-Garcia Exc1/2 UP distance slope before and after per-cell PHF1-intensity adjustment,
# tangle-free Exc-IT-L2-3-CBLN2-HOPX neurons. Redrawn from the source data of
# R/phf1_module_distance_intensity_adjust.R -- no refit, no Seurat object.
# make_forest_plot(), save_fixed_panel() and the colours are copied verbatim from that script;
# only the two Otero-Garcia rows are drawn.
#
#   Rscript phf1_v2/R/replot_intensity_adjust_forest_otero.R

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

CT      <- "Exc-IT-L2-3-CBLN2-HOPX"
IN_DIR  <- "plots/phf1_distance_intensity_adjust_1000um"
LAB     <- "Otero-Garcia Exc1/2 UP"
COL     <- "#0072B2"                     # SIG_DEF$otero_up$colour
HEIGHT_IN <- 1.3                         # canvas height (2 rows)
WIDTH_CM  <- 6.4                         # total canvas width
SHORT     <- "Exc1/2 UP"                 # y labels drop "Otero-Garcia"

# --- verbatim from R/phf1_module_distance_intensity_adjust.R ---------------------------------
darken <- function(col, f = 0.45) {
  v <- as.numeric(grDevices::col2rgb(col)) * f
  grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
}
make_forest_plot <- function(fdf, key_cols, xlab) {
  lv <- names(key_cols)
  fdf$series <- factor(fdf$series, levels = rev(lv))       # first entry at the TOP
  ggplot(fdf, aes(beta, series, colour = series)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_segment(aes(x = ci_lo, xend = ci_hi, y = series, yend = series), linewidth = 0.4) +
    geom_point(shape = 16, size = 2.2) +
    scale_colour_manual(values = key_cols, guide = "none", drop = FALSE) +
    labs(x = xlab, y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_text(size = 6),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
}
save_fixed_panel <- function(p, path, height = 2.1) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  ggsave(path, g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = height, units = "in", device = "pdf")
}
# ---------------------------------------------------------------------------------------------

src <- read.delim(file.path(IN_DIR, sprintf("source_data_distance_intensity_adjust_forest_%s.tsv", CT)),
                  check.names = FALSE, stringsAsFactors = FALSE)
fdf <- src %>% filter(startsWith(series, LAB))
stopifnot(nrow(fdf) == 2, setequal(fdf$series, paste(LAB, c("(unadjusted)", "(adjusted)"))))
# display labels without "Otero-Garcia"; the original series name is kept in the source data
fdf$series_original <- fdf$series
fdf$series <- sub(paste0("^", LAB), SHORT, fdf$series)
key_cols <- setNames(c(COL, darken(COL)), paste(SHORT, c("(unadjusted)", "(adjusted)")))

stem <- sprintf("distance_intensity_adjust_forest_%s_otero_only", CT)
# three x ticks (ASCII minus) and a right-aligned x title so they fit the narrow canvas
p_s4a <- make_forest_plot(fdf, key_cols, expression(paste(beta, ", per S.D. log-distance"))) +
  scale_x_continuous(breaks = c(-0.004, -0.002, 0), labels = c("-0.004", "-0.002", "0")) +
  theme(axis.title.x = element_text(hjust = 1))
ggsave(file.path(IN_DIR, paste0("plot_", stem, ".pdf")), p_s4a,
       width = WIDTH_CM, height = HEIGHT_IN * 2.54, units = "cm", device = "pdf")
write.table(fdf, file.path(IN_DIR, paste0("source_data_", stem, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summ <- read.delim(file.path(IN_DIR, "stats_distance_intensity_adjust_summary.tsv"), stringsAsFactors = FALSE)
s <- summ[summ$celltype == CT & summ$sig_label == LAB, ]
stopifnot(nrow(s) == 1,
          abs(s$dist_slope_unadj - fdf$beta[grepl("unadjusted", fdf$series_original)]) < 1e-12,
          abs(s$dist_slope_adj   - fdf$beta[grepl("\\(adjusted", fdf$series_original)]) < 1e-12)

sink(file.path(IN_DIR, paste0("stats_", stem, ".txt")))
cat("Fig. S1A, Otero-Garcia Exc1/2 UP only -- distance slope before vs after PHF1-intensity adjustment\n")
cat("Redrawn from source_data_distance_intensity_adjust_forest_", CT, ".tsv (R/phf1_module_distance_intensity_adjust.R);\n",
    "no refit. Full model, covariates and cell set: stats_distance_intensity_adjust.txt.\n\n", sep = "")
cat(sprintf("n = %d tangle-free %s neurons (%d dropped for missing intensity), 9 donors\n",
            s$n_cells, CT, s$n_intensity_dropped))
cat(sprintf("Unadjusted: beta = %.3g per SD log-distance [%.3g, %.3g], p = %.2g\n",
            s$dist_slope_unadj, s$dist_slope_unadj_ci_lo, s$dist_slope_unadj_ci_hi, s$dist_p_unadj))
cat(sprintf("Adjusted:   beta = %.3g per SD log-distance [%.3g, %.3g], p = %.2g\n",
            s$dist_slope_adj, s$dist_slope_adj_ci_lo, s$dist_slope_adj_ci_hi, s$dist_p_adj))
cat("Points, estimate; bars, 95% Wald CI (two-sided). Negative = higher score nearer tangles.\n\n")
print(sessionInfo())
sink()
cat("Done\n")
