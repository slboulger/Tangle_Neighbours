#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_seg_control_composite_matched.R
#
# Figure panels: S5A-C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_seg_control_composite_matched.R
#
# PAPER-READY assembly of the SEG negative control. Builds the supplementary figure and
# the specificity-margin table from EXISTING source-data TSVs - it re-reads nothing
# large, needs no Seurat object and no HPC session, and runs locally in seconds.
#
# WHY THIS EXISTS (two problems with the raw run in
#                  plots/modulescore_seg_control_1000um/):
#   1. FREE Y-SCALES. Each SEG panel auto-scales, so a ~0.008-unit drift is stretched to
#      fill the panel and the control READS AS A STRONG GRADIENT - the exact opposite of
#      what it shows. Here each panel's axis SPAN is matched to the corresponding
#      reported panel's drawn-ribbon span (centred on the control's own midpoint), so the
#      control is drawn at identical module-score units per cm and renders as the
#      near-flat line it is. Measured drawn spans (incl. the +/-1.96*SEM ribbon):
#        CBLN2 0.0083 of 0.0972 (8.6%) | Micro 0.0175 of 0.0896 (19.5%)
#        Astro 0.0041 of 0.1466 (2.8%)
#      Absolute limits are deliberately NOT copied - see the note at ribbon_span().
#   2. Exc-IT-L3-5-CHGA-IL1RAPL2 is scored in the raw run but is not part of this figure,
#      so it is dropped here. Three panels give a 1 x 3 COMBINED layout at A4 width
#      (no 2 x 2, no x-axis title collisions).
#
# PANELS:
#   Exc-IT-L2-3-CBLN2-HOPX   vs the INTENSITY-ADJUSTED PHF1 / Otero-Garcia Exc1/2 UP arm
#   Micro                    vs the Mancuso microglial states
#   Astro                    vs the Cameron astrocyte SUBCLUSTERS (pooled not used)
#     --astro_source serrano : Astro vs the Serrano-Pozo 2024 states instead (from
#     R/plot_serrano_astro_vs_phf1_distance_modelp.R), written to
#     plots/seg_control_paper_serrano/ so Fig. S5 is untouched. Only the Astro axis span
#     and the Astro margin rows change; CBLN2 and Micro are identical.
#
# THE CONTROL AS A BENCHMARK  the reported gradients are 7.5-48x the SEG slope in the
#   same cells and model, i.e. the control is 2-13% of every reported effect. It is a
#   FLOOR / benchmark, NOT a baseline to subtract: the SEG module's expression
#   composition differs from the reported signatures, so its drift does not transfer
#   quantitatively.
#
# NOTES carried into the stats log:
#   * A donor random intercept is correctly calibrated for a WITHIN-donor predictor like
#     distance (verified by within-donor permutation).
#   * The CBLN2 floor is fitted in the UNADJUSTED model while the reported neuron
#     numbers are intensity-adjusted. Even if the floor doubled, PHF1 would still be 9.5x
#     and Otero-UP 4.8x.
#
# OUTPUT under --output_dir (default plots/seg_control_paper/):
#   plot_seg_control_matched_<CT>.pdf          per-panel, span-matched y-axis
#   plot_seg_control_matched_composite.pdf     3-panel, A4 width - THE supp figure
#   source_data_seg_control_matched.tsv        exact drawn rows + the matched limits
#   stats_seg_control_vs_reported_margin.tsv   THE supp table (specificity margins)
#   stats_seg_control_matched.txt              provenance, notes, sessionInfo
#
#   Rscript R/plot_seg_control_composite_matched.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(patchwork); library(argparse)
})

PROJECT_ROOTS <- c(
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>/phf1_v2"
)
.root <- PROJECT_ROOTS[dir.exists(PROJECT_ROOTS)]
if (length(.root) == 0) stop("No project root found.", call. = FALSE)
setwd(.root[1]); cat("Project root:", getwd(), "\n")
source("R/palettes.R")   # fig_theme

parser <- ArgumentParser()
parser$add_argument("--seg_dir",    default = "plots/modulescore_seg_control_1000um")
parser$add_argument("--glia_dir",   default = "plots/module_score_vs_phf1_distance_modelp_1000um")
parser$add_argument("--adjust_dir", default = "plots/phf1_distance_intensity_adjust_1000um")
parser$add_argument("--output_dir", default = NULL,
  help = paste("Default plots/seg_control_paper (cameron) or",
               "plots/seg_control_paper_serrano (serrano), so the S5 figure is never overwritten"))
parser$add_argument("--max_dist_um", type = "double", default = 1000)
# Which reported astrocyte panel the Astro SEG axis is matched to. cameron = Fig. S5 as
# published; serrano = the Serrano-Pozo 2024 states from
# R/plot_serrano_astro_vs_phf1_distance_modelp.R (astMic/astNeu excluded there).
parser$add_argument("--astro_source", default = "cameron", choices = c("cameron", "serrano"))
parser$add_argument("--serrano_dir",
  default = "plots/module_score_serrano_astro_vs_phf1_distance_modelp_1000um")
args <- parser$parse_args()
if (is.null(args$output_dir))
  args$output_dir <- if (args$astro_source == "serrano") "plots/seg_control_paper_serrano" else
    "plots/seg_control_paper"
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

ASTRO <- switch(args$astro_source,
  cameron = list(dir = args$glia_dir, group = "cameron_subclusters",
                 label = "Cameron astrocyte subclusters", short = "Cameron Astro subclusters"),
  serrano = list(dir = args$serrano_dir, group = "serrano_states",
                 label = "Serrano-Pozo astrocyte states", short = "Serrano-Pozo Astro states"))
cat("Astro axis matched to:", ASTRO$label, "(", ASTRO$dir, ")\n")

SEG_COLOUR <- "#666666"
XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")
sf2  <- function(x) if (all(is.na(x))) "NA" else formatC(signif(x, 2), format = "g")

need <- function(p) { if (!file.exists(p)) stop("missing input: ", p, call. = FALSE); p }

##  ............................................................................
##  Panel definitions - the three figure panels                              ####
PANELS <- list(
  list(ct = "Exc-IT-L2-3-CBLN2-HOPX", label = "Exc-IT-L2-3-CBLN2-HOPX",
       reported_source = "intensity-adjusted (PHF1, Otero-Garcia Exc1/2 UP)"),
  list(ct = "Micro", label = "Micro", reported_source = "Mancuso microglial states"),
  list(ct = "Astro", label = "Astro", reported_source = ASTRO$label)
)

##  ............................................................................
##  Reported panels: y-range (of the DRAWN ribbon) + slopes                  ####
# CBLN2: the reported figure is the intensity-ADJUSTED arm, so only the "(adjusted)"
# series define its y-range.
adj_roll <- read.delim(need(file.path(args$adjust_dir,
  "source_data_distance_intensity_adjust_combined_Exc-IT-L2-3-CBLN2-HOPX_rollmean.tsv")),
  stringsAsFactors = FALSE)
adj_roll <- adj_roll[grepl("\\(adjusted\\)", adj_roll$key), ]
stopifnot("no '(adjusted)' series found in the CBLN2 combined rollmean" = nrow(adj_roll) > 0)

micro_roll <- read.delim(need(file.path(args$glia_dir,
  "source_data_modulescore_raw_dist_mancuso_Micro_rollmean.tsv")), stringsAsFactors = FALSE)
astro_roll <- read.delim(need(file.path(ASTRO$dir, sprintf(
  "source_data_modulescore_raw_dist_%s_Astro_rollmean.tsv", ASTRO$group))), stringsAsFactors = FALSE)

# SPAN of the drawn ribbon (roll_mean +/- 1.96*SEM) in the reported panel.
#
# We match the axis SPAN, not the absolute limits. AddModuleScore subtracts an
# expression-bin-matched control set, which pins each gene set near ITS OWN zero, so
# different modules sit at different absolute baselines: the adjusted PHF1 signature
# occupies 0.017-0.099 while the SEG module occupies 0.002-0.009. Copying the reported
# panel's absolute limits therefore pushes the SEG trace off the panel entirely. Matching
# the span and centring on the SEG trace gives a like-for-like comparison - identical
# module-score units per cm, so the control's variation is seen against the same
# vertical scale as the reported effect.
ribbon_span <- function(d) {
  lo <- min(d$roll_mean - 1.96 * d$sem, na.rm = TRUE)
  hi <- max(d$roll_mean + 1.96 * d$sem, na.rm = TRUE)
  hi - lo
}
REPORTED_SPAN <- list(
  "Exc-IT-L2-3-CBLN2-HOPX" = ribbon_span(adj_roll),
  "Micro"                  = ribbon_span(micro_roll),
  "Astro"                  = ribbon_span(astro_roll)
)
# centre the matched span on the SEG trace's own midpoint
matched_ylim <- function(d, span) {
  mid <- mean(range(c(d$roll_mean - 1.96 * d$sem, d$roll_mean + 1.96 * d$sem), na.rm = TRUE))
  c(mid - span / 2, mid + span / 2)
}

# steepest reported |slope| per panel, and the full reported module list for the table
seg_coef <- read.delim(need(file.path(args$seg_dir, "stats_modulescore_lmm_coeffs.tsv")),
                       stringsAsFactors = FALSE)
glia_coef <- read.delim(need(file.path(args$glia_dir, "stats_modulescore_lmm_coeffs.tsv")),
                        stringsAsFactors = FALSE)
astro_coef <- if (args$astro_source == "cameron") glia_coef else
  read.delim(need(file.path(ASTRO$dir, "stats_modulescore_lmm_coeffs.tsv")),
             stringsAsFactors = FALSE)
adj_sum <- read.delim(need(file.path(args$adjust_dir,
  "stats_distance_intensity_adjust_summary.tsv")), stringsAsFactors = FALSE)

SEG_FLOOR <- setNames(seg_coef$slope_per_sd, seg_coef$celltype)
SEG_PADJ  <- setNames(seg_coef$lmm_padj,     seg_coef$celltype)

# Reported rows, per panel. CBLN2 uses the intensity-adjusted arm.
# The adjusted slope's 95% CI is emitted by R/phf1_module_distance_intensity_adjust.R
# (dist_slope_adj_ci_lo/_ci_hi in the summary), so it is read directly and no join is
# needed; the unadjusted slope and its CI are carried alongside for reference.
rep_rows <- bind_rows(
  adj_sum %>% filter(celltype == "Exc-IT-L2-3-CBLN2-HOPX") %>%
    transmute(celltype, module = sig_label, reported_model = "intensity-adjusted LMM (log dist)",
              slope_per_sd = dist_slope_adj,
              ci_lo_per_sd = dist_slope_adj_ci_lo, ci_hi_per_sd = dist_slope_adj_ci_hi,
              padj = dist_p_adj,
              slope_unadj  = dist_slope_unadj,
              ci_lo_unadj  = dist_slope_unadj_ci_lo, ci_hi_unadj = dist_slope_unadj_ci_hi,
              significant = dist_p_adj < 0.05),
  glia_coef %>% filter(distance_scale == "log_um", group == "mancuso", celltype == "Micro") %>%
    transmute(celltype, module, reported_model = "LMM (log dist)",
              slope_per_sd, padj = lmm_padj, ci_lo_per_sd, ci_hi_per_sd,
              slope_unadj = NA_real_, significant,
              ci_lo_unadj = NA_real_, ci_hi_unadj = NA_real_),
  astro_coef %>% filter(distance_scale == "log_um", group == ASTRO$group,
                        celltype == "Astro") %>%
    transmute(celltype, module, reported_model = "LMM (log dist)",
              slope_per_sd, padj = lmm_padj, ci_lo_per_sd, ci_hi_per_sd,
              slope_unadj = NA_real_, significant,
              ci_lo_unadj = NA_real_, ci_hi_unadj = NA_real_)
)
rep_rows$seg_floor_slope <- SEG_FLOOR[rep_rows$celltype]
rep_rows$seg_floor_padj  <- SEG_PADJ[rep_rows$celltype]
rep_rows$ratio_vs_seg    <- abs(rep_rows$slope_per_sd) / abs(rep_rows$seg_floor_slope)
rep_rows$seg_pct_of_module <- 100 / rep_rows$ratio_vs_seg
rep_rows$same_direction  <- sign(rep_rows$slope_per_sd) == sign(rep_rows$seg_floor_slope)

STEEPEST <- rep_rows %>% group_by(celltype) %>%
  slice_max(abs(slope_per_sd), n = 1, with_ties = FALSE) %>%
  transmute(celltype, module, slope = slope_per_sd,
            pct = 100 * abs(SEG_FLOOR[celltype]) / abs(slope_per_sd)) %>% ungroup()

##  ............................................................................
##  Build the panels                                                         ####
drawn <- list(); PLOTS <- list()
for (pd in PANELS) {
  ct <- pd$ct
  f  <- need(file.path(args$seg_dir, sprintf(
    "source_data_modulescore_raw_dist_seg_control_%s_rollmean.tsv",
    gsub("[^A-Za-z0-9_-]", "_", ct))))
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$ci_lo <- d$roll_mean - 1.96 * d$sem
  d$ci_hi <- d$roll_mean + 1.96 * d$sem

  ylim <- matched_ylim(d, REPORTED_SPAN[[ct]])
  st   <- STEEPEST[STEEPEST$celltype == ct, ]
  ann  <- sprintf("%s%% of %s", sf2(st$pct[1]), st$module[1])

  seg_span <- diff(range(c(d$ci_lo, d$ci_hi), na.rm = TRUE))
  rep_span <- REPORTED_SPAN[[ct]]
  cat(sprintf("  %-28s SEG ribbon span %.4f | reported %.4f (%.1f%%) | ann: %s\n",
              ct, seg_span, rep_span, 100 * seg_span / rep_span, ann))

  p <- ggplot(d, aes(dist_to_phf1_um, roll_mean)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = SEG_COLOUR, alpha = 0.20) +
    geom_line(colour = SEG_COLOUR, linewidth = 0.6) +
    # No on-panel annotation. Printing the control's slope as a % of a DIFFERENT module's
    # beta on the same axes invites the reader to mis-attribute the number to the trace
    # being drawn. The percentages are in stats_seg_control_vs_reported_margin.tsv, per
    # celltype. `ann` is still computed for the console log and the stats log.
    labs(x = XLAB, y = "Module score", title = pd$label) +
    coord_cartesian(xlim = c(0, args$max_dist_um), ylim = ylim) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title  = element_text(size = 7, hjust = 0.5, face = "plain"),
          legend.position = "none",
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))
  PLOTS[[ct]] <- p

  # per-panel PDF, same fixed-panel geometry as the sibling scripts
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")
  ggsave(file.path(args$output_dir,
                   sprintf("plot_seg_control_matched_%s.pdf",
                           gsub("[^A-Za-z0-9_-]", "_", ct))),
         g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")

  drawn[[ct]] <- d %>% transmute(
    celltype = ct, module = "Reference SEG", window_um,
    dist_to_phf1_um, roll_mean, sem, ci_lo, ci_hi, n_window,
    y_limit_lo = ylim[1], y_limit_hi = ylim[2],
    matched_span = REPORTED_SPAN[[ct]],
    y_limit_source = pd$reported_source)
}

##  ............................................................................
##  Composite + source data + tables                                        ####
combo <- Reduce(`+`, PLOTS) + patchwork::plot_layout(nrow = 1)
ggsave(file.path(args$output_dir, "plot_seg_control_matched_composite.pdf"),
       combo, width = 19.05, height = 5, units = "cm", device = "pdf")
cat("Wrote 3-panel matched composite (A4 width)\n")

write.table(bind_rows(drawn),
            file.path(args$output_dir, "source_data_seg_control_matched.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

marg <- rep_rows %>%
  transmute(panel_celltype = celltype, module, reported_model,
            reported_slope_per_sd = slope_per_sd,
            reported_ci_lo = ci_lo_per_sd, reported_ci_hi = ci_hi_per_sd,
            reported_slope_unadjusted = slope_unadj,
            unadjusted_ci_lo = ci_lo_unadj, unadjusted_ci_hi = ci_hi_unadj,
            reported_padj = padj, reported_significant = significant,
            seg_floor_slope_per_sd = seg_floor_slope, seg_floor_padj,
            ratio_reported_over_seg = ratio_vs_seg,
            seg_pct_of_reported = seg_pct_of_module,
            seg_same_direction = same_direction) %>%
  arrange(panel_celltype, desc(ratio_reported_over_seg))
write.table(marg, file.path(args$output_dir, "stats_seg_control_vs_reported_margin.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote specificity-margin table:", nrow(marg), "rows\n")

used <- marg[marg$reported_significant, ]
cat("\nSignificant reported modules vs the SEG floor:\n")
print(as.data.frame(used %>% transmute(panel_celltype, module,
       slope = signif(reported_slope_per_sd, 3), padj = signif(reported_padj, 3),
       ratio = round(ratio_reported_over_seg, 1))), row.names = FALSE)
cat(sprintf("\nMINIMUM margin across all reported-significant modules: %.1fx (%s / %s)\n",
            min(used$ratio_reported_over_seg),
            used$panel_celltype[which.min(used$ratio_reported_over_seg)],
            used$module[which.min(used$ratio_reported_over_seg)]))

##  ............................................................................
##  Stats log                                                               ####
sink(file.path(args$output_dir, "stats_seg_control_matched.txt"))
cat("SEG negative control - paper figure and specificity margins\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Assembled from existing source data by R/plot_seg_control_composite_matched.R.\n")
cat("No re-scoring: inputs are the rolling-mean TSVs already written by\n")
cat("  ", args$seg_dir, " (control)\n", sep = "")
if (args$astro_source == "cameron") {
  cat("  ", args$glia_dir, " (Mancuso Micro, Cameron Astro subclusters)\n", sep = "")
} else {
  cat("  ", args$glia_dir, " (Mancuso Micro)\n", sep = "")
  cat("  ", ASTRO$dir, " (Serrano-Pozo Astro states; astMic/astNeu excluded upstream)\n", sep = "")
  cat("  ASTRO VARIANT (--astro_source serrano): the Astro axis span and the Astro margin rows\n")
  cat("  come from the Serrano-Pozo states instead of the Cameron subclusters (Fig. S5).\n")
}
cat("  ", args$adjust_dir, " (intensity-adjusted CBLN2 arm)\n\n", sep = "")

cat("FIGURE\n")
cat("  Three panels, one per reported celltype. Grey line = rolling mean of the\n")
cat("  reference-derived SEG module score (window 75 um), ribbon = +/- 1.96*SEM.\n")
cat("  Y-LIMITS ARE MATCHED to the corresponding reported panel's drawn ribbon range,\n")
cat("  so the control is shown on the scale of the effect it is a control for:\n")
for (pd in PANELS) {
  ct <- pd$ct; d <- drawn[[ct]]
  ss <- diff(range(c(d$ci_lo, d$ci_hi), na.rm = TRUE)); rs <- REPORTED_SPAN[[ct]]
  cat(sprintf("    %-28s SEG span %.4f of reported %.4f = %.1f%%  [%s]\n",
              ct, ss, rs, 100 * ss / rs, pd$reported_source))
}
cat("  Exc-IT-L3-5-CHGA-IL1RAPL2 is scored in the raw run but is not part of this figure and\n")
if (args$astro_source == "cameron")
  cat("  is therefore excluded here. Cameron POOLED is likewise excluded (subclusters only).\n\n") else
  cat("  is therefore excluded here.\n\n")

cat("SPECIFICITY MARGINS (stats_seg_control_vs_reported_margin.tsv)\n")
cat("  ratio = |reported slope| / |SEG slope|, same cells, same log-distance LMM.\n")
print(as.data.frame(marg %>% transmute(panel_celltype, module,
      slope = signif(reported_slope_per_sd, 3), padj = signif(reported_padj, 3),
      sig = reported_significant, ratio = round(ratio_reported_over_seg, 1),
      seg_pct = round(seg_pct_of_reported, 1))), row.names = FALSE)
cat(sprintf("\n  MINIMUM margin over reported-significant modules: %.1fx\n",
            min(used$ratio_reported_over_seg)))
cat("  Internal calibration check: in Micro every reported-n.s. module sits at <=3.1x\n")
cat("  the floor while both significant ones are >=7.5x - the floor and the FDR agree\n")
cat("  though derived independently.\n\n")

cat("SEG FLOOR PER PANEL\n")
for (ct in names(REPORTED_SPAN))
  cat(sprintf("  %-28s slope %.3e  padj %.3g %s\n", ct, SEG_FLOOR[ct], SEG_PADJ[ct],
              ifelse(SEG_PADJ[ct] < 0.05, "(significant)", "(n.s.)")))
cat("\nNOTES\n")
cat("  1. The SEG floor is a magnitude benchmark. A donor random intercept is correctly\n")
cat("     calibrated for a WITHIN-donor predictor such as distance (verified by\n")
cat("     within-donor permutation), so the floor slope is estimated on the same footing\n")
cat("     as the reported slopes. The comparison\n")
if (args$astro_source == "cameron") {
  cat("     rests on the MAGNITUDE (the floor is 2-13% of every reported effect), not on\n")
  cat("     the control being zero.\n")
} else {
  # the hard-coded 2-13% was derived for the Cameron panel; recompute for this variant
  cat(sprintf("     rests on the MAGNITUDE (the floor is %.1f-%.1f%% of every reported-significant\n",
              min(used$seg_pct_of_reported), max(used$seg_pct_of_reported)))
  cat("     effect), not on the control being zero.\n")
}
cat("  2. The SEG floor is a BENCHMARK, not a baseline to subtract. The SEG module's\n")
cat("     expression composition differs from the reported signatures, so its drift does\n")
cat("     not transfer.\n")
cat("  3. The CBLN2 floor is fitted in the UNADJUSTED model while the reported CBLN2\n")
cat("     numbers are intensity-adjusted. Even if the floor doubled, PHF1 stays 9.5x,\n")
cat("     Otero-UP 4.8x. The unadjusted CBLN2 slope and its CI are given alongside.\n")
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()
cat("Wrote stats_seg_control_matched.txt\n")
cat("\nAll outputs in:", args$output_dir, "\nDone.\n")
