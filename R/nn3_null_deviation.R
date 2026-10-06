#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# nn3_null_deviation.R
#
# Figure panels: S3C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# nn3_null_deviation.R
#
# HOW FAR THE OBSERVED nn3-DISTANCE COEFFICIENT SITS FROM ITS GEOMETRIC NULL.
#
# This is a PURE RE-PLOT of published null output. It fits NOTHING: every value is read
# from the two tables R/null_replica_nn3.R already wrote, so this figure cannot disagree
# with the null replica, and re-running it after the null is regenerated is the only way its
# numbers ever change.
#   plots/nn3_neuron_spacing/source_data_nn3_null_distance.tsv    1000 permutation betas x 10 subtypes
#   plots/nn3_neuron_spacing/stats_nn3_null_replica_summary.tsv   observed, null median/quantiles, excess, padj_emp
#
# WHY THIS FIGURE EXISTS. dist_to_phf1_um is itself a nearest-neighbour distance over the
# same neuron point pattern that nn3 measures, so the raw coefficient has a geometric
# component. The quantity of interest is the DEVIATION FROM THE PERMUTATION NULL, which this
# plot shows.
#
# TWO PANELS, WRITTEN AS TWO FIGURES:
#   plot_nn3_null_deviation_CBLN2.pdf  the 1000 null coefficients for Exc-IT-L2-3-CBLN2-HOPX
#                                      as a histogram, with the observed value marked and the
#                                      excess drawn as the gap from the null median. A second
#                                      x axis restates the coefficient as a % change in
#                                      spacing over a tenfold change in distance (50 vs
#                                      500 um).
#   plot_nn3_null_deviation_all.pdf    every subtype's excess with the null's OWN spread as
#                                      the error bar, so "outside the null" is read directly.
#                                      The two sparse subtypes (SPON1-FGD4 n=72,
#                                      CTXN1-ERC2 n=282) have wide nulls; they are drawn,
#                                      not hidden.
#
# THE SECOND X AXIS IS A RESCALING, NOT A SECOND ESTIMATE. 1 s.d. of log-distance is
# DIST_SD = 0.76977 within this subtype (from source_data_nn3_distance_coef.tsv), so a
# tenfold change in distance spans log(10)/DIST_SD = 2.9913 s.d. and a coefficient b maps to
# 100*(exp(b*2.9913)-1) %. Monotone, so the axis is well defined.
#
# The observed value plotted is beta_obs_fast, the within-donor projection estimate the
# permutations themselves are computed with -- NOT the lmer coefficient (-0.02958). The two
# agree to 3.6e-4 and null_replica_nn3.R asserts that before reporting, but mixing them here
# would put the observed marker on a different scale from the histogram it sits in.
#
# Run: Rscript R/nn3_null_deviation.R      (seconds; no model is fitted)
# Outputs the standard TRIPLE under plots/nn3_null_deviation/.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/nn3_null_hist_panel.R")   # nn3_null_hist(): ONE definition, drawn here and as a
                                   # panel of R/nn3_dist_phf1_panel_null.R

FOCUS   <- "Exc-IT-L2-3-CBLN2-HOPX"
FDR     <- 0.05
OBS_COL <- "#BD0026"    # the PHF1 red used across the nn3 family for the observed value
NULL_COL <- "grey78"
D_LO    <- 50           # the two distances used for the second x axis
D_HI    <- 500

src_dir <- "plots/nn3_neuron_spacing"
perm_f  <- file.path(src_dir, "source_data_nn3_null_distance.tsv")
sum_f   <- file.path(src_dir, "stats_nn3_null_replica_summary.tsv")
coef_f  <- file.path(src_dir, "source_data_nn3_distance_coef.tsv")
for (f in c(perm_f, sum_f, coef_f))
  if (!file.exists(f)) stop("Missing: ", f, "\nRun R/null_replica_nn3.R first.")

out_dir <- "plots/nn3_null_deviation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

perm <- read.delim(perm_f, sep = "\t", check.names = FALSE)
sm   <- read.delim(sum_f,  sep = "\t", check.names = FALSE) %>% filter(analysis == "B_distance")
cf   <- read.delim(coef_f, sep = "\t", check.names = FALSE)
stopifnot(nrow(sm) > 0, FOCUS %in% sm$celltype, FOCUS %in% perm$celltype)

# The observed column is repeated on every permutation row; it must be constant per subtype
# or the file has been edited by hand.
chk <- perm %>% group_by(celltype) %>% summarise(k = n_distinct(beta_closer_obs_fast))
if (any(chk$k != 1L)) stop("beta_closer_obs_fast is not constant within a subtype in ", perm_f)

DIST_SD <- cf$dist_sd[cf$celltype == FOCUS]
stopifnot(length(DIST_SD) == 1L, is.finite(DIST_SD))
K_SPAN  <- log(D_HI / D_LO) / DIST_SD            # s.d. of log-distance spanned by D_LO..D_HI
to_pct  <- function(b) 100 * (exp(b * K_SPAN) - 1)

## ---------------------------------------------------------------------------
## Panel 1 -- the focus subtype's null distribution with the observed value in it
## ---------------------------------------------------------------------------
pf <- perm %>% filter(celltype == FOCUS)
sf <- sm   %>% filter(celltype == FOCUS)
obs <- sf$beta_obs_fast; nmed <- sf$null_median
n_as_shallow <- sum(pf$beta_closer_null >= obs)   # permutations at least as shallow as observed

# n_as_shallow is computed and reported in the stats log; it is not drawn on the figure.
p1 <- nn3_null_hist(beta_null = pf$beta_closer_null, beta_obs = obs,
                    null_median = nmed, null_q025 = sf$null_q025, null_q975 = sf$null_q975,
                    dist_sd = DIST_SD, d_lo = D_LO, d_hi = D_HI,
                    obs_col = OBS_COL, null_col = NULL_COL)
ggsave(file.path(out_dir, "plot_nn3_null_deviation_CBLN2.pdf"), p1,
       width = 9.6, height = 6.4, units = "cm", device = "pdf")

## ---------------------------------------------------------------------------
## Panel 2 -- every subtype's excess against the null's own spread
## ---------------------------------------------------------------------------
# Centring the null's 2.5-97.5% range on zero turns it into the width of the noise the excess
# has to clear: a point outside its own bar is a subtype whose observed coefficient falls
# outside the central 95% of its permutations.
lv <- c(neuron_order, "Unassigned Neuron")
ex <- sm %>%
  transmute(celltype = factor(celltype, levels = rev(lv)),
            excess = excess_over_null,
            null_lo = null_q025 - null_median, null_hi = null_q975 - null_median,
            n_cells, padj_emp,
            sig = ifelse(!is.na(padj_emp) & padj_emp < FDR, "padj < 0.05", "n.s.")) %>%
  arrange(celltype)

p2 <- ggplot(ex, aes(excess, celltype)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = null_lo, xmax = null_hi), orientation = "y", width = 0,
                linewidth = 2.2, colour = NULL_COL) +
  geom_point(aes(colour = sig), size = 1.7) +
  scale_colour_manual(values = c("padj < 0.05" = OBS_COL, "n.s." = "grey55"), name = NULL) +
  # Short axis title: the full form ("... log um per s.d. closer to a PHF1+ neuron") is
  # wider than the figure and was clipped at the canvas edge. mu via plotmath.
  labs(x = expression("Excess over the geometric null (log " * mu * "m)"), y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top",
        plot.margin = margin(4, 10, 4, 4, unit = "pt"))
ggsave(file.path(out_dir, "plot_nn3_null_deviation_all.pdf"), p2,
       width = 10.4, height = 6.0, units = "cm", device = "pdf")

## ---------------------------------------------------------------------------
## Source data + stats log
## ---------------------------------------------------------------------------
wr(pf %>% transmute(celltype, perm, beta_closer_null,
                    beta_closer_obs_fast, pct_change_null = to_pct(beta_closer_null)),
   "source_data_nn3_null_deviation_CBLN2.tsv")
wr(sm %>% transmute(celltype, n_cells, beta_obs_fast, null_median, null_q025, null_q975,
                    excess_over_null, inside_null_ci, p_emp, padj_emp,
                    pct_obs = to_pct(beta_obs_fast), pct_null = to_pct(null_median),
                    pct_excess = 100 * (exp(excess_over_null * K_SPAN) - 1)),
   "source_data_nn3_null_deviation_all.tsv")

sink(file.path(out_dir, "stats_nn3_null_deviation.txt"))
ttl <- "Observed nn3-distance coefficient vs its geometric permutation null"
cat(ttl, "\n"); cat(strrep("=", nchar(ttl)), "\n"); cat("Date:", format(Sys.time()), "\n\n")
cat("NOTHING IS FITTED HERE. Every value is read from:\n")
cat("  ", perm_f, "\n  ", sum_f, "\n  ", coef_f, " (dist_sd only)\n", sep = "")
cat("written by R/null_replica_nn3.R and R/nn3_over_phf1_distance.R.\n\n")
cat("WHY: dist_to_phf1_um is itself a nearest-neighbour distance over the same neuron point\n")
cat("pattern nn3 measures, so the raw coefficient has a geometric component. The quantity of\n")
cat("interest is the deviation from the null.\n\n")
cat("NULL CONSTRUCTION (from R/generate_phf1_null_labels.r, seed 42): 1000 labellings, PHF1\n")
cat("shuffled within celltype x sample with per-group PHF1+ counts held fixed. nn3_um is\n")
cat("invariant to PHF1 relabelling, so only the regressor changes.\n\n")
cat("THE NULL'S 2.5-97.5% RANGE IS THE SPREAD OF THE PERMUTATION DISTRIBUTION.\n")
cat("  It is NOT a confidence interval and must not be described as one.\n\n")
cat(sprintf("SECOND X AXIS: dist_sd = %.5f within %s, so %d vs %d um spans %.4f s.d. of\n",
            DIST_SD, FOCUS, D_LO, D_HI, K_SPAN))
cat(sprintf("  log-distance and b maps to 100*(exp(b*%.4f)-1) %%. A rescaling of the same\n", K_SPAN))
cat("  coefficient, not a second estimate.\n\n")
cat(sprintf("FOCUS SUBTYPE: %s\n", FOCUS))
cat(sprintf("  observed (projection est.)   %+.5f   -> %+.2f%% over %d vs %d um\n",
            obs, to_pct(obs), D_LO, D_HI))
cat(sprintf("  null median                  %+.5f   -> %+.2f%%\n", nmed, to_pct(nmed)))
cat(sprintf("  null central 95%%             %+.5f to %+.5f\n", sf$null_q025, sf$null_q975))
cat(sprintf("  excess over null             %+.5f   -> %+.2f%% wider than geometry predicts\n",
            sf$excess_over_null, 100 * (exp(sf$excess_over_null * K_SPAN) - 1)))
cat(sprintf("  permutations >= observed     %d of %d\n", n_as_shallow, nrow(pf)))
cat(sprintf("  empirical p %.4g | BH-adjusted %.4g | inside null 95%%: %s\n",
            sf$p_emp, sf$padj_emp, sf$inside_null_ci))
cat("  NOTE the observed value is beta_obs_fast, the within-donor projection estimate the\n")
cat("  permutations use; the lmer coefficient is -0.02958 and the two agree to 3.6e-4.\n")
cat("  Do not mix them: the histogram is on the projection scale.\n")

cat("\n=== ALL SUBTYPES ===\n")
print(as.data.frame(sm %>%
        transmute(celltype, n_cells, beta_obs_fast, null_median, excess_over_null,
                  inside_null_ci, p_emp, padj_emp)),
      row.names = FALSE, digits = 4)

cat("\nNOTES:\n")
cat(" - Exc-ET-L5-SPON1-FGD4 (n=72) and Exc-IT-L6-CTXN1-ERC2 (n=282) have very wide nulls;\n")
cat("   they are drawn rather than hidden.\n")
cat(" - The empirical p is two-tailed against the null's OWN distribution, not against zero,\n")
cat("   and uses the (1+r)/(1+n) convention, so it is bounded below by 2/(1+1000).\n")
cat("\n=== sessionInfo() ===\n"); print(sessionInfo())
sink()

cat("Wrote", length(list.files(out_dir)), "files to", out_dir, "\n")
cat(sprintf("%s: observed %+.5f vs null %+.5f, excess %+.5f (%+.2f%% over %d vs %d um), %d/%d perms, padj_emp %.4g\n",
            FOCUS, obs, nmed, sf$excess_over_null,
            100 * (exp(sf$excess_over_null * K_SPAN) - 1), D_LO, D_HI,
            n_as_shallow, nrow(pf), sf$padj_emp))
