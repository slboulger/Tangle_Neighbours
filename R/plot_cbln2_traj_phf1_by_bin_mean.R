#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_cbln2_traj_phf1_by_bin_mean.R
#
# Figure panels: 2C
# Reads the per-bin table written by cbln2_slingshot_trajectory.R.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_cbln2_traj_phf1_by_bin_mean.R
#
# % tangle-bearing (PHF1+) CBLN2 neurons in each equal-width pseudotime bin, drawn as the DONOR
# MEAN with its 95% confidence band (the per-donor lines are plot_cbln2_traj_phf1_by_bin.pdf,
# written by R/cbln2_slingshot_trajectory.R).
#
# Donor is the unit: for each bin, the mean of the per-donor % PHF1+ over the donors with
# >= min_pos_donor PHF1+ cells (the same donors as the trajectory script's per-donor lines),
# with a t-based 95% CI (df = n donors - 1). The band is clipped at 0 for display only (a
# percentage cannot be negative); the unclipped limits are in the source data and stats.
#
# Reads plots/cbln2_slingshot_trajectory/source_data_cbln2_traj_phf1_by_bin.tsv (written by the
# trajectory script), so it needs no refit:
#   Rscript R/plot_cbln2_traj_phf1_by_bin_mean.R                              # CBLN2
#   Rscript R/plot_cbln2_traj_phf1_by_bin_mean.R Exc-IT-L3-5-CHGA-IL1RAPL2    # another celltype
# Annotation: a bar over bins 1..10 with the Friedman p across bins (unit = donor x bin %, bins
# within donor); Kruskal-Wallis and the selection-adjusted slope p are in the stats log only.
# OUTPUT (same folder): plot_/source_data_/stats_cbln2_traj_phf1_by_bin_mean.*

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # fig_theme
source("R/cbln2_embedding_utils.R")  # ct_tag()
source("R/fixed_panel_utils.R")      # BIN_PANEL_W_CM: same x-axis length as the Otero bin plot

CT   <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else CBLN2_CT
TAG  <- ct_tag(CT)
DIR  <- file.path("plots", paste0(TAG, "_slingshot_trajectory"))
NAME <- paste0(TAG, "_traj_phf1_by_bin_mean")
src  <- read.delim(file.path(DIR, paste0("source_data_", TAG, "_traj_phf1_by_bin.tsv")),
                   colClasses = c(sample_id = "character"))
donor <- src %>% filter(level == "donor")
n_bins <- max(src$bin)

summ <- donor %>% group_by(bin) %>%
  summarise(n_donors = n(), mean_pct = mean(pct_phf1), sd_pct = sd(pct_phf1),
            se = sd_pct / sqrt(n_donors), tcrit = qt(0.975, n_donors - 1),
            ci_lo = mean_pct - tcrit * se, ci_hi = mean_pct + tcrit * se,
            n_cells = sum(n), n_phf1 = sum(n_phf1), .groups = "drop") %>%
  mutate(ci_lo_plot = pmax(ci_lo, 0))
stopifnot(all(summ$n_donors == n_distinct(donor$sample_id)))   # every donor in every bin

# Selection-adjusted permutation p for the rise (R/cbln2_trajectory_selection_null.R), if run.
null_file <- file.path("plots", paste0(TAG, "_trajectory_optimise"), paste0("stats_", TAG, "_trajopt_null.tsv"))
null_res  <- if (file.exists(null_file)) read.delim(null_file) else NULL
p_lab <- NULL
if (!is.null(null_res) && "slope" %in% null_res$statistic) {
  r <- null_res[null_res$statistic == "slope", ]
  p_lab <- if (r$p_selection_adjusted <= 1 / (1 + r$n_perm)) sprintf("p < %.3g", 1 / r$n_perm) else sprintf("p = %.2g", r$p_selection_adjusted)
}

# Tests across bins: unit = donor x bin %. Friedman (bins within donor) is annotated;
# Kruskal-Wallis (treats each donor's 10 bins as independent) is kept in the log.
kw <- kruskal.test(pct_phf1 ~ factor(bin), data = donor)
kw_eps2 <- unname(kw$statistic) / (nrow(donor) - 1)               # epsilon^2 = H / (n - 1)
fr <- friedman.test(pct_phf1 ~ bin | sample_id, data = donor)
# Annotated test = Friedman, label just "p = " (test named in the caption); KW in the log.
# plotmath label: small p as "1.6 x 10^-7", otherwise plain; drawn with parse = TRUE.
p_plotmath <- function(p) {
  if (p >= 1e-3) return(sprintf('p == "%s"', formatC(p, format = "g", digits = 2)))
  m <- strsplit(formatC(p, format = "e", digits = 1), "e")[[1]]
  sprintf('p == "%s" %%*%% 10^%d', m[1], as.integer(m[2]))
}
kw_txt <- p_plotmath(fr$p.value)
y_top <- max(summ$ci_hi); yb <- y_top * 1.05; y_lim <- c(0, y_top * 1.16)

p <- ggplot(summ, aes(bin, mean_pct)) +
  geom_ribbon(aes(ymin = ci_lo_plot, ymax = ci_hi), fill = "grey75", alpha = 0.6) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  annotate("segment", x = 1, xend = n_bins, y = yb, yend = yb, linewidth = 0.3) +
  # no end ticks: a plain bar spans all bins, so it does not read as a first-vs-last contrast
  annotate("text", x = (1 + n_bins) / 2, y = yb, label = kw_txt, parse = TRUE, vjust = -0.4, size = 2.1) +
  scale_x_continuous(breaks = seq_len(n_bins)) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0))) +
  labs(x = "Pseudotime bin", y = "% Tangle-bearing") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        plot.margin = margin(4, 4, 4, 4))
ggsave(file.path(DIR, paste0("plot_", NAME, ".pdf")), p, width = 7, height = 5.5, units = "cm", device = "pdf")
# Square figure-panel version: 5.08 x 5.08 cm (2 x 2 in), same plot and source data.
message(sprintf("Natural panel width at 5.08 cm: %.3f cm (fixed to BIN_PANEL_W_CM = %.2f); left %.3f right %.3f", panel_width_cm(p, 5.08), BIN_PANEL_W_CM, side_widths_cm(p)[["left"]], side_widths_cm(p)[["right"]]))
save_fixed_panel(p, file.path(DIR, paste0("plot_", NAME, "_5.08cm.pdf")), BIN_PANEL_W_CM, 5.08, 5.08)
write.table(summ %>% select(-tcrit), file.path(DIR, paste0("source_data_", NAME, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Effect size: first vs last bin, donor-paired (each donor contributes to both bins).
first_last <- donor %>% filter(bin %in% c(1, n_bins)) %>%
  select(sample_id, bin, pct_phf1) %>% tidyr::pivot_wider(names_from = bin, values_from = pct_phf1)
dd <- first_last[[as.character(n_bins)]] - first_last[["1"]]
tt <- t.test(dd)

sink(file.path(DIR, paste0("stats_", NAME, ".txt")))
cat("% tangle-bearing", CT, "cells per pseudotime bin: donor mean with 95% t CI\n", strrep("=", 72), "\n", sep = "")
cat("Run:", format(Sys.time()), "\nInput:", file.path(DIR, paste0("source_data_", TAG, "_traj_phf1_by_bin.tsv")), "\n")
cat("Donors (as used by the trajectory script):", paste(sort(unique(donor$sample_id)), collapse = ", "),
    sprintf("(n = %d; CI df = %d)\n", n_distinct(donor$sample_id), n_distinct(donor$sample_id) - 1))
cat("Band clipped at 0 for display only; unclipped limits below.\n\n")
print(as.data.frame(summ %>% select(bin, n_donors, n_cells, n_phf1, mean_pct, sd_pct, ci_lo, ci_hi)),
      digits = 3, row.names = FALSE)
cat(sprintf("\nKruskal-Wallis across the %d pseudotime bins (log only; unit = donor x bin %%, n = %d):\n  H = %.2f, df = %d, p = %.3g, epsilon^2 = %.2f\n",
            n_bins, nrow(donor), unname(kw$statistic), unname(kw$parameter), kw$p.value, kw_eps2))
cat(sprintf("  Friedman, bins within donor (THE p ANNOTATED ON THE BIN CURVE; KW ignores that each donor contributes every bin):\n  chi2 = %.2f, df = %d, p = %.3g, Kendall's W = %.2f\n",
            unname(fr$statistic), unname(fr$parameter), fr$p.value, unname(fr$statistic) / (n_distinct(donor$sample_id) * (n_bins - 1))))
cat("  The selection-adjusted slope p below repeats the trajectory selection on every permutation.\n")
cat(sprintf("\nEffect size, last (bin %d) minus first (bin 1), donor-paired: %+.2f percentage points [%.2f, %.2f], dz = %.2f, %d/%d donors higher\n",
            n_bins, mean(dd), tt$conf.int[1], tt$conf.int[2], mean(dd) / sd(dd), sum(dd > 0), length(dd)))
if (!is.null(null_res)) {
  cat("\nSelection-adjusted permutation test (R/cbln2_trajectory_selection_null.R; whole trajectory search re-run on\n")
  cat("PHF1 labels shuffled within donor):\n")
  print(null_res[, c("statistic", "observed", "null_mean", "null_q95", "p_selection_adjusted", "n_perm", "n_donors")],
        digits = 4, row.names = FALSE)
  cat("The plot carries the Friedman p, not this 'slope' p.\n")
} else cat("No selection-adjusted p yet: run R/cbln2_trajectory_selection_null.R --celltype", CT, "\n")
cat("\n")
print(sessionInfo())
sink()
message("Written: ", file.path(DIR, paste0("plot_", NAME, ".pdf")))
