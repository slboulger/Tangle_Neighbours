# ---------------------------------------------------------------------------
# imc_tangle_depth_deciles_plot_utils.R
#
# Figure panels: 4C, 4D - helper
# Sourced by imc_tangle_depth_density_shift.R; plot_density_shift() draws 4C and 4D.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_tangle_depth_deciles_plot_utils.R
#
# Plot helpers for the IMC depth-decile panels.
#
#   plot_density_shift() -- one point per depth decile: x = % difference in columnar density,
#     AD vs control; y = the PHF1 measure. Used for 4C and 4D.
#   plot_tangle_vs_density_shift() / save_tangle_vs_density_shift() /
#   replot_tangle_vs_density_shift() -- an alternative layout (x = % PHF1+ neurons, y = neuron
#     share ratio) drawn from a source_data_tangle_vs_density_shift.tsv table; not used for any
#     figure panel.
#
# Design follows plot_cluster_loss_vs_tangle_distance_raw.pdf (imc_cluster_loss_vs_tangle_
# distance.R, mk()): outlined crosshairs (wider grey25 segment under the fill colour, square caps),
# shape-21 points with a thin dark outline, legend on the right, no in-panel text.
#
# Colour: one sequential colour per decile ordered by its PHF1 RANK, so the most tangle-laden decile
# is the most severe colour. viridis family, inferno, in its native order so rank 1 is the darkest
# (inferno() runs dark -> pale). begin/end trim the near-black and near-white ends; the outlines keep
# the palest entries visible on white.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

plot_tangle_vs_density_shift <- function(link) {
  link <- link %>% arrange(rank)
  lev  <- sprintf("%s (#%d)", link$depth, as.integer(link$rank))   # legend in rank order
  link$decile_lab <- factor(lev, levels = lev)
  cols <- setNames(viridisLite::inferno(nrow(link), begin = 0.12, end = 0.92), lev)

  ggplot(link, aes(x = 100 * prop, y = share_ratio)) +
    geom_hline(yintercept = 1, linewidth = 0.25, linetype = "dashed", colour = "grey60") +
    # all outlines first so an outline never cuts across another decile's bar
    geom_segment(aes(x = 100 * prop, xend = 100 * prop, y = ratio_low, yend = ratio_high),
                 linewidth = 0.75, colour = "grey25", lineend = "square") +
    geom_segment(aes(x = 100 * prop_low, xend = 100 * prop_high, y = share_ratio, yend = share_ratio),
                 linewidth = 0.75, colour = "grey25", lineend = "square") +
    geom_segment(aes(x = 100 * prop, xend = 100 * prop, y = ratio_low, yend = ratio_high,
                     colour = decile_lab), linewidth = 0.45, lineend = "square") +
    geom_segment(aes(x = 100 * prop_low, xend = 100 * prop_high, y = share_ratio, yend = share_ratio,
                     colour = decile_lab), linewidth = 0.45, lineend = "square") +
    geom_point(aes(fill = decile_lab), shape = 21, size = 2, stroke = 0.3, colour = "grey15") +
    scale_fill_manual(values = cols, name = "Depth from pia\n(PHF1+ rank)", drop = FALSE) +
    scale_colour_manual(values = cols, guide = "none", drop = FALSE) +
    labs(x = "PHF1+ neurons (%, pooled AD)", y = "Neuron share, Braak 6 / Braak 0-1") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
          plot.margin = margin(4, 4, 4, 4), legend.position = "right",
          legend.key.size = grid::unit(8, "pt"), legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5))
}

save_tangle_vs_density_shift <- function(link, out_dir) {
  ggsave(file.path(out_dir, "plot_tangle_vs_density_shift.pdf"), plot_tangle_vs_density_shift(link),
         width = 3.6, height = 2.6, device = "pdf")
}

replot_tangle_vs_density_shift <- function(out_dir) {
  link <- read.delim(file.path(out_dir, "source_data_tangle_vs_density_shift.tsv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
  save_tangle_vs_density_shift(link, out_dir)
}

# -------------------------------------------------------------------
# Used by imc_tangle_depth_density_shift.R (4C, 4D): x = % difference in columnar density,
# AD vs control -- the x axis of plot_cluster_loss_vs_tangle_distance_raw.pdf -- and y = the PHF1
# measure. One point per depth decile, coloured and labelled by its PHF1 rank only (rank 1 = most
# PHF1 = darkest). `d` needs rank, depth, pct_change, pct_lo, pct_hi, y, y_lo, y_hi.
# -------------------------------------------------------------------
plot_density_shift <- function(d, ylab, legend_title, hline = NULL, log_y = FALSE,
                               tags = c("highest", "lowest"),
                               xlab = "Difference in columnar density,\nAD vs control (%)") {
  d <- d %>% arrange(rank)
  # "rank (decile depth from pia)", e.g. "1 (10-20%, highest)"; `depth` is the decile label
  tag <- d$depth
  tag[1] <- paste0(tag[1], ", ", tags[1]); tag[length(tag)] <- paste0(tag[length(tag)], ", ", tags[2])
  lev <- sprintf("%d (%s)", as.integer(d$rank), tag)
  d$rank_lab <- factor(lev, levels = lev)
  cols <- setNames(viridisLite::inferno(nrow(d), begin = 0.12, end = 0.92), lev)

  g <- ggplot(d, aes(pct_change, y)) +
    geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey60")
  if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linewidth = 0.25,
                                           linetype = "dashed", colour = "grey60")
  if (log_y) g <- g + scale_y_log10()   # as plot_cluster_loss_vs_tangle_distance_raw.pdf
  g + geom_segment(aes(x = pct_change, xend = pct_change, y = y_lo, yend = y_hi),
                   linewidth = 0.75, colour = "grey25", lineend = "square") +
    geom_segment(aes(x = pct_lo, xend = pct_hi, y = y, yend = y),
                 linewidth = 0.75, colour = "grey25", lineend = "square") +
    geom_segment(aes(x = pct_change, xend = pct_change, y = y_lo, yend = y_hi, colour = rank_lab),
                 linewidth = 0.45, lineend = "square") +
    geom_segment(aes(x = pct_lo, xend = pct_hi, y = y, yend = y, colour = rank_lab),
                 linewidth = 0.45, lineend = "square") +
    geom_point(aes(fill = rank_lab), shape = 21, size = 2, stroke = 0.3, colour = "grey15") +
    scale_fill_manual(values = cols, name = legend_title, drop = FALSE) +
    scale_colour_manual(values = cols, guide = "none", drop = FALSE) +
    # x title wrapped as in imc_cluster_loss_vs_tangle_distance.R
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
          plot.margin = margin(4, 4, 4, 4), legend.position = "right",
          legend.key.size = grid::unit(8, "pt"), legend.text = element_text(size = 6),
          legend.title = element_text(size = 6.5))
}
