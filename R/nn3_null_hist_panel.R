#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# nn3_null_hist_panel.R
#
# Shared utility - sourced by the scripts above, no panel of its own
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# nn3_null_hist_panel.R
#
# SHARED BUILDER for the permutation-null histogram: the 1000 null coefficients for one
# subtype with the observed value marked and the excess drawn as the gap from the null median.
#
# It lives in its own file because TWO scripts draw it -- R/nn3_null_deviation.R (standalone)
# and R/nn3_dist_phf1_panel_null.R (as a panel of the assembled figure) -- and a figure that
# appears twice must have ONE definition. Duplicating the ggplot code is how the two copies
# start disagreeing after the first edit to either.
#
# Fits nothing. Every value is passed in, having been read from the published null tables.

nn3_null_hist <- function(beta_null, beta_obs, null_median, null_q025, null_q975,
                          dist_sd, d_lo = 50, d_hi = 500,
                          obs_col = "#BD0026", null_col = "grey78", nbin = 40) {
  stopifnot(length(beta_null) > 0, is.finite(beta_obs), is.finite(dist_sd), dist_sd > 0)

  # Second x axis: a RESCALING of the same coefficient, not a second estimate. A tenfold
  # change in distance spans log(10)/dist_sd s.d. of log-distance, so b -> 100*(exp(b*k)-1)%.
  k_span <- log(d_hi / d_lo) / dist_sd
  to_pct <- function(b) 100 * (exp(b * k_span) - 1)

  # ymax from the SAME equal-width binning ggplot uses, so the annotations cannot land on top
  # of the tallest bar (hist()'s pretty() breaks give a different count).
  ymax <- max(table(cut(beta_null, breaks = nbin)))
  xr   <- range(c(beta_null, beta_obs))
  pad  <- 0.06 * diff(xr)

  ggplot2::ggplot(data.frame(b = beta_null), ggplot2::aes(b)) +
    ggplot2::geom_histogram(bins = nbin, fill = null_col, colour = NA) +
    # the null's central 95%. This is the SPREAD OF THE NULL DISTRIBUTION, not a confidence
    # interval, and must never be described as one.
    ggplot2::geom_vline(xintercept = c(null_q025, null_q975), linetype = 3,
                        linewidth = 0.25, colour = "grey45") +
    ggplot2::geom_vline(xintercept = null_median, linetype = 2,
                        linewidth = 0.35, colour = "grey25") +
    ggplot2::geom_vline(xintercept = beta_obs, linewidth = 0.7, colour = obs_col) +
    ggplot2::annotate("segment", x = null_median, xend = beta_obs,
                      y = 1.04 * ymax, yend = 1.04 * ymax, linewidth = 0.4, colour = obs_col,
                      arrow = grid::arrow(length = grid::unit(3, "pt"), type = "closed")) +
    ggplot2::annotate("text", x = mean(c(null_median, beta_obs)), y = 1.10 * ymax,
                      size = 2.2, colour = obs_col,
                      label = sprintf("excess %+.3f", beta_obs - null_median)) +
    ggplot2::annotate("text", x = beta_obs, y = 0.5 * ymax, angle = 90, vjust = -0.5,
                      size = 2.3, fontface = "bold", colour = obs_col, label = "observed") +
    ggplot2::annotate("text", x = null_median, y = 0.5 * ymax, angle = 90, vjust = 1.3,
                      size = 2.2, colour = "grey25", label = "null median") +
    ggplot2::scale_x_continuous(
      limits = c(xr[1] - pad, xr[2] + pad),
      sec.axis = ggplot2::sec_axis(~ to_pct(.), breaks = c(-16, -14, -12, -10, -8),
        name = sprintf("Change in spacing, %d vs %d um (%%)", d_lo, d_hi))) +
    # headroom for the arrow and its label only
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.18))) +
    ggplot2::labs(x = expression("Coefficient per s.d. closer to a PHF1+ neuron (log " * mu * "m)"),
                  y = "Permutations") +
    ggplot2::theme_classic(base_size = 8) + fig_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                   plot.margin = ggplot2::margin(t = 4, r = 8, b = 4, l = 5))
}
