# ---------------------------------------------------------------------------
# fixed_panel_utils.R
#
# Figure panels: 2C, 2D
# Shared utility - sourced by plot_cbln2_traj_phf1_by_bin_mean.R (2C) and
# plot_cbln2_traj_otero_by_bin.R (2D); no panel of its own.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# fixed_panel_utils.R
#
# Give ggplots that sit side by side in a figure the SAME physical panel (x-axis) width, whatever
# their y-axis labels take up. ggsave() sizes the whole canvas, so a wider y-axis label shortens the
# x-axis; fixing the panel cell of the gtable instead makes the axes line up exactly.
#
# Used by plot_cbln2_traj_phf1_by_bin_mean.R and plot_cbln2_traj_otero_by_bin.R so the % tangle and
# Otero-score bin plots share one x-axis length: BIN_PANEL_W_CM.

# Shared x-axis length of the pseudotime-bin plots (cm): the % tangle plot's natural panel width on
# its 5.08 cm canvas (4.200 cm). The % plot keeps its 5.08 cm canvas; the Otero plot's canvas grows
# to fit a panel of this width.
BIN_PANEL_W_CM <- 4.20

# Widths (cm) to the left and right of the panel of a single-panel ggplot.
side_widths_cm <- function(p) {
  pdf(NULL); on.exit(invisible(dev.off()))
  g <- ggplot2::ggplotGrob(p)
  pcol <- g$layout$l[g$layout$name == "panel"]
  cm <- function(i) grid::convertWidth(g$widths[i], "cm", valueOnly = TRUE)
  c(left = sum(vapply(seq_len(pcol - 1), cm, numeric(1))),
    right = sum(vapply(setdiff(seq_along(g$widths), seq_len(pcol)), cm, numeric(1))))
}

# Natural panel width (cm) of a single-panel ggplot drawn on a canvas `width_cm` wide.
panel_width_cm <- function(p, width_cm) {
  pdf(NULL); on.exit(invisible(dev.off()))
  g <- ggplot2::ggplotGrob(p)
  pcol <- g$layout$l[g$layout$name == "panel"]
  other <- sum(vapply(seq_along(g$widths)[-pcol], function(i)
    grid::convertWidth(g$widths[i], "cm", valueOnly = TRUE), numeric(1)))
  width_cm - other
}

# Save `p` with its panel exactly `panel_w_cm` wide starting `left_cm` from the canvas edge (padding
# is added on the left when the plot's own axis/labels are narrower), on a width_cm x height_cm
# canvas; spare width goes to the right. Stops if the plot cannot fit, rather than clipping it.
save_fixed_panel <- function(p, file, panel_w_cm, width_cm, height_cm, left_cm = NULL) {
  pdf(NULL)
  g <- ggplot2::ggplotGrob(p)
  pcol <- g$layout$l[g$layout$name == "panel"]
  cm <- function(i) grid::convertWidth(g$widths[i], "cm", valueOnly = TRUE)
  left  <- sum(vapply(seq_len(pcol - 1), cm, numeric(1)))
  right <- sum(vapply(setdiff(seq_along(g$widths), seq_len(pcol)), cm, numeric(1)))
  invisible(dev.off())
  if (!is.null(left_cm)) {
    if (left > left_cm + 1e-6) stop(sprintf("Left of panel needs %.2f cm > left_cm %.2f", left, left_cm))
    g <- gtable::gtable_add_cols(g, grid::unit(left_cm - left, "cm"), pos = 0)   # left padding
    pcol <- pcol + 1; left <- left_cm
  }
  if (left + panel_w_cm + right > width_cm + 1e-3)
    stop(sprintf("Panel %.2f cm + axes/labels %.2f cm exceed the %.2f cm canvas", panel_w_cm, left + right, width_cm))
  g$widths[pcol] <- grid::unit(panel_w_cm, "cm")
  g <- gtable::gtable_add_cols(g, grid::unit(1, "null"))      # spare width to the right
  ggplot2::ggsave(file, g, width = width_cm, height = height_cm, units = "cm", device = "pdf")
  invisible(c(left = left, right = right))
}

# Save a patchwork whose panel widths are fixed in cm (plot_layout(widths = unit(..., "cm"))) on a
# canvas exactly as wide as it needs, so nothing is centred or padded. Returns the width used (cm).
save_patchwork_exact_width <- function(pw, file, height_cm) {
  pdf(NULL)
  g <- patchwork::patchworkGrob(pw)
  w <- sum(vapply(seq_along(g$widths), function(i) grid::convertWidth(g$widths[i], "cm", valueOnly = TRUE), numeric(1)))
  invisible(dev.off())
  ggplot2::ggsave(file, g, width = w, height = height_cm, units = "cm", device = "pdf")
  w
}
