# ---------------------------------------------------------------------------
# rollmean_by_donor.R
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
# rollmean_by_donor.R
#
# Per-donor sliding-window mean of a module score over distance to the nearest
# PHF1+ neuron, plus the matching figure: ONE LINE PER DONOR, coloured by Braak.
#
# WHY THIS EXISTS. The calling scripts' raw panels average all 9 donors into a
# single rolling-mean curve, which cannot show whether a gradient is shared across
# the cohort or carried by a subset. This view draws the donors separately, with no
# model imposed.
#
# Sourced by plot_phf1_module_vs_phf1_distance_modelp.R and
# plot_modulescore_vs_phf1_distance_modelp.R. Needs braak_levels and fig_theme
# from R/palettes.R.

# A window holding fewer cells than this is left NA rather than drawn: splitting
# by donor thins the far tail badly, and a "mean" of 3 cells is a spike, not a
# trend. Kept generous because the resulting gaps are informative -- a line that
# stops early tells you that donor has no far-field cells.
ROLLMEAN_MIN_IN_WINDOW  <- 20
ROLLMEAN_MIN_CELLS_DONOR <- 100

# braak_palette in R/palettes.R is built for FILLS: its stage-1 "#FFEDA0" is far
# too pale to read as a thin line on white (the same reason the palettes park the
# pale yellow last for lines). These are darker members of the same warm
# sequential ramp, so the stage ordering and hue family are unchanged, and stage 6
# is the canonical colour exactly. Use ONLY for line geometries; if a third script
# ever needs it, move it into R/palettes.R.
braak_line_palette <- c("1" = "#D4A017", "4" = "#F16913", "6" = "#BD0026")

# Long-format per-donor rolling means for one or more score columns.
# Returns ALL grid rows, including those left NA by min_in_window, so the source
# data records what was too sparse to draw; the plotter filters them.
rollmean_by_donor <- function(dat, score_cols, x_cap, half_window,
                              dist_col = "dist_to_phf1_um",
                              donor_col = "sample_id", braak_col = "Braak",
                              grid_n = 200,
                              min_in_window = ROLLMEAN_MIN_IN_WINDOW,
                              min_cells_donor = ROLLMEAN_MIN_CELLS_DONOR) {
  score_cols <- intersect(score_cols, colnames(dat))
  if (!length(score_cols)) return(NULL)
  grid   <- seq(0, x_cap, length.out = grid_n)
  donors <- sort(unique(as.character(dat[[donor_col]])))
  out <- list()
  for (s in donors) {
    i <- which(as.character(dat[[donor_col]]) == s)
    if (length(i) < min_cells_donor) next
    d  <- as.numeric(dat[[dist_col]])[i]
    bk <- if (!is.null(braak_col) && braak_col %in% colnames(dat))
            as.character(dat[[braak_col]])[i][1] else NA_character_
    for (m in score_cols) {
      y  <- as.numeric(dat[[m]])[i]
      ok <- is.finite(d) & is.finite(y)
      dd <- d[ok]; yy <- y[ok]
      st <- vapply(grid, function(g) {
        k <- which(dd >= g - half_window & dd <= g + half_window)
        n <- length(k)
        if (n < min_in_window) return(c(NA_real_, NA_real_, n))
        c(mean(yy[k]), stats::sd(yy[k]) / sqrt(n), n)   # n >= 20, so sd is defined
      }, numeric(3))
      out[[paste(s, m)]] <- data.frame(
        module = m, sample_id = s, Braak = bk, dist_to_phf1_um = grid,
        roll_mean = st[1, ], sem = st[2, ], n_window = st[3, ],
        n_cells_donor = length(dd), window_um = 2 * half_window,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# One line per donor, coloured by Braak. `facet_col` (optional) facets the panel,
# normally one facet per module/signature with a free y scale because different
# gene sets are not on a common AddModuleScore scale.
plot_rollmean_by_donor <- function(long, xlab, ylab, x_cap,
                                   facet_col = NULL, facet_levels = NULL,
                                   ncol = NULL) {
  d <- long[is.finite(long$roll_mean), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  has_braak <- any(!is.na(d$Braak))
  if (has_braak) d$Braak <- factor(as.character(d$Braak), levels = braak_levels)
  if (!is.null(facet_col) && !is.null(facet_levels))
    d[[facet_col]] <- factor(d[[facet_col]], levels = facet_levels)

  p <- ggplot2::ggplot(d, ggplot2::aes(dist_to_phf1_um, roll_mean, group = sample_id))
  p <- p + if (has_braak)
    ggplot2::geom_line(ggplot2::aes(colour = Braak), linewidth = 0.4, alpha = 0.9)
  else
    ggplot2::geom_line(ggplot2::aes(colour = sample_id), linewidth = 0.4, alpha = 0.9)
  p <- p + if (has_braak)
    ggplot2::scale_colour_manual(values = braak_line_palette, name = "Braak", drop = FALSE)
  else ggplot2::scale_colour_discrete(name = "Donor")
  if (!is.null(facet_col))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)),
                                 scales = "free_y", ncol = ncol)
  p +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::coord_cartesian(xlim = c(0, x_cap)) +
    ggplot2::theme_classic(base_size = 8) + fig_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.width = grid::unit(12, "pt"),
      legend.box.spacing = grid::unit(3, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 6, colour = "black"),
      plot.margin = ggplot2::margin(t = 4, r = 6, b = 4, l = 5))
}

# Shared stats-log body for a per-donor rolling-mean figure.
rollmean_by_donor_stats <- function(long, x_cap, half_window, n_donors_total,
                                    min_in_window = ROLLMEAN_MIN_IN_WINDOW,
                                    min_cells_donor = ROLLMEAN_MIN_CELLS_DONOR) {
  cat("Model-free sliding-window mean of the per-cell AddModuleScore over distance to the\n")
  cat("nearest PHF1+ neuron, computed SEPARATELY WITHIN EACH DONOR. One line per donor,\n")
  cat("coloured by Braak stage. No model, no test: this is a descriptive figure showing the\n")
  cat("between-donor spread.\n\n")
  cat(sprintf("Window width: %g um (half-window %g) | grid: 200 points over 0-%g um\n",
              2 * half_window, half_window, x_cap))
  cat(sprintf("A window with < %d cells is left blank rather than drawn, so a line that stops\n",
              min_in_window))
  cat("  early means that donor has too few cells at that distance -- the gaps are data.\n")
  cat(sprintf("Donors with < %d cells in the celltype are dropped entirely.\n\n", min_cells_donor))
  dn <- unique(long[, c("sample_id", "Braak", "n_cells_donor")])
  dn <- dn[order(dn$Braak, dn$sample_id), , drop = FALSE]
  cat(sprintf("Donors drawn: %d of %d\n", nrow(dn), n_donors_total))
  print(dn, row.names = FALSE)
  drawn <- long[is.finite(long$roll_mean), , drop = FALSE]
  if (nrow(drawn)) {
    reach <- stats::aggregate(dist_to_phf1_um ~ sample_id, data = drawn, FUN = max)
    names(reach)[2] <- "max_distance_drawn_um"
    cat("\nHow far each donor's line reaches (limited by cell density, not by biology):\n")
    print(reach[order(reach$max_distance_drawn_um), ], row.names = FALSE)
  }
  cat("\nNOTE: per-donor curves are noisier than the cohort curve by construction; read the\n")
  cat("  SPREAD and the SIGN agreement between donors, not the wiggle of any single line.\n")
  invisible(NULL)
}
