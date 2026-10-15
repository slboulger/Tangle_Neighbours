#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_traj_umap_publication.R
#
# Figure panels: 2B
# Draws the publication UMAP from the per-cell table written by cbln2_slingshot_trajectory.R.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_traj_umap_publication.R
#
# Publication UMAP of a celltype's tangle-associated trajectory: every cell coloured by
# pseudotime, tangle-bearing (PHF1+) cells drawn ON TOP as enlarged circles with a black
# outline and filled with their own pseudotime colour, no frame or axes -- only short UMAP 1 /
# UMAP 2 arrows at the bottom-left corner (Foiani et al. 2026 style).
#
# Reads the trajectory script's per-cell table (plots/<tag>_slingshot_trajectory/
# source_data_<tag>_traj_percell.tsv), so no refit:
#   Rscript R/plot_traj_umap_publication.R                              # CBLN2
#   Rscript R/plot_traj_umap_publication.R Exc-IT-L3-5-CHGA-IL1RAPL2    # L3-5
# OUTPUT (same folder): plot_/source_data_/stats_<tag>_traj_umap_publication.*
#
# Colour = cividis, the pseudotime scale used by every trajectory figure (the PHF1-p95 UMAP
# uses viridis, so the two never get confused). SHOW_CURVES = TRUE adds the trajectory paths.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/cbln2_embedding_utils.R")   # ct_tag(), ct_percell_file()

# ---- CONFIG -----------------------------------------------------------------------------
CT          <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else CBLN2_CT
W_CM        <- 7.0      # canvas width (cm); height follows the UMAP's aspect ratio (no empty bands)
LEGEND_CM   <- 2.0      # width taken by the legend
CELL_SIZE   <- 0.25     # background cells
TANGLE_SIZE <- 1.6      # tangle-bearing circles
TANGLE_STROKE <- 0.3    # their outline
ARROW_FRAC  <- 0.25     # arrow length as a fraction of the smaller UMAP range
SHOW_CURVES <- FALSE    # TRUE to draw slingshot's own lineage curves (busy with several lineages)
SHOW_MEAN_PATH <- FALSE # one smoothed path through pseudotime
PATH_ARROW  <- FALSE    # arrowhead at the tangle-rich end of the mean path
PATH_BINS   <- 40       # equal-count pseudotime bins for the mean path
PATH_DF     <- 6        # smoothing-spline df (lower = smoother)

TAG <- ct_tag(CT)
DIR <- file.path("plots", paste0(TAG, "_slingshot_trajectory"))
NAME <- paste0(TAG, "_traj_umap_publication")
pc <- read.delim(ct_percell_file(CT), colClasses = c(sample_id = "character"))
stopifnot(all(c("UMAP_1", "UMAP_2", "pt", "PHF1") %in% names(pc)))
pc$PHF1 <- as.logical(pc$PHF1)

bg  <- pc %>% arrange(pt)                      # bright (late) cells drawn last
tng <- pc %>% filter(PHF1) %>% arrange(pt)

# Mean trajectory path: median UMAP position of each equal-count pseudotime bin, smoothed with a
# spline in x(pt) and y(pt). With branching lineages it is an AVERAGE route, so it can pass between
# branches; it shows the direction of the axis, not any one lineage.
pb <- pc %>% mutate(b = ntile(pt, PATH_BINS)) %>% group_by(b) %>%
  summarise(pt = median(pt), UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = n(), .groups = "drop")
pt_grid <- seq(min(pb$pt), max(pb$pt), length.out = 200)
mean_path <- tibble(pt = pt_grid,
                    UMAP_1 = predict(smooth.spline(pb$pt, pb$UMAP_1, df = PATH_DF), pt_grid)$y,
                    UMAP_2 = predict(smooth.spline(pb$pt, pb$UMAP_2, df = PATH_DF), pt_grid)$y)

## One builder for both sizes; geometry (corner arrows, limits) and sizes are passed in.
curves_file <- file.path(DIR, paste0("source_data_", TAG, "_traj_curves.tsv"))
build_umap <- function(x0, y0, L, gap, sz, legend = c("right", "top"), xlim = NULL, ylim = NULL) {
  legend <- match.arg(legend)
  cbar <- if (legend == "right")
    guide_colourbar(order = 1, barheight = grid::unit(18, "mm"), barwidth = grid::unit(2.5, "mm")) else
    guide_colourbar(order = 1, direction = "horizontal", title.position = "top",
                    barwidth = grid::unit(sz$barw, "mm"), barheight = grid::unit(sz$barh, "mm"))
  ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2, colour = pt), size = sz$cell, stroke = 0) +
    { if (SHOW_CURVES && file.exists(curves_file))
        geom_path(data = read.delim(curves_file), aes(UMAP_1, UMAP_2, group = lineage), colour = "black", linewidth = 0.4) } +
    geom_point(data = tng, aes(UMAP_1, UMAP_2, fill = pt, shape = "Tangle-bearing"),
               colour = "black", size = sz$tangle, stroke = sz$stroke) +
    { if (SHOW_MEAN_PATH) list(
        geom_path(data = mean_path, aes(UMAP_1, UMAP_2), colour = "white", linewidth = sz$halo, lineend = "round"),
        geom_path(data = mean_path, aes(UMAP_1, UMAP_2), colour = "black", linewidth = sz$path, lineend = "round",
                  arrow = if (PATH_ARROW) grid::arrow(length = grid::unit(sz$path_head, "mm"), type = "closed", angle = 25))) } +
    annotate("segment", x = x0, xend = x0 + L, y = y0, yend = y0, linewidth = sz$axis,
             arrow = grid::arrow(length = grid::unit(sz$axis_head, "mm"), type = "closed")) +
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + L, linewidth = sz$axis,
             arrow = grid::arrow(length = grid::unit(sz$axis_head, "mm"), type = "closed")) +
    # Labels start at the arrows' shared origin and run along them, so they never meet in the corner.
    annotate("text", x = x0, y = y0 - gap, label = "UMAP 1", size = sz$txt, hjust = 0, vjust = 1) +
    annotate("text", x = x0 - gap, y = y0, label = "UMAP 2", size = sz$txt, hjust = 0, vjust = 0, angle = 90) +
    scale_colour_viridis_c(option = "cividis", aesthetics = c("colour", "fill"), name = "Pseudotime", guide = cbar) +
    scale_shape_manual(values = c("Tangle-bearing" = 21), name = NULL,
                       guide = guide_legend(order = 2, override.aes = list(fill = "grey80", colour = "black",
                                                                           size = sz$tangle, stroke = sz$stroke))) +
    coord_equal(xlim = xlim, ylim = ylim, expand = is.null(xlim), clip = "off") +
    theme_void(base_size = 8) +
    theme(legend.title = element_text(size = sz$leg_title), legend.text = element_text(size = sz$leg_text),
          legend.key.size = grid::unit(3, "mm"), legend.position = legend,
          legend.box = if (legend == "top") "horizontal" else "vertical",
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(sz$leg_gap, "pt"),
          plot.margin = sz$margin)
}

## ---- 7 cm version (legend right) ----
dx <- diff(range(pc$UMAP_1)); dy <- diff(range(pc$UMAP_2))
H_CM <- max(4, (W_CM - LEGEND_CM) * (dy * 1.12) / (dx * 1.12) + 0.9)   # height from the aspect ratio
sz_main <- list(cell = CELL_SIZE, tangle = TANGLE_SIZE, stroke = TANGLE_STROKE, path = 0.6, halo = 1.3,
                path_head = 1.6, axis = 0.35, axis_head = 1.2, txt = 2.1, leg_title = 7, leg_text = 6,
                leg_gap = 5.5, margin = margin(2, 2, 12, 12))
## Axis arrows are made LABEL_PAD cm longer than their "UMAP 1"/"UMAP 2" label, and the label sits
## LABEL_GAP cm from the arrow, so text never runs into the arrowhead or the line.
LABEL_PAD <- 0.15
LABEL_GAP <- 0.10
text_cm <- function(label, size_mm) {   # rendered width of a label, in cm
  pdf(NULL); on.exit(invisible(dev.off()))
  grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = size_mm * ggplot2::.pt))),
                     "cm", valueOnly = TRUE)
}
s_main <- (W_CM - LEGEND_CM - 0.8) / (dx * 1.06)   # approx. cm per UMAP unit in the 7 cm figure
L <- max(ARROW_FRAC * min(dx, dy), (text_cm("UMAP 2", sz_main$txt) + LABEL_PAD) / s_main)
p <- build_umap(x0 = min(pc$UMAP_1) - 0.03 * dx, y0 = min(pc$UMAP_2) - 0.03 * dy, L = L,
                gap = LABEL_GAP / s_main, sz = sz_main, legend = "right")
ggsave(file.path(DIR, paste0("plot_", NAME, ".pdf")), p, width = W_CM, height = H_CM, units = "cm", device = "pdf")

## ---- 5.08 x 5.08 cm version: legend ABOVE, tight; everything enlarged ----
## The UMAP keeps its true aspect ratio, so on a square canvas there is spare room. Instead of
## letting it fall between legend and UMAP, the legend's real height is measured, the UMAP is scaled
## to the largest size that fits, the plot sits directly under the legend, and the spare room becomes
## a band at the bottom-left that holds the UMAP 1 / UMAP 2 arrows (never over the cells).
SQ_CM <- 5.08
sz_sq <- list(cell = 0.35, tangle = 1.5, stroke = 0.3, path = 0.75, halo = 1.6, path_head = 1.8,
              axis = 0.4, axis_head = 1.3, txt = 2.5, leg_title = 7.5, leg_text = 6.5, barw = 20, barh = 2.2,
              leg_gap = 2, margin = margin(2, 2, 2, 2))
pt2cm <- 2.54 / 72.27
pdf(NULL)   # text metrics for the legend measurement
g_tmp <- ggplotGrob(build_umap(0, 0, 1, 0.1, sz_sq, legend = "top"))
leg_i <- which(g_tmp$layout$name %in% c("guide-box-top", "guide-box"))[1]
h_leg <- grid::convertHeight(grid::grobHeight(g_tmp$grobs[[leg_i]]), "cm", valueOnly = TRUE)
invisible(dev.off())
W_av <- SQ_CM - 4 * pt2cm                                  # minus left/right margins
H_av <- SQ_CM - 4 * pt2cm - h_leg - sz_sq$leg_gap * pt2cm  # minus margins, legend, legend gap
BAND_CM <- 0.42                                            # left + bottom band: arrow labels only (arrows hug the cells)
s_cm <- min((W_av - BAND_CM) / dx, (H_av - BAND_CM) / dy)  # cm per UMAP unit (largest that fits)
Rx <- W_av / s_cm; Ry <- H_av / s_cm
extra_x <- Rx - (dx + BAND_CM / s_cm)
xlim <- min(pc$UMAP_1) - BAND_CM / s_cm - extra_x / 2 + c(0, Rx)
ylim <- max(pc$UMAP_2) - c(Ry, 0)                          # top flush under the legend; spare room at the bottom
# Axes start just outside the cells' bottom-left corner and run along the cell bounding box, so
# only the labels need the reserved band and the UMAP can be as large as possible.
L_sq_cm <- text_cm("UMAP 2", sz_sq$txt) + LABEL_PAD
p_sq <- build_umap(x0 = min(pc$UMAP_1) - 0.06 / s_cm, y0 = min(pc$UMAP_2) - 0.06 / s_cm, L = L_sq_cm / s_cm,
                   gap = LABEL_GAP / s_cm, sz = sz_sq, legend = "top", xlim = xlim, ylim = ylim)
ggsave(file.path(DIR, paste0("plot_", NAME, "_5.08cm.pdf")), p_sq, width = SQ_CM, height = SQ_CM,
       units = "cm", device = "pdf")

write.table(pc %>% transmute(cell_id, sample_id, UMAP_1, UMAP_2, pt, PHF1,
                             drawn_as = ifelse(PHF1, "tangle-bearing circle", "cell")),
            file.path(DIR, paste0("source_data_", NAME, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
if (SHOW_MEAN_PATH) write.table(mean_path, file.path(DIR, paste0("source_data_", NAME, "_mean_path.tsv")),
                                sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(DIR, paste0("stats_", NAME, ".txt")))
cat(CT, "trajectory UMAP (publication style)\n", strrep("=", 60), "\n", sep = "")
cat("Run:", format(Sys.time()), "\nInput:", ct_percell_file(CT), "\n")
cat(sprintf("Trajectory: %s, k %s, root %s\n", unique(pc$fit_space), unique(pc$k), unique(pc$root)))
cat(sprintf("Cells %d; tangle-bearing %d (circles, fill = own pseudotime); pseudotime %.2f-%.2f\n",
            nrow(pc), nrow(tng), min(pc$pt), max(pc$pt)))
cat(sprintf("Canvas %.1f x %.1f cm; cells size %.2f; circles size %.1f, outline %.2f; curves %s\n",
            W_CM, H_CM, CELL_SIZE, TANGLE_SIZE, TANGLE_STROKE, SHOW_CURVES))
if (SHOW_MEAN_PATH) cat(sprintf("Mean path: median UMAP position of %d equal-count pseudotime bins, smooth.spline df %d in x(pt) and y(pt); runs from early to late pseudotime\n",
                                PATH_BINS, PATH_DF))
cat("Descriptive figure, no test.\n\n"); print(sessionInfo())
sink()
message("Written: ", file.path(DIR, paste0("plot_", NAME, ".pdf")))
