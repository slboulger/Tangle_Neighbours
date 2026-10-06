# ---------------------------------------------------------------------------
# palettes.R
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
# Shared palettes, orderings, and ggplot theme.
# Sourced by every plotting script in the project so colours / styling are
# centrally controlled and reused identically across figures.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# -------------------------------------------------------------------
# Cell-type colour palette
# -------------------------------------------------------------------
# Exc family: sequential warm (light -> dark) following cortical layer depth
# Inh family: sequential cool
# Glia / vascular: distinct hues from a separate qualitative pool
# Unassigned / Other: neutral greys
celltype_palette <- c(
  # Excitatory neurons (renamed in rename_celltypes.R, only present in seu_PHF1.rds$celltype)
  "Exc-IT-L2-3-CBLN2-HOPX"     = "#FDD0A2",
  "Exc-IT-L3-5-CHGA-IL1RAPL2"  = "#FDAE6B",
  "Exc-ET-L5-SPON1-FGD4"       = "#FD8D3C",
  "Exc-IT-L6-CTXN1-ERC2"       = "#E6550D",
  "Exc-CT-L6-SYNPO2-SEMA3E"    = "#A63603",

  # Inhibitory neurons
  "Inh-PVALB"                  = "#9ECAE1",
  "Inh-SST"                    = "#6BAED6",
  "Inh-VIP"                    = "#3182BD",
  "Inh-LAMP5"                  = "#08519C",

  # Glia & vascular
  "Astro"                      = "#74C476",
  "Oligo"                      = "#238B45",
  "OPC"                        = "#A1D99B",
  "Micro"                      = "#9C755F",
  "Endo"                       = "#BC80BD",
  "VLMC"                       = "#B15928",

  # Broad labels (in case used in a broad heatmap / panel)
  "Glutamatergic"              = "#D94801",
  "GABAergic"                  = "#2171B5",

  # Catch-alls
  "Other"                      = "#999999",
  "Unassigned"                 = "#CCCCCC",
  "Unassigned Neuron"          = "#BDBDBD"
)

# Canonical orderings (use to factor() the celltype column before plotting)
neuron_order <- c(
  "Exc-IT-L2-3-CBLN2-HOPX",
  "Exc-IT-L3-5-CHGA-IL1RAPL2",
  "Exc-ET-L5-SPON1-FGD4",
  "Exc-IT-L6-CTXN1-ERC2",
  "Exc-CT-L6-SYNPO2-SEMA3E",
  "Inh-PVALB",
  "Inh-SST",
  "Inh-VIP",
  "Inh-LAMP5"
)

glia_order <- c("Astro", "Oligo", "OPC", "Micro", "Endo", "VLMC")

# -------------------------------------------------------------------
# Broad celltype grouping (collapses the fine `celltype` column)
# -------------------------------------------------------------------
# Eight broad populations for composition-overview figures. Colours reuse the
# matching fine colour when the mapping is 1-to-1, and a representative shade of
# the converging family when several fine types collapse into one.
broad_order <- c(
  "Excitatory Neurons", "Inhibitory Neurons",
  "Astrocytes", "Microglia", "Oligo/OPC", "Vascular",
  "Unassigned Neuron", "Unassigned"
)

# Map fine celltype labels -> broad population (vectorised, base R only).
to_broad_celltype <- function(x) {
  x   <- as.character(x)
  out <- rep(NA_character_, length(x))
  out[grepl("^Exc", x)]         <- "Excitatory Neurons"
  out[grepl("^Inh", x)]         <- "Inhibitory Neurons"
  out[x == "Astro"]             <- "Astrocytes"
  out[x == "Micro"]             <- "Microglia"
  out[x %in% c("Oligo", "OPC")] <- "Oligo/OPC"
  out[x %in% c("Endo", "VLMC")] <- "Vascular"
  out[x == "Unassigned Neuron"] <- "Unassigned Neuron"
  out[x == "Unassigned"]        <- "Unassigned"
  out[is.na(out)]               <- "Other"   # safety catch-all (flag if hit)
  out
}

broad_palette <- c(
  "Excitatory Neurons" = "#D94801",   # Glutamatergic family
  "Inhibitory Neurons" = "#2171B5",   # GABAergic family
  "Astrocytes"         = "#74C476",   # Astro (1-to-1)
  "Microglia"          = "#9C755F",   # Micro (1-to-1)
  "Oligo/OPC"          = "#238B45",   # Oligo (representative)
  "Vascular"           = "#BC80BD",   # Endo (representative)
  "Unassigned Neuron"  = "#BDBDBD",   # 1-to-1
  "Unassigned"         = "#CCCCCC",   # 1-to-1
  "Other"              = "#999999"
)

# -------------------------------------------------------------------
# Braak palette (3 stages present: 1, 4, 6) — sequential YlOrRd-ish
# -------------------------------------------------------------------
braak_palette <- c(
  "1" = "#FFEDA0",
  "4" = "#FEB24C",
  "6" = "#BD0026"
)
braak_levels <- c("1", "4", "6")

# Darkened Braak palette for LINE and POINT geometries. The canonical stage-1
# fill #FFEDA0 is a pale yellow that is effectively invisible as a thin line or a
# small point on white. These are darker members of the same warm sequential ramp,
# so the stage ordering and hue family are unchanged and stage 6 is the canonical
# colour exactly. Use braak_palette for fills (boxes, bars, ribbons) and this for
# lines/points.
#
# R/rollmean_by_donor.R carries an identical local definition (the same three values).
braak_line_palette <- c("1" = "#D4A017", "4" = "#F16913", "6" = "#BD0026")

# -------------------------------------------------------------------
# Nature-compact ggplot theme
# Sized for ~half-page subfigures; readable when shrunk.
# -------------------------------------------------------------------
fig_theme <- theme(
  axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7,
                              colour = "black"),
  axis.text.y  = element_text(size = 7, colour = "black"),
  axis.title   = element_text(size = 8),
  axis.ticks   = element_line(linewidth = 0.3),
  axis.line    = element_line(linewidth = 0.3),
  legend.text  = element_text(size = 6),
  legend.title = element_text(size = 7),
  legend.key.height = grid::unit(8, "pt"),
  legend.key.width  = grid::unit(6, "pt"),
  # Extra left + bottom padding so 45deg long celltype labels
  # (e.g. "Exc-IT-L3-5-CHGA-IL1RAPL2") aren't clipped at the PDF edge.
  plot.margin  = margin(t = 4, r = 6, b = 10, l = 38, unit = "pt")
)

# Clean, high-impact forest-plot theme (used by the LR forest scripts). No panel
# border box, no vertical gridlines, a single thin black x-axis; faint horizontal
# guides to read estimates across rows; y axis line/ticks dropped (categorical);
# blank facet-strip backgrounds; compact bottom legend. Deliberately NOT theme_bw
# (the grey box/gridlines read as dated).
forest_theme <- theme_classic(base_size = 8) +
  theme(
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.x        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.x       = element_line(linewidth = 0.3, colour = "black"),
    axis.text          = element_text(colour = "black", size = 7),
    axis.text.x        = element_text(size = 7, colour = "black"),
    axis.title.x       = element_text(size = 8, colour = "black"),
    axis.title.y       = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 6, colour = "black"),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 6.5, colour = "black"),
    legend.key.size    = grid::unit(7, "pt"),
    legend.box.spacing = grid::unit(2, "pt"),
    legend.margin      = margin(0, 0, 0, 0),
    plot.title         = element_text(size = 8, face = "bold", hjust = 0, colour = "black"),
    plot.margin        = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt")
  )

##  ............................................................................
##  Legend completeness                                                     ####

#' Carry one UNDRAWN row per missing factor level, so a guide key gets its glyph.
#'
#' A ggplot guide key is built from LAYER DATA, so a level that never occurs -- every
#' series n.s., or every series significant -- yields a key with its label and an EMPTY
#' glyph. `drop = FALSE` and `override.aes` do NOT fix this: the scale keeps the key,
#' but there is no data row to draw the line or point from. Adding one row whose plotted
#' values are all NA gives the key something to draw from while putting nothing on the
#' panel. Verified: an all-n.s. linetype guide goes from 1 drawn key to 2,
#' matching the both-levels-present case; the same holds for shape guides.
#'
#' NOT needed for COLOUR/FILL guides -- a colour key can be drawn from the scale value
#' alone, so `drop = FALSE` is sufficient there. Only linetype, linewidth and shape keys,
#' which additionally need the colour/size that only layer data carries, require this.
#'
#' Generalises the local copies in R/phf1_channel_intensity_exc.R,
#' R/imc_phf1_glia_distance.R, R/imc_phf1_dna_dist_phf1_panel.R and
#' R/imc_phf1_morphology_dist_phf1_panel.R (which each hard-code their own level set and
#' column list); a local definition shadows this one.
#'
#' @param d Data frame; NULL or zero-row input is returned unchanged.
#' @param col Column holding the factor whose levels drive the guide.
#' @param levels_want Full set of levels the guide should show.
#' @param keep Numeric columns to leave alone. Default none: every numeric in the filler
#'   row is NA'd, which is what keeps it off the panel.
#' @return `d` with at most `length(levels_want)` extra all-NA rows appended.
pad_levels <- function(d, col, levels_want, keep = character(0)) {
  if (is.null(d) || !nrow(d)) return(d)
  d[[col]] <- factor(as.character(d[[col]]), levels = levels_want)
  miss <- setdiff(levels_want, as.character(d[[col]]))
  if (!length(miss)) return(d)
  filler <- d[rep(1L, length(miss)), , drop = FALSE]
  filler[[col]] <- factor(miss, levels = levels_want)
  for (v in setdiff(names(filler)[vapply(filler, is.numeric, logical(1))], keep))
    filler[[v]] <- NA_real_
  # List columns would be CLONED from the donor row, so a caller that carries one
  # (e.g. per-donor points) must extract it before padding -- see the forest plot in
  # R/plot_reactome_stress_death_vs_phf1_distance.R.
  dplyr::bind_rows(d, filler)
}
