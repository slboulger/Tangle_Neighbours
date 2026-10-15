#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_cbln2_umap_phf1_p95.R
#
# Figure panels: 2B - upstream step
# Builds and caches the Exc-IT-L2-3-CBLN2-HOPX SCT/PCA/Harmony/UMAP embedding that the
# trajectory uses (via cbln2_embedding_utils.R); its own UMAP plot is not a paper panel.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_cbln2_umap_phf1_p95.R
#
# Exc-IT-L2-3-CBLN2-HOPX neurons on a UMAP, each cell coloured by its normalised
# 95th-percentile PHF1 mask intensity (phf1_intensity_p95_z), with the manually
# annotated tangle-bearing (PHF1+) neurons ringed on top. A visualisation of where
# the immunofluorescence signal and the manual calls sit on the expression manifold.
#
# INPUT
#   celltype_sce/Exc-IT-L2-3-CBLN2-HOPX_sce.qs -- 26,585 cells, 188 PHF1+. Not
#   celltype_sce_neighbours/*.qs, which excludes the PHF1+ cells, and not
#   seu_PHF1.rds: this SCE already carries counts, PHF1 and phf1_intensity_p95_z.
#
# COLOUR
#   phf1_intensity_p95_z is the Methods definition (add_phf1_intensity.R): robust z
#   (median/MAD) within sample, over ALL cells of the sample. It is NOT recomputed
#   within CBLN2, so the scale matches Fig. 1F and plot_phf1_intensity_validation*.R.
#   The distribution has a heavy right tail, so the colour scale is squished to the
#   [--cap_lo_q, --cap_hi_q] quantiles over the plotted cells; otherwise a handful of
#   tangles would use the whole colour range. The top colourbar label reads ">= cap".
#   The cap values and how many cells (and PHF1+ cells) sit beyond them are in the
#   stats log. Cells with NA intensity (mask < 50% inside the PHF1 image) are drawn
#   grey, underneath. Cells are drawn in ascending intensity so the bright ones are
#   not buried.
#
# EMBEDDING
#   CBLN2 is re-embedded on its own (gene filter, LogNormalize + vst HVGs, SCTransform
#   glmGamPoi regressing nCount_RNA + percent.neg + percent.false, PCA, harmony on
#   sample_id, UMAP n.epochs 200, seed 42). PCs by the project's min(co1, co2) elbow
#   rule unless --dims is given. Coordinates are cached in subcluster_phf1/. The recipe
#   and cache live in R/cbln2_embedding_utils.R, shared with the trajectory scripts so
#   the trajectory sits on exactly this UMAP.
#   The SCE's stored reducedDim("UMAP") is not used: it is the whole-object embedding,
#   in which PHF1+ cells group partly by library size; the re-embedding regresses
#   nCount_RNA out.
#
# STATS
#   Descriptive (how well the colour separates the ringed cells). They mirror
#   plot_phf1_intensity_validation_neurons.R (rank AUC, paired donor-mean difference +
#   t CI, dz, paired Wilcoxon), restricted to CBLN2. Donor is the inferential unit; a
#   row restricted to donors with >= 10 PHF1+ cells is reported as well.
#
# OUTPUT (the standard triple), plots/cbln2_umap_phf1_p95/, <name> = cbln2_umap_phf1_p95:
#   plot_<name>.pdf                    -- the UMAP
#   plot_<name>_interactive.html       -- the same UMAP in plotly (WebGL): zoom/pan, hover gives
#                                         cell_id, donor, Braak, PHF1 call, p95 (raw + z), nCount;
#                                         click a legend entry to hide/show it. Not self-contained
#                                         unless pandoc is on the PATH: keep the lib/ folder beside it.
#                                         Skipped with a warning if plotly is not installed (this
#                                         script does not install it).
#                                         When run inside RStudio the widget also opens in the Viewer.
#   source_data_<name>.tsv             -- one row per drawn cell (coords, raw + plotted colour, PHF1)
#   stats_<name>.txt                   -- log: inputs, embedding provenance, colour caps, effect sizes, sessionInfo
#   stats_<name>_auc.tsv               -- per-donor AUC / medians
#   stats_<name>_effectsize.tsv        -- headline numbers + donor-paired effect size
#
# RUN (an R with Seurat, glmGamPoi and harmony):
#   Rscript R/plot_cbln2_umap_phf1_p95.R                # builds the embedding cache first time
#   Rscript R/plot_cbln2_umap_phf1_p95.R --dims 1:15    # override the elbow rule
# Nothing is written back to any canonical object.

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(SummarizedExperiment)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(argparse)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # fig_theme, braak_levels
source("R/cbln2_embedding_utils.R")  # the shared CBLN2 embedding (also used by the trajectory scripts)

parser <- ArgumentParser(description = "Celltype UMAP coloured by PHF1 p95 z, PHF1+ ringed")
parser$add_argument("--celltype", default = CBLN2_CT,
                    help = "celltype to embed (default CBLN2; outputs/caches named by ct_tag(), CBLN2 paths unchanged)")
parser$add_argument("--dims", default = "elbow",
                    help = "PCs for harmony/UMAP: 'elbow' (project rule) or N or A:B")
parser$add_argument("--norm", default = "sct", choices = c("sct", "lognorm"),
                    help = "normalisation for the embedding; sct is the project recipe")
parser$add_argument("--cap_lo_q", type = "double", default = 0.01,
                    help = "lower colour-scale quantile (squished below)")
parser$add_argument("--cap_hi_q", type = "double", default = 0.99,
                    help = "upper colour-scale quantile (squished above)")
parser$add_argument("--force_rebuild", action = "store_true", default = FALSE,
                    help = "ignore the cached embedding")
args <- parser$parse_args()
stopifnot(args$cap_lo_q >= 0, args$cap_hi_q <= 1, args$cap_lo_q < args$cap_hi_q)

CT          <- args$celltype
SCE_PATH    <- ct_sce_path(CT)
COLOUR_COL  <- "phf1_intensity_p95_z"
OUTPUT_DIR  <- file.path("plots", paste0(ct_tag(CT), "_umap_phf1_p95"))
NAME        <- paste0(ct_tag(CT), "_umap_phf1_p95")
SEED        <- EMB_SEED

MIN_POS_FOREST <- 10           # donors with >= 10 PHF1+ cells, reported as a second row
PHF1_COLOUR    <- "#D55E00"    # vermillion, the project's PHF1 colour
NA_COLOUR      <- "#D9D9D9"    # grey85 as hex: plotly does not read R colour names

set.seed(SEED)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

##  ............................................................................
##  Helpers                                                                 ####
##  ............................................................................

# rank-based ROC-AUC (no pROC dependency) - identical to plot_phf1_intensity_validation*.R
auc_fn <- function(score, pos) {
  ok <- !is.na(score)
  score <- score[ok]; pos <- as.logical(pos[ok])
  n1 <- sum(pos); n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

write_tsv_ <- function(x, file) {
  write.table(x, file.path(OUTPUT_DIR, file), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

##  ............................................................................
##  Load                                                                    ####
##  ............................................................................
message("Reading ", SCE_PATH)
sce <- qs::qread(SCE_PATH)
cd  <- as.data.frame(colData(sce))
cd$cell_id <- colnames(sce)
cd$PHF1    <- as.logical(as.character(cd$PHF1))
stopifnot(all(as.character(cd$celltype) == CT), !anyNA(cd$PHF1), COLOUR_COL %in% colnames(cd))

##  ............................................................................
##  Embedding                                                               ####
##  ............................................................................
cached <- load_or_build_cbln2_embedding(
  get_inputs = function() {
    meta <- cd[, setdiff(colnames(cd), "cell_id"), drop = FALSE]
    meta$sample_id <- factor(as.character(meta$sample_id))
    list(counts = assay(sce, "counts"), meta = meta)
  },
  norm = args$norm, dims = args$dims, need_full = FALSE, force_rebuild = args$force_rebuild, ct = CT)
coords <- cached$coords
prep   <- cached$prep
rm(sce); invisible(gc())

##  ............................................................................
##  Plot frame                                                              ####
##  ............................................................................
df <- cd %>%
  transmute(cell_id, sample_id = as.character(sample_id),
            Braak = factor(as.character(Braak), levels = braak_levels),
            PHF1, nCount_RNA, phf1_intensity_p95, phf1_intensity_p95_z = .data[[COLOUR_COL]]) %>%
  inner_join(coords, by = "cell_id")
stopifnot(nrow(df) == nrow(cd))   # every cell has coordinates

z_ok   <- df$phf1_intensity_p95_z[!is.na(df$phf1_intensity_p95_z)]
cap_lo <- unname(quantile(z_ok, args$cap_lo_q))
cap_hi <- unname(quantile(z_ok, args$cap_hi_q))
df$colour_value <- pmin(pmax(df$phf1_intensity_p95_z, cap_lo), cap_hi)
df$drawn_as     <- ifelse(is.na(df$phf1_intensity_p95_z), "NA_grey",
                          ifelse(df$PHF1, "coloured+ringed", "coloured"))
df <- df %>% arrange(!is.na(phf1_intensity_p95_z), phf1_intensity_p95_z)   # NA first, then ascending

# Colourbar: pretty breaks inside the caps, plus the caps themselves marked as squished.
span <- cap_hi - cap_lo
br   <- pretty(c(cap_lo, cap_hi), n = 4)
br   <- br[br > cap_lo + 0.12 * span & br < cap_hi - 0.12 * span]
brks <- c(cap_lo, br, cap_hi)
labs_ <- c(sprintf("<=%.1f", cap_lo), format(br, trim = TRUE), sprintf(">=%.1f", cap_hi))

pos <- df[df$PHF1, ]
p <- ggplot() +
  geom_point(data = df[is.na(df$phf1_intensity_p95_z), ], aes(UMAP_1, UMAP_2),
             colour = NA_COLOUR, size = 0.25, stroke = 0) +
  geom_point(data = df[!is.na(df$phf1_intensity_p95_z), ],
             aes(UMAP_1, UMAP_2, colour = colour_value), size = 0.25, stroke = 0) +
  # PHF1+ cells redrawn larger in their own colour, so the ring does not hide it.
  geom_point(data = pos[!is.na(pos$phf1_intensity_p95_z), ],
             aes(UMAP_1, UMAP_2, colour = colour_value), size = 0.8, stroke = 0) +
  geom_point(data = pos, aes(UMAP_1, UMAP_2, shape = "Tangle-bearing"),
             colour = PHF1_COLOUR, fill = NA, size = 1.3, stroke = 0.35) +
  scale_colour_viridis_c(limits = c(cap_lo, cap_hi), breaks = brks, labels = labs_,
                         name = "PHF1 p95\n(per-sample z)") +
  scale_shape_manual(values = c("Tangle-bearing" = 21), name = NULL) +
  guides(colour = guide_colourbar(order = 1),
         shape  = guide_legend(order = 2, override.aes = list(colour = PHF1_COLOUR,
                                                             size = 1.6, stroke = 0.4))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.key.height = grid::unit(14, "pt"),
        plot.margin = margin(t = 4, r = 4, b = 4, l = 4, unit = "pt"))
ggsave(file.path(OUTPUT_DIR, paste0("plot_", NAME, ".pdf")), p,
       width = 8, height = 6.2, units = "cm", device = "pdf")

write_tsv_(df %>% select(cell_id, sample_id, Braak, PHF1, UMAP_1, UMAP_2,
                         phf1_intensity_p95, phf1_intensity_p95_z, colour_value, drawn_as),
           paste0("source_data_", NAME, ".tsv"))
if (interactive()) print(p)

##  ............................................................................
##  Interactive UMAP (plotly, WebGL)                                        ####
##  ............................................................................
## Built with plot_ly directly, not ggplotly(): ggplotly on ~26k points is slow and
## drops the squished colourbar labels. Same encoding as the PDF: viridis on the
## capped value, NA grey underneath, ascending draw order, PHF1+ ringed vermillion.
## Hover reports the UNCAPPED z, so a saturated cell still shows its true value.
html_file <- file.path(OUTPUT_DIR, paste0("plot_", NAME, "_interactive.html"))
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  hover <- with(df, sprintf(
    "%s<br>donor %s (Braak %s)<br>PHF1 call: %s<br>p95 z %.2f (raw %.0f)<br>nCount %d",
    cell_id, sample_id, Braak, ifelse(PHF1, "tangle-bearing", "-"),
    phf1_intensity_p95_z, phf1_intensity_p95, as.integer(nCount_RNA)))
  df$hover <- hover
  na_df  <- df[is.na(df$phf1_intensity_p95_z), ]
  col_df <- df[!is.na(df$phf1_intensity_p95_z), ]
  pos_df <- df[df$PHF1, ]
  viridis_scale <- lapply(seq(0, 1, length.out = 11), function(v)
    list(v, viridisLite::viridis(11)[round(v * 10) + 1]))

  pi <- plotly::plot_ly(type = "scattergl", mode = "markers") %>%
    plotly::add_trace(data = na_df, x = ~UMAP_1, y = ~UMAP_2, text = ~hover,
                      hoverinfo = "text", name = "NA intensity",
                      marker = list(color = NA_COLOUR, size = 3)) %>%
    plotly::add_trace(data = col_df, x = ~UMAP_1, y = ~UMAP_2, text = ~hover,
                      hoverinfo = "text", name = "PHF1 p95 z", showlegend = FALSE,
                      marker = list(color = ~colour_value, colorscale = viridis_scale,
                                    cmin = cap_lo, cmax = cap_hi, size = 3,
                                    colorbar = list(title = "PHF1 p95<br>(per-sample z)",
                                                    tickvals = brks, ticktext = labs_))) %>%
    plotly::add_trace(data = pos_df, x = ~UMAP_1, y = ~UMAP_2, text = ~hover,
                      hoverinfo = "text", name = "Tangle-bearing",
                      marker = list(symbol = "circle-open", color = PHF1_COLOUR, size = 9,
                                    line = list(width = 1.5, color = PHF1_COLOUR))) %>%
    plotly::layout(xaxis = list(title = "UMAP 1", zeroline = FALSE, showgrid = FALSE),
                   yaxis = list(title = "UMAP 2", zeroline = FALSE, showgrid = FALSE,
                                scaleanchor = "x"),
                   legend = list(orientation = "h", x = 0, y = -0.12), plot_bgcolor = "white")
  pi <- plotly::config(pi, toImageButtonOptions = list(format = "svg", filename = NAME))
  # Self-contained needs pandoc; otherwise the JS libraries go in lib/ beside the HTML.
  htmlwidgets::saveWidget(pi, normalizePath(html_file, mustWork = FALSE),
                          selfcontained = requireNamespace("rmarkdown", quietly = TRUE) &&
                                          isTRUE(rmarkdown::pandoc_available()),
                          libdir = "lib", title = paste(CT, "UMAP - PHF1 p95"))
  message("Interactive UMAP: ", html_file)
  if (interactive()) print(pi)
} else {
  warning("plotly/htmlwidgets not installed in this R -- interactive HTML skipped. ",
          "The PDF and source data are complete.", call. = FALSE)
  html_file <- NA_character_
}

##  ............................................................................
##  Stats: how well the colour separates the ringed cells                   ####
##  ............................................................................
dz_df <- df[!is.na(df$phf1_intensity_p95_z), ]
overall_auc <- auc_fn(dz_df$phf1_intensity_p95_z, dz_df$PHF1)

per_donor <- dz_df %>%
  group_by(sample_id, Braak) %>%
  summarise(n = n(), n_phf1_pos = sum(PHF1),
            median_z_pos = if (any(PHF1)) median(phf1_intensity_p95_z[PHF1]) else NA_real_,
            median_z_neg = median(phf1_intensity_p95_z[!PHF1]),
            mean_z_pos   = if (any(PHF1)) mean(phf1_intensity_p95_z[PHF1]) else NA_real_,
            mean_z_neg   = mean(phf1_intensity_p95_z[!PHF1]),
            auc = auc_fn(phf1_intensity_p95_z, PHF1), .groups = "drop") %>%
  arrange(desc(n_phf1_pos))
write_tsv_(per_donor, paste0("stats_", NAME, "_auc.tsv"))

paired_row <- function(min_pos, label) {
  d  <- per_donor %>% filter(n_phf1_pos >= min_pos)
  dd <- d$mean_z_pos - d$mean_z_neg
  if (length(dd) < 3) return(tibble(donor_set = label, n_donors = length(dd)))
  tt <- t.test(dd)
  ww <- suppressWarnings(wilcox.test(d$mean_z_pos, d$mean_z_neg, paired = TRUE))
  tibble(donor_set = label, n_donors = length(dd), n_phf1_pos = sum(d$n_phf1_pos),
         paired_diff_mean = mean(dd), paired_diff_ci_lo = tt$conf.int[1],
         paired_diff_ci_hi = tt$conf.int[2], dz = mean(dd) / sd(dd),
         n_donors_same_direction = sum(sign(dd) == sign(mean(dd))),
         auc_donor_min = min(d$auc), auc_donor_max = max(d$auc),
         wilcox_paired_V = unname(ww$statistic), wilcox_paired_p = ww$p.value)
}
eff <- bind_rows(paired_row(1, "donors with >=1 PHF1+"),
                 paired_row(MIN_POS_FOREST, sprintf("donors with >=%d PHF1+", MIN_POS_FOREST))) %>%
  mutate(n_cells = nrow(dz_df),
         n_tangle_bearing = sum(dz_df$PHF1), auc_pooled = overall_auc,
         colour_cap_lo = cap_lo, colour_cap_hi = cap_hi, .before = 1)
write_tsv_(eff, paste0("stats_", NAME, "_effectsize.tsv"))

##  ............................................................................
##  Stats log                                                               ####
##  ............................................................................
old_width <- options(width = 160)
sink(file.path(OUTPUT_DIR, paste0("stats_", NAME, ".txt")))
cat(CT, "on UMAP, coloured by", COLOUR_COL, "- annotated PHF1+ ringed\n")
cat(strrep("=", 78), "\n\n")
cat("Script   : R/plot_cbln2_umap_phf1_p95.R\n")
cat("Run      :", format(Sys.time()), "\n")
cat("Input    :", SCE_PATH, "\n")
cat("HTML     :", html_file, "\n")
cat("Args     :", paste(names(args), unlist(lapply(args, paste, collapse = ",")), sep = "=", collapse = "  "), "\n\n")

cat("--- Cells ---\n")
cat(sprintf("Drawn: %d cells (%d PHF1+, %d PHF1-); %d with NA intensity drawn grey (%d of them PHF1+)\n\n",
            nrow(df), sum(df$PHF1), sum(!df$PHF1), sum(is.na(df$phf1_intensity_p95_z)),
            sum(is.na(df$phf1_intensity_p95_z) & df$PHF1)))

cat("--- Embedding ---\n")
str(prep, give.attr = FALSE, no.list = TRUE)
cat("\n")

cat("--- Colour scale ---\n")
cat(sprintf("%s squished to [%.3f, %.3f] = its q%g / q%g quantiles over the drawn cells.\n",
            COLOUR_COL, cap_lo, cap_hi, args$cap_lo_q, args$cap_hi_q))
cat(sprintf("Cells above the upper cap: %d (%d of %d PHF1+ = %.0f%%; %d PHF1-)\n",
            sum(dz_df$phf1_intensity_p95_z > cap_hi), sum(dz_df$phf1_intensity_p95_z > cap_hi & dz_df$PHF1),
            sum(dz_df$PHF1), 100 * mean(dz_df$phf1_intensity_p95_z[dz_df$PHF1] > cap_hi),
            sum(dz_df$phf1_intensity_p95_z > cap_hi & !dz_df$PHF1)))
cat(sprintf("Cells below the lower cap: %d\n", sum(dz_df$phf1_intensity_p95_z < cap_lo)))
cat("Quantiles of", COLOUR_COL, "by manual call:\n")
print(do.call(rbind, tapply(dz_df$phf1_intensity_p95_z, ifelse(dz_df$PHF1, "PHF1+", "PHF1-"),
                            quantile, probs = c(0, .05, .25, .5, .75, .95, 1))))
cat(sprintf("\nPHF1- cells brighter than the median PHF1+ cell (%.2f): %d of %d (%.2f%%)\n\n",
            median(dz_df$phf1_intensity_p95_z[dz_df$PHF1]),
            sum(!dz_df$PHF1 & dz_df$phf1_intensity_p95_z > median(dz_df$phf1_intensity_p95_z[dz_df$PHF1])),
            sum(!dz_df$PHF1),
            100 * mean(dz_df$phf1_intensity_p95_z[!dz_df$PHF1] > median(dz_df$phf1_intensity_p95_z[dz_df$PHF1]))))

cat("--- Effect size: does the colour separate the ringed cells? ---\n")
cat("Descriptive. Mirrors plot_phf1_intensity_validation_neurons.R.\n")
cat("Donor is the inferential unit; the pooled AUC is cell-level.\n\n")
cat(sprintf("Pooled rank ROC-AUC (p95 z vs manual call): %.4f\n\n", overall_auc))
cat("Per donor:\n")
print(as.data.frame(per_donor), digits = 3, row.names = FALSE)
cat("\nDonor-paired difference in mean p95 z (PHF1+ - PHF1-), t 95% CI, dz, paired Wilcoxon:\n")
print(as.data.frame(eff %>% select(donor_set, n_donors, any_of(c(
  "n_phf1_pos", "paired_diff_mean", "paired_diff_ci_lo", "paired_diff_ci_hi", "dz",
  "n_donors_same_direction", "wilcox_paired_V", "wilcox_paired_p")))), digits = 4, row.names = FALSE)
cat("\n")

cat("--- sessionInfo ---\n")
print(sessionInfo())
sink()
options(old_width)

message(sprintf("%d cells (%d PHF1+) | colour cap [%.2f, %.2f] | pooled AUC %.3f -> %s",
                nrow(df), sum(df$PHF1), cap_lo, cap_hi, overall_auc, OUTPUT_DIR))
