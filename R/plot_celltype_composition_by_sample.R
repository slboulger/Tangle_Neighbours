#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_celltype_composition_by_sample.R
#
# Figure panels: 1D
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_celltype_composition_by_sample.R
#
# Celltype composition across samples: one 100%-stacked proportional bar per
# sample (sample on x, celltype proportion on y), coloured by celltype_palette.
# Samples ordered by Braak stage (1 -> 4 -> 6) then sample_id.
#
# Triple-output convention, output under plots/celltype_composition_by_sample/.
# Source of truth is seu_PHF1.rds@meta.data. Read-only on the object.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(ggplot2)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # celltype_palette, neuron_order, glia_order, braak_levels, fig_theme

outdir <- "plots/celltype_composition_by_sample"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

WIDTH_CM  <- 5.9    # figure width  (incl. legend)
HEIGHT_CM <- 6.14   # figure height

# ---- load + slim frame ----------------------------------------------------
pick <- function(md, candidates, what, required = TRUE) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) {
    msg <- sprintf("Could not find a column for '%s'. Tried: %s.", what,
                   paste(candidates, collapse = ", "))
    if (required) stop(msg) else { warning(msg); return(NA_character_) }
  }
  hit[1]
}

seu <- readRDS("seu_PHF1.rds")
md  <- seu@meta.data

col_samp <- pick(md, c("sample_id", "sampleID"), "sample_id")
col_ct   <- pick(md, c("celltype"), "celltype")
col_brk  <- pick(md, c("Braak"), "Braak", required = FALSE)

df <- tibble(
  sample_id = as.character(md[[col_samp]]),
  celltype  = as.character(md[[col_ct]])
) %>% filter(!is.na(sample_id), !is.na(celltype))
if (!is.na(col_brk)) df$Braak <- as.character(md[[col_brk]])

# ---- per-sample celltype proportions (base table = robust) -----------------
ctab <- as.data.frame(table(sample_id = df$sample_id, celltype = df$celltype),
                      stringsAsFactors = FALSE)
tot  <- tapply(ctab$Freq, ctab$sample_id, sum)
ctab$n_sample <- tot[ctab$sample_id]
ctab$prop     <- ctab$Freq / ctab$n_sample

# celltype fill order (legend consistency; first level at top of each stack)
ct_levels <- c(neuron_order, glia_order, "Other", "Unassigned", "Unassigned Neuron")
ct_levels <- ct_levels[ct_levels %in% unique(ctab$celltype)]
ct_all    <- c(ct_levels, setdiff(unique(ctab$celltype), ct_levels))
ctab$celltype <- factor(ctab$celltype, levels = ct_all)

# sample order: by Braak (1 -> 4 -> 6) then sample_id
if (!is.na(col_brk)) {
  brk_map <- df %>% distinct(sample_id, Braak)
  brk_map$Braak <- factor(brk_map$Braak, levels = braak_levels)
  brk_map <- brk_map[order(brk_map$Braak, brk_map$sample_id), ]
  samp_order <- brk_map$sample_id
  ctab <- merge(ctab, brk_map, by = "sample_id", all.x = TRUE)
} else {
  samp_order <- sort(unique(ctab$sample_id))
}
ctab$sample_id <- factor(ctab$sample_id, levels = samp_order)

# ---- plot: 100%-stacked proportional bars ---------------------------------
p <- ggplot(ctab, aes(x = sample_id, y = prop, fill = celltype)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = celltype_palette, drop = FALSE, name = "Cell type") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +  # no side padding
  labs(x = "Sample", y = "Proportion of cells") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x  = element_blank(),   # no individual sample names
        axis.ticks.x = element_blank(),
        legend.box.spacing = grid::unit(2, "pt"),  # pull legend to the panel
        legend.margin      = margin(0, 0, 0, 0),
        plot.margin        = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")) +
  guides(fill = guide_legend(ncol = 1, keyheight = grid::unit(7, "pt")))

ggsave(file.path(outdir, "plot_celltype_composition_by_sample.pdf"),
       p, width = WIDTH_CM / 2.54, height = HEIGHT_CM / 2.54, units = "in", device = "pdf")

# ---- source data ----------------------------------------------------------
src <- ctab[, c("sample_id", intersect("Braak", colnames(ctab)),
                "celltype", "Freq", "n_sample", "prop")]
names(src)[names(src) == "Freq"] <- "n"
src <- src[order(src$sample_id, src$celltype), ]
write.table(src, file.path(outdir, "source_data_celltype_composition_by_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- stats log ------------------------------------------------------------
sink(file.path(outdir, "stats_celltype_composition_by_sample.txt"))
cat("Celltype composition across samples - descriptive (no statistical test)\n")
cat("=======================================================================\n\n")
cat(sprintf("Samples: %d  |  celltypes: %d\n", length(samp_order), length(ct_all)))
cat(sprintf("Sample order (x-axis): %s\n\n", paste(samp_order, collapse = ", ")))
cat("Cells per sample:\n"); print(tot)
cat("\nCelltype counts per sample (rows = celltype, cols = sample):\n")
print(table(celltype = factor(df$celltype, levels = ct_all), sample_id = df$sample_id))
cat("\n"); print(sessionInfo())
sink()

cat(sprintf("Done. %d samples, %d celltypes. Outputs in %s\n",
            length(samp_order), length(ct_all), outdir))
