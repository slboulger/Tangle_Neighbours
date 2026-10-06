# ---------------------------------------------------------------------------
# plot_rrho2_du_annotated.R
#
# Figure panels: 4E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# =============================================================================
# plot_rrho2_du_annotated.R
#
# Scatter of the du discordant genes (DOWN in PHF1+ / UP near tangles), positioned by
# their two moderated t-statistics (x = cell-autonomous t, - = down in PHF1+;
# y = proximity t, + = up near tangles) and coloured by the MANUAL functional module
# annotation in rrho2_du_annotated.csv. Key genes are labelled.
#
# Sibling of R/plot_rrho2_ud_top100_annotated.R (same grammar, opposite discordant
# direction). Shared modules (proteostasis, metabolic, adhesion) keep the same colours
# across the two panels. Reads only existing files -- no RRHO2 / network needed.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel)
})
proj_hpc   <- "<PROJECT_ROOT>/phf1_v2"
proj_local <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(proj_hpc)) proj_hpc else proj_local)
source("R/palettes.R")   # fig_theme

# ------------------------------- config --------------------------------------
ct        <- "Exc-IT-L2-3-CBLN2-HOPX"
rrho2_dir <- file.path("plots/decompose_tangle_response_phf1", ct, "rrho2")
annot_f   <- file.path(rrho2_dir, "rrho2_du_annotated.csv")
tval_f    <- file.path(rrho2_dir, paste0("source_data_rrho2_tstat_quadrant_", ct, ".tsv"))
out_dir   <- rrho2_dir
slug      <- paste0("rrho2_du_annotated_", gsub("[^A-Za-z0-9]+", "_", ct))

LABEL_GENES <- c("CTSD", "PLD3", "GAP43", "STMN2", "NEFL", "UBB",
                 "HSP90AB1", "GRIN1", "ATP1A3", "DNM1")

# Module palette (Paul Tol CB-safe 'muted'); shared modules keep the SAME hue as the ud panel
# (proteostasis rose, metabolic sand, adhesion grey) so colours read consistently across figures.
module_levels <- c("lysosomal", "synaptic_vesicle", "receptor_channel", "adhesion",
                   "axonal_growth", "proteostasis", "redox", "metabolic", "signalling", "transcription")
module_pal <- c(
  lysosomal        = "#332288",   # indigo
  synaptic_vesicle = "#44AA99",   # teal
  receptor_channel = "#88CCEE",   # cyan
  adhesion         = "grey55",    # (shared with ud panel)
  axonal_growth    = "#117733",   # green
  proteostasis     = "#CC6677",   # rose  (shared with ud panel)
  redox            = "#999933",   # olive
  metabolic        = "#DDCC77",   # sand  (shared with ud panel)
  signalling       = "#AA4499",   # purple
  transcription    = "#882255")   # wine
module_labels <- c(
  lysosomal = "Lysosomal", synaptic_vesicle = "Synaptic vesicle",
  receptor_channel = "Receptor channel", adhesion = "Adhesion",
  axonal_growth = "Axonal growth", proteostasis = "Proteostasis",
  redox = "Redox", metabolic = "Metabolic", signalling = "Signalling",
  transcription = "Transcription")

# ------------------------------- load ----------------------------------------
stopifnot(file.exists(annot_f), file.exists(tval_f))
ann <- read.csv(annot_f, stringsAsFactors = FALSE)
tv  <- read.delim(tval_f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

d <- ann %>% left_join(tv[, c("gene", "ca_t", "prox_t")], by = "gene")
miss_t <- d$gene[!is.finite(d$ca_t) | !is.finite(d$prox_t)]
if (length(miss_t)) message("WARNING: no t-values for ", length(miss_t), " gene(s): ",
                            paste(miss_t, collapse = ", "))
d <- d %>% filter(is.finite(ca_t), is.finite(prox_t))
unknown_mod <- setdiff(unique(d$module), module_levels)
if (length(unknown_mod)) stop("modules not in the du palette: ", paste(unknown_mod, collapse = ", "),
     "
  These would render as grey points with no legend entry. Either add them to",
     "
  module_levels/module_pal/module_labels, or fold them into an explicit",
     "
  category as plot_rrho2_ud_top100_annotated.R does for CONTAMINATION_RISK",
     "
  and IMPLAUSIBLE.")
d$module <- factor(d$module, levels = module_levels)

lab <- d %>% filter(gene %in% LABEL_GENES)
missing_lab <- setdiff(LABEL_GENES, d$gene)
if (length(missing_lab)) message("NOTE: label genes absent from the du set (cannot be plotted): ",
                                 paste(missing_lab, collapse = ", "))

# ------------------------------- plot ----------------------------------------
# common axis spans (from BOTH the ud and du gene sets) so the ud/du panels are a matched pair --
# same scale and size, each centred on its own data. Identical block to the ud script.
ud_g <- readLines(file.path(rrho2_dir, paste0("rrho2_genes_ud_top100_", ct, ".txt")))
du_g <- read.csv(file.path(rrho2_dir, "rrho2_du_annotated.csv"))$gene
.span <- function(g, col) diff(range(tv[[col]][match(g, tv$gene)], na.rm = TRUE))
S_x <- max(.span(ud_g, "ca_t"),   .span(du_g, "ca_t"))
S_y <- max(.span(ud_g, "prox_t"), .span(du_g, "prox_t"))
xr <- mean(range(d$ca_t))   + c(-1, 1) * S_x / 2
yr <- mean(range(d$prox_t)) + c(-1, 1) * S_y / 2
p <- ggplot(d, aes(ca_t, prox_t, colour = module)) +
  geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
  ggrepel::geom_text_repel(data = lab, aes(ca_t, prox_t, label = gene), inherit.aes = FALSE,
             size = 2, colour = "grey15", max.overlaps = Inf,
             min.segment.length = 0, segment.size = 0.2, box.padding = 0.3) +
  scale_colour_manual(values = module_pal, labels = module_labels, breaks = module_levels,
                      drop = FALSE, name = "Module") +
  labs(x = expression("Cell-autonomous " * italic(t) * "  (- = down in PHF1+ vs PHF1-)"),
       y = expression("Proximity " * italic(t) * "  (+ = up near tangles)")) +
  coord_fixed(ratio = 1, xlim = xr, ylim = yr) +   # equal scale + same spans as ud -> matched pair
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_classic(base_size = 8) + fig_theme +   # house standard (axis.title 8, legend 6/7 via fig_theme)
  theme(legend.position = "right", legend.key.size = grid::unit(7, "pt"),
        plot.margin = margin(2, 2, 2, 2))

# identical canvas to the ud panel so the pair renders at the same size and scale
ggsave(file.path(out_dir, paste0("plot_", slug, ".pdf")),
       p, width = 9.18, height = 6.0, units = "cm", device = "pdf")

# ------------------------------- source data + stats -------------------------
keep <- intersect(c("gene", "module", "ca_t", "prox_t", "annotation", "baseline_abundance", "flag"),
                  names(d))
sd <- d[order(d$module, d$ca_t), keep]
write.table(sd, file.path(out_dir, paste0("source_data_", slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug, ".txt")))
cat("du discordant genes (down in PHF1+ / up near tangles) --", ct, "\n")
cat("x = cell-autonomous t (limma; - down in PHF1+); y = proximity t (dream, negated; + up near tangles)\n")
cat("Colour = manual functional module (rrho2_du_annotated.csv).\n\n")
cat("n plotted:", nrow(d), " (of", nrow(ann), "annotated)\n")
if (length(miss_t)) cat("dropped (no t-value):", paste(miss_t, collapse = ", "), "\n")
cat("\nGenes per module:\n"); print(as.data.frame(table(module = d$module)))
if ("flag" %in% names(d)) { cat("\nFlag tally:\n"); print(as.data.frame(table(flag = d$flag))) }
cat("\nSpearman r(ca_t, prox_t) over the du set:",
    round(suppressWarnings(cor(d$ca_t, d$prox_t, method = "spearman")), 3), "\n")
cat("\nLabelled genes present:", paste(sort(lab$gene), collapse = ", "), "\n")
if (length(missing_lab)) cat("Labelled genes NOT in du set (skipped):", paste(missing_lab, collapse = ", "), "\n")
cat("\n"); print(sessionInfo())
sink()

message("Done: plot_", slug, ".pdf (+ source_data + stats) in ", out_dir)
