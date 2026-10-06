# ---------------------------------------------------------------------------
# plot_rrho2_ud_top100_annotated.R
#
# Figure panels: 4D
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# =============================================================================
# plot_rrho2_ud_top100_annotated.R
#
# Scatter of the TOP-100 ud discordant genes (up in PHF1+ / down near tangles),
# positioned by their two moderated t-statistics (x = cell-autonomous t, + = up in
# PHF1+; y = proximity t, + = up near tangles) and coloured by the MANUAL functional
# module annotation in rrho2_ud_top100_annotated.csv. Key genes are labelled.
#
# Companion to R/rrho2_decomposition_phf1.r. Reads only existing files (the annotation
# CSV + the t-stat quadrant source data), so it needs neither RRHO2 nor the network.
# Triple output into the same rrho2/ folder.
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
annot_f   <- file.path(rrho2_dir, "rrho2_ud_top100_annotated.csv")
tval_f    <- file.path(rrho2_dir, paste0("source_data_rrho2_tstat_quadrant_", ct, ".tsv"))
out_dir   <- rrho2_dir
slug      <- paste0("rrho2_ud_top100_annotated_", gsub("[^A-Za-z0-9]+", "_", ct))

LABEL_GENES <- c("PSMA5", "UBE2S", "UBE2B", "ATG10", "PICALM", "CHMP2A", "LAPTM4B",
                 "LONP1", "PDPR", "NDUFA2", "NDUFB6", "GRIN2B", "NTRK2", "ADAM22",
                 "RAPGEF4", "PRKRA", "ANK1")

# Module palette: 9 functional modules on Paul Tol's CB-safe 'muted' set; minor ('adhesion',
# 'other') and flagged ('CONTAMINATION_RISK','IMPLAUSIBLE') categories de-emphasised in grey.
module_levels <- c("chromatin", "dna_repair", "proteostasis", "translation", "metabolic",
                   "cell_cycle_brake", "cytoskeleton_actin", "stress_isr", "synaptic_neuronal",
                   "adhesion", "other")
module_pal <- c(
  chromatin          = "#332288",   # indigo
  dna_repair         = "#88CCEE",   # cyan
  proteostasis       = "#CC6677",   # rose
  translation        = "#AA4499",   # purple
  metabolic          = "#DDCC77",   # sand
  cell_cycle_brake   = "#882255",   # wine
  cytoskeleton_actin = "#117733",   # green
  stress_isr         = "#999933",   # olive
  synaptic_neuronal  = "#44AA99",   # teal
  adhesion           = "grey55",
  other              = "grey72")    # 'Other' also absorbs CONTAMINATION_RISK + IMPLAUSIBLE (folded below)
# pretty legend labels (capitalised, underscores -> spaces; acronyms kept upper-case)
module_labels <- c(
  chromatin = "Chromatin", dna_repair = "DNA repair", proteostasis = "Proteostasis",
  translation = "Translation", metabolic = "Metabolic", cell_cycle_brake = "Cell cycle brake",
  cytoskeleton_actin = "Cytoskeleton actin", stress_isr = "Stress ISR",
  synaptic_neuronal = "Synaptic neuronal", adhesion = "Adhesion", other = "Other")

# ------------------------------- load ----------------------------------------
stopifnot(file.exists(annot_f), file.exists(tval_f))
ann <- read.csv(annot_f, stringsAsFactors = FALSE)
tv  <- read.delim(tval_f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

d <- ann %>% left_join(tv[, c("gene", "ca_t", "prox_t")], by = "gene")
miss_t <- d$gene[!is.finite(d$ca_t) | !is.finite(d$prox_t)]
if (length(miss_t)) message("WARNING: no t-values for ", length(miss_t), " gene(s): ",
                            paste(miss_t, collapse = ", "))
d <- d %>% filter(is.finite(ca_t), is.finite(prox_t))
# fold flagged categories into 'Other' for the figure (originals kept in source data)
d$module_grp <- ifelse(d$module %in% c("CONTAMINATION_RISK", "IMPLAUSIBLE", "other"),
                       "other", d$module)
unknown_mod <- setdiff(unique(d$module_grp), module_levels)
if (length(unknown_mod)) { message("NOTE: modules not in palette (folded to 'Other'): ",
                                    paste(unknown_mod, collapse = ", "))
  d$module_grp[d$module_grp %in% unknown_mod] <- "other" }
d$module_grp <- factor(d$module_grp, levels = module_levels)

lab <- d %>% filter(gene %in% LABEL_GENES)
missing_lab <- setdiff(LABEL_GENES, lab$gene)
if (length(missing_lab)) message("NOTE: label genes absent from annotated top-100: ",
                                 paste(missing_lab, collapse = ", "))

# ------------------------------- plot ----------------------------------------
# common axis spans (from BOTH the ud and du gene sets) so the ud/du panels are a matched pair --
# same scale and size, each centred on its own data. ud fills its axes; du (smaller effects) does not.
ud_g <- readLines(file.path(rrho2_dir, paste0("rrho2_genes_ud_top100_", ct, ".txt")))
du_g <- read.csv(file.path(rrho2_dir, "rrho2_du_annotated.csv"))$gene
.span <- function(g, col) diff(range(tv[[col]][match(g, tv$gene)], na.rm = TRUE))
S_x <- max(.span(ud_g, "ca_t"),   .span(du_g, "ca_t"))
S_y <- max(.span(ud_g, "prox_t"), .span(du_g, "prox_t"))
xr <- mean(range(d$ca_t))   + c(-1, 1) * S_x / 2
yr <- mean(range(d$prox_t)) + c(-1, 1) * S_y / 2
p <- ggplot(d, aes(ca_t, prox_t, colour = module_grp)) +
  geom_point(size = 1.3, alpha = 0.85, stroke = 0) +
  ggrepel::geom_text_repel(data = lab, aes(ca_t, prox_t, label = gene), inherit.aes = FALSE,
             size = 2, colour = "grey15", max.overlaps = Inf,
             min.segment.length = 0, segment.size = 0.2, box.padding = 0.3) +
  scale_colour_manual(values = module_pal, labels = module_labels, breaks = module_levels,
                      drop = FALSE, name = "Module") +
  labs(x = expression("Cell-autonomous " * italic(t) * "  (+ = up in PHF1+ vs PHF1-)"),
       y = expression("Proximity " * italic(t) * "  (- = down near tangles)")) +
  coord_fixed(ratio = 1, xlim = xr, ylim = yr) +   # equal unit scaling: 2 units on x == 2 units on y
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_classic(base_size = 8) + fig_theme +   # house standard (axis.title 8, legend 6/7 via fig_theme)
  theme(legend.position = "right", legend.key.size = grid::unit(7, "pt"),
        plot.margin = margin(2, 2, 2, 2))

# canvas sized so the square (coord_fixed) panel fills it: device width = panel + y-title + legend,
# device height = panel + x-title, i.e. width ~= height + legend/title width -> no reserved margin.
ggsave(file.path(out_dir, paste0("plot_", slug, ".pdf")),
       p, width = 9.18, height = 6.0, units = "cm", device = "pdf")

# ------------------------------- source data + stats -------------------------
sd <- d[, c("gene", "module", "ca_t", "prox_t", "annotation", "flag")]
sd <- sd[order(sd$module, -sd$ca_t), ]
write.table(sd, file.path(out_dir, paste0("source_data_", slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug, ".txt")))
cat("Top-100 ud discordant genes (up in PHF1+ / down near tangles) --", ct, "\n")
cat("x = cell-autonomous t (limma; + up in PHF1+); y = proximity t (dream, negated; + up near tangles)\n")
cat("Colour = manual functional module (rrho2_ud_top100_annotated.csv).\n\n")
cat("n plotted:", nrow(d), " (of", nrow(ann), "annotated)\n")
if (length(miss_t)) cat("dropped (no t-value):", paste(miss_t, collapse = ", "), "\n")
cat("\nGenes per module:\n"); print(as.data.frame(table(module = d$module)))
cat("\nFlag tally:\n"); print(as.data.frame(table(flag = d$flag)))
cat("\nSpearman r(ca_t, prox_t) over the top-100:",
    round(suppressWarnings(cor(d$ca_t, d$prox_t, method = "spearman")), 3), "\n")
cat("\nLabelled genes present:", paste(sort(lab$gene), collapse = ", "), "\n")
if (length(missing_lab)) cat("Labelled genes NOT found:", paste(missing_lab, collapse = ", "), "\n")
cat("\n"); print(sessionInfo())
sink()

message("Done: plot_", slug, ".pdf (+ source_data + stats) in ", out_dir)
