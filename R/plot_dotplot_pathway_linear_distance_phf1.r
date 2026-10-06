# ---------------------------------------------------------------------------
# plot_dotplot_pathway_linear_distance_phf1.r
#
# Figure panels: 3B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# =====================================================================
# Dotplot of distance-DEG enrichR pathway enrichment, ONE celltype.
#
# Reads the per-celltype enrichR ORA results produced downstream of the
# canonical cell-level linear distance-to-PHF1 DEG
# (R/deg_dream_linear_distance_phf1.r), and draws a single combined dotplot:
#   y axis  = pathway (description), grouped into categories (one facet per
#             group, with a small break between groups)
#   x axis  = signed odds ratio  (UP pathways +, DOWN pathways -)
#   size    = overlap  (n genes from the DEG list hitting the pathway)
#   colour  = FDR
#
# Direction note (distance DEG): the enrichR UP/DOWN folders encode the
# distance-coefficient sign (NOT PHF1+/-): UP = up with distance (higher FAR
# from PHF1, lower near); DOWN = down with distance (higher NEAR PHF1).
# This plot is drawn in PROXIMITY terms to match nDEG_by_celltype_linear.R: the x-axis sign is
# FLIPPED so positive odds ratio = UPREGULATED NEAR NFTs (enrichR DOWN folder)
# and negative = downregulated near NFTs (enrichR UP folder).
#
# RUN INTERACTIVELY, ONE CELLTYPE AT A TIME:
#   1. Edit the CONFIG block, then source PHASE 1. It prints the candidate
#      pathway table and writes a selection TSV with a `plot` column (all 0).
#   2. Open that TSV and code the `plot` column:
#        (blank) = exclude;  1, 2, 3, ... = the category group to plot it in.
#      Name the categories via `category_labels` in the CONFIG block. Save.
#   3. Source PHASE 2 to draw the triple (PDF + source_data TSV + stats TXT).
#
# Only 9 celltypes have distance enrichR output:
#   Astro, Endo, Exc-IT-L2-3-CBLN2-HOPX, Exc-IT-L3-5-CHGA-IL1RAPL2, Micro,
#   OPC, Oligo, Unassigned Neuron, Unassigned
# =====================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(stringr); library(scales)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # exports: celltype_palette, braak_palette, neuron_order, glia_order, braak_levels, fig_theme

# ---------------------------------------------------------------------
# CONFIG  -- edit per run
# ---------------------------------------------------------------------
celltype        <- "Exc-IT-L2-3-CBLN2-HOPX"   # set ONE celltype (see header for the 9 valid options)
fdr_cap         <- 0.05      # FDR colour-scale upper limit (log10 scale; values above are squished)
plot_w_cm       <- 9.38      # fixed plot width (cm)
row_cm          <- 0.42      # per-pathway height; total height auto-scales with n selected
wrap_width      <- 42        # wrap pathway labels to this many characters
reset_selection <- FALSE    # TRUE to (re)write/overwrite an existing selection TSV; FALSE keeps your edits

# Category grouping: code the `plot` column of the selection TSV with integers
# (blank = exclude; 1, 2, 3, ... = category group). Map those codes -> the group
# labels shown on the plot here. Codes with no entry fall back to "Group <n>".
category_labels <- c(
  "1" = "Category 1",
  "2" = "Category 2",
  "3" = "Category 3"
)
cat_wrap <- 16               # wrap category group labels to this many characters

# Overlap (size) scale -- FIXED across celltypes so a given overlap maps to the
# same circle size in every plot (global overlap span across the distance
# celltypes is 3-10). Bump the upper limit if a future enrichR rerun exceeds it.
overlap_limits <- c(3, 10)
overlap_breaks <- c(3, 5, 10)

slug    <- paste0("pathway_linear_distance_", gsub("[^A-Za-z0-9]+", "_", celltype))
out_dir <- file.path("plots/dotplot_pathway_linear_distance", celltype)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sel_file <- file.path(out_dir, paste0("selection_", slug, ".tsv"))

# Columns common to both the *_merged.tsv files and the per-database files.
KEEP_COLS <- c("geneset", "description", "database",
               "size", "overlap", "odds_ratio", "pval", "FDR", "genes")

# Read a celltype/direction enrichR result. Prefers the *_merged.tsv (all
# databases combined); falls back to rbind-ing the per-database TSVs.
read_enrichr_dir <- function(ct, direction) {
  d <- file.path("deg/de_linear_distance", ct, "enrichr", direction)
  if (!dir.exists(d)) return(NULL)
  merged <- file.path(d, "merged_enrichr.tsv")
  if (file.exists(merged)) {
    df <- read.delim(merged, sep = "\t", header = TRUE, quote = "",
                     stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    fs <- list.files(d, pattern = "\\.tsv$", full.names = TRUE)
    fs <- fs[!grepl("merged_enrichr\\.tsv$", fs)]
    if (length(fs) == 0) return(NULL)
    df <- do.call(rbind, lapply(fs, function(f)
      read.delim(f, sep = "\t", header = TRUE, quote = "",
                 stringsAsFactors = FALSE, check.names = FALSE)[, KEEP_COLS, drop = FALSE]))
  }
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df[, KEEP_COLS, drop = FALSE]
  df$direction <- direction
  df
}

# =====================================================================
# PHASE 1 -- build candidate table + write selection TSV
# =====================================================================
up <- read_enrichr_dir(celltype, "UP")
dn <- read_enrichr_dir(celltype, "DOWN")

if (is.null(up) && is.null(dn)) {
  stop("No distance enrichR output for celltype '", celltype,
       "'. Valid: Astro, Endo, Exc-IT-L2-3-CBLN2-HOPX, Exc-IT-L3-5-CHGA-IL1RAPL2, ",
       "Micro, OPC, Oligo, Unassigned Neuron, Unassigned.")
}

cand <- bind_rows(up, dn) %>%
  mutate(
    description      = trimws(description),
    celltype         = celltype,
    # PROXIMITY sign (flipped vs the distance coefficient, matches nDEG_by_celltype_linear.R):
    # enrichR DOWN = up near NFTs -> positive; UP = down near NFTs -> negative.
    signed_oddsratio = ifelse(direction == "UP", -odds_ratio, odds_ratio),
    plot             = ""
  ) %>%
  arrange(direction, FDR) %>%
  select(plot, celltype, direction, description, database, geneset,
         overlap, odds_ratio, signed_oddsratio, FDR, pval, genes)

message("\n", celltype, ": ", sum(cand$direction == "UP"), " UP / ",
        sum(cand$direction == "DOWN"), " DOWN candidate pathways.\n")
print(cand[, c("plot", "direction", "description", "overlap", "odds_ratio", "FDR")])

if (!file.exists(sel_file) || reset_selection) {
  write.table(cand, sel_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("\nWrote selection template:\n  ", sel_file,
          "\n-> Code the `plot` column (blank = exclude; 1, 2, 3, ... = category",
          " group), SAVE, then run PHASE 2.\n")
} else {
  message("\nSelection file already exists (edits preserved):\n  ", sel_file,
          "\n-> Code the `plot` column (blank = exclude; 1, 2, 3, ... = category",
          " group), SAVE, then run PHASE 2.",
          "\n   Set reset_selection <- TRUE to regenerate a blank template.\n")
}

# ===========================  >>> STOP <<<  ==========================
# Edit `sel_file`: in the `plot` column leave blank to exclude a pathway, or
# put an integer category code (1, 2, 3, ...) to plot it in that group. Save,
# then source everything below. Name the groups via `category_labels` (CONFIG).
# =====================================================================

# =====================================================================
# PHASE 2 -- read selection + draw triple
# =====================================================================
sel_all <- read.delim(sel_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)

# Robustness: a selection file edited in a spreadsheet can pick up a duplicate
# `plot` column. Coalesce any/all `plot` columns into one, taking the FIRST
# positive integer category code found per row (blank / 0 = exclude).
plot_cols <- which(names(sel_all) == "plot")
if (length(plot_cols) > 1) {
  coalesced <- apply(sel_all[, plot_cols, drop = FALSE], 1, function(v) {
    iv <- suppressWarnings(as.integer(trimws(as.character(v))))
    iv <- iv[!is.na(iv) & iv > 0]
    if (length(iv)) iv[1] else NA_integer_
  })
  sel_all <- sel_all[, -plot_cols, drop = FALSE]
  sel_all$plot <- coalesced
  message("Note: selection file had ", length(plot_cols),
          " `plot` columns; merged them (first positive category code per row wins).")
}

# `plot` column codes the category: blank/NA = exclude, positive integer = group.
sel_all$cat_id <- suppressWarnings(as.integer(sel_all$plot))
sel <- sel_all[!is.na(sel_all$cat_id) & sel_all$cat_id > 0, , drop = FALSE]

if (nrow(sel) == 0) {
  stop("No pathways selected in ", sel_file,
       " -- set the `plot` column to a positive integer category code on at ",
       "least one row (blank = exclude), save, and re-run PHASE 2.")
}

# Category groups: ordered (facets, top -> bottom) by the MAX odds ratio within
# each category -- strongest-enrichment group at the top -- not by the numeric
# code. `category` defines the facet groups (labels not drawn on the plot);
# `category_name` (mapped via category_labels, fallback "Group <n>") goes to the
# source-data table.
cat_max_or <- tapply(sel$odds_ratio, sel$cat_id, max, na.rm = TRUE)
cat_ids    <- as.integer(names(sort(cat_max_or, decreasing = TRUE)))
cat_nm  <- ifelse(as.character(cat_ids) %in% names(category_labels),
                  category_labels[as.character(cat_ids)], paste0("Group ", cat_ids))
sel$category_name <- cat_nm[match(sel$cat_id, cat_ids)]
sel$category <- factor(sel$cat_id, levels = cat_ids, labels = str_wrap(cat_nm, cat_wrap))

# Order rows in display order (category by max-OR, then signed odds ratio).
sel <- sel[order(match(sel$cat_id, cat_ids), sel$signed_oddsratio), , drop = FALSE]

# Wrapped y-axis label, ordered along the signed odds-ratio axis within each group.
sel$label <- str_wrap(sel$description, wrap_width)
sel$label <- reorder(sel$label, sel$signed_oddsratio)

# FDR colour scale is log10 (never reaches 0). The lower limit is rounded DOWN
# to a clean power of 10 so the bottom of the bar is a round number. Breaks
# ALWAYS include both ends (fdr_lo and fdr_cap) so readers can read the extremes,
# plus the decade ticks in between.
fdr_min <- suppressWarnings(min(sel$FDR[sel$FDR > 0], na.rm = TRUE))
if (!is.finite(fdr_min) || fdr_min >= fdr_cap) fdr_min <- fdr_cap / 10
fdr_lo  <- 10^floor(log10(fdr_min))
fdr_inner  <- 10^seq(floor(log10(fdr_lo)) + 1, ceiling(log10(fdr_cap)) - 1)
fdr_inner  <- fdr_inner[fdr_inner > fdr_lo & fdr_inner < fdr_cap]
fdr_breaks <- sort(unique(c(fdr_lo, fdr_inner, fdr_cap)))
fdr_fmt    <- function(x) formatC(x, format = "fg", digits = 2, drop0trailing = TRUE)

# Size MAPPING stays fixed (overlap_limits) so circles are comparable across
# celltypes, but only show the fixed breaks that ENCOMPASS this plot's overlap
# range -- bracket from the last break <= min to the first break >= max.
ov_min <- min(sel$overlap); ov_max <- max(sel$overlap)
lo_i <- suppressWarnings(max(which(overlap_breaks <= ov_min)))
hi_i <- suppressWarnings(min(which(overlap_breaks >= ov_max)))
if (!is.finite(lo_i)) lo_i <- 1L
if (!is.finite(hi_i)) hi_i <- length(overlap_breaks)
size_breaks <- overlap_breaks[lo_i:hi_i]

p <- ggplot(sel, aes(x = signed_oddsratio, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
  geom_point(aes(fill = FDR, size = overlap),
             shape = 21, colour = "black", stroke = 0.3, alpha = 0.9) +
  # Size mapping FIXED via limits (circles comparable across celltypes); legend
  # shows only the breaks that encompass this plot's range (size_breaks).
  # order = 1 keeps this legend above the FDR bar on every plot.
  scale_size(name = "Overlap", range = c(2, 6), limits = overlap_limits,
             breaks = size_breaks, guide = guide_legend(order = 1)) +
  scale_fill_gradient(low = "navy", high = "gold", name = "FDR", transform = "log10",
                      limits = c(fdr_lo, fdr_cap), oob = scales::squish,
                      breaks = fdr_breaks, labels = fdr_fmt,
                      guide = guide_colourbar(
                        reverse = TRUE, order = 2,
                        theme = theme(legend.key.height = grid::unit(1.6, "cm"),
                                      legend.key.width  = grid::unit(0.30, "cm")))) +
  # Group pathways into categories: one panel per group, each sized to its own
  # n pathways (space = "free_y"), with a small break between groups. Categories
  # are NOT labelled on the plot -- the y axis shows only pathway names.
  facet_grid(rows = vars(category), scales = "free_y", space = "free_y") +
  labs(title = celltype, x = "Odds ratio (up near NFTs +, down near NFTs -)", y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(plot.title       = element_text(size = 8, face = "bold", hjust = 0.5),  # celltype only
        axis.text.x      = element_text(angle = 0, hjust = 0.5, vjust = 1),  # numeric x
        strip.background = element_blank(),
        strip.text       = element_blank(),  # group, but don't label the categories
        panel.spacing.y  = grid::unit(4.5, "pt"),  # the small break between groups
        # override fig_theme's wide left margin (meant for rotated x labels) so
        # there is no empty gap at the left of the figure.
        plot.margin      = margin(t = 4, r = 6, b = 4, l = 2, unit = "pt"))

# Height: per-pathway rows + a little overhead per group (panel break).
h_cm <- max(4, nrow(sel) * row_cm + length(cat_ids) * 0.35 + 2.2)

# --- TRIPLE OUTPUT -------------------------------------------------
# 1. vector PDF
ggsave(file.path(out_dir, paste0("plot_", slug, ".pdf")), p,
       width = plot_w_cm, height = h_cm, units = "cm", device = "pdf")

# 2. source data -- exact rows that produced the dots (long-format)
source_df <- sel %>%
  select(celltype, category_code = cat_id, category = category_name, direction,
         description, database, geneset, overlap, odds_ratio, signed_oddsratio,
         FDR, pval, genes)
write.table(source_df, file.path(out_dir, paste0("source_data_", slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 3. human-readable stats log
sink(file.path(out_dir, paste0("stats_", slug, ".txt")))
cat("Dotplot -- distance-DEG (linear distance-to-PHF1) enrichR pathway enrichment\n")
cat("Celltype:", celltype, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("enrichR folders (distance sign): UP = up with distance (down near NFTs);",
    " DOWN = down with distance (up near NFTs).\n")
cat("PLOTTED in PROXIMITY terms (sign flipped, matches nDEG_by_celltype_linear.R):",
    " positive signed_oddsratio = UPREGULATED NEAR NFTs.\n")
cat("Inputs:\n")
cat("  ", file.path("deg/de_linear_distance", celltype, "enrichr/UP/merged_enrichr.tsv"), "\n")
cat("  ", file.path("deg/de_linear_distance", celltype, "enrichr/DOWN/merged_enrichr.tsv"), "\n")
cat("  selection:", sel_file, "\n\n")
cat("Candidate pathways: ", nrow(sel_all),
    "  (UP ", sum(sel_all$direction == "UP"),
    " / DOWN ", sum(sel_all$direction == "DOWN"), ")\n", sep = "")
cat("Selected (plotted): ", nrow(sel),
    "  (UP ", sum(sel$direction == "UP"),
    " / DOWN ", sum(sel$direction == "DOWN"), ")\n", sep = "")
cat("Categories (", length(cat_ids), ", ordered top->bottom by max odds ratio):\n", sep = "")
for (i in seq_along(cat_ids)) {
  cat("  [", cat_ids[i], "] ", cat_nm[i], "  (n = ", sum(sel$cat_id == cat_ids[i]),
      ", max OR = ", round(max(sel$odds_ratio[sel$cat_id == cat_ids[i]]), 2), ")\n", sep = "")
}
cat("\nAesthetics: y = pathway (grouped into categories via facets), ",
    "x = signed odds ratio (+ up near NFTs / - down near NFTs), size = overlap, colour = FDR.\n", sep = "")
cat("Plot size: ", plot_w_cm, " x ", round(h_cm, 2), " cm.\n\n", sep = "")
cat("Plotted pathways:\n")
print(source_df[, c("category_code", "category", "direction", "description",
                     "overlap", "odds_ratio", "signed_oddsratio", "FDR")])
cat("\n")
print(sessionInfo())
sink()

message("Done. Triple written to: ", out_dir)
