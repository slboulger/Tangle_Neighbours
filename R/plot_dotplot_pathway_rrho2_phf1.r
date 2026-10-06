# ---------------------------------------------------------------------------
# plot_dotplot_pathway_rrho2_phf1.r
#
# Figure panels: S6
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# =====================================================================
# Dotplot of RRHO2 quadrant enrichR pathway enrichment, ONE celltype.
#
# Fork of R/plot_dotplot_pathway_cat_phf1.r that plots the GO tables produced by
# R/rrho2_decomposition_phf1.r (plots/decompose_tangle_response_phf1/{ct}/rrho2/
# go_rrho2_{uu,dd,ud,du}_{ct}.tsv). Same grammar and two-phase workflow:
#   y axis  = pathway (description), grouped into categories (facets)
#   x axis  = signed odds ratio (sign set per contrast, see CONFIG)
#   size    = overlap ; colour = FDR
#
# You CHOOSE which RRHO2 contrasts to plot and their axis sign in the CONFIG
# `contrasts` list (each RRHO2 quadrant is a contrast):
#   uu = up in PHF1+  & up near tangles     (shared/convergent, up)
#   dd = down in PHF1+ & down near tangles   (shared/convergent, down)
#   ud = up in PHF1+  & down near tangles    (discordant)
#   du = down in PHF1+ & up near tangles     (discordant)
#
# RUN INTERACTIVELY, ONE CELLTYPE AT A TIME:
#   1. Edit CONFIG (celltype + contrasts), source PHASE 1 -> prints candidates,
#      writes a selection TSV with a blank `plot` column.
#   2. Code the `plot` column (blank = exclude; 1,2,3,... = category group). Save.
#   3. Source PHASE 2 -> draws a SEPARATE dotplot per group (concordant = uu/dd,
#      discordant = ud/du), reusing your `plot`-column annotations. Each group axis is
#      signed only if both its directions are present, else one-sided positive OR.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(stringr); library(scales)
})
# HPC canonical path, local-mount fallback (runs in either)
proj_hpc   <- "<PROJECT_ROOT>/phf1_v2"
proj_local <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(proj_hpc)) proj_hpc else proj_local)
source("R/palettes.R")  # celltype_palette, fig_theme, ...

# ---------------------------------------------------------------------
# CONFIG  -- edit per run
# ---------------------------------------------------------------------
celltype        <- "Exc-IT-L2-3-CBLN2-HOPX"   # one subtype with an rrho2/ folder
rrho2_dir       <- file.path("plots/decompose_tangle_response_phf1", celltype, "rrho2")

# CHOOSE WHICH CONTRASTS TO PLOT. Each entry = one RRHO2 quadrant, its axis `sign`
# (+1 plotted right, -1 plotted left) and a `direction` label. Comment out any you
# don't want. Default = the concordant/shared program (uu right, dd left).
contrasts <- list(
  list(quad = "uu", sign =  1, direction = "up in both"),
  list(quad = "dd", sign = -1, direction = "down in both"), 
  list(quad = "ud_top100", sign =  1, direction = "up PHF1+/down near"),  # discordant (top-100 ORA)
  list(quad = "du", sign = -1, direction = "down PHF1+/up near")   # discordant
)
x_lab           <- "Odds ratio (up in both +, down in both -)"  # match your chosen contrasts

fdr_cap         <- 0.05      # FDR colour-scale upper limit (log10; values above squished)
plot_w_cm       <- 9.38      # fixed plot width (cm)
row_cm          <- 0.42      # per-pathway height; total auto-scales
wrap_width      <- 34        # wrap pathway labels
reset_selection <- FALSE     # TRUE to overwrite an existing selection TSV; FALSE keeps edits

category_labels <- c("1" = "Category 1", "2" = "Category 2", "3" = "Category 3")
cat_wrap <- 16

overlap_limits <- c(3, 40)   # FIXED size scale (comparable across figures) -- matches cat dotplot
overlap_breaks <- c(3, 5, 10, 20, 30)

slug    <- paste0("pathway_rrho2_phf1_", gsub("[^A-Za-z0-9]+", "_", celltype))
out_dir <- file.path("plots/dotplot_pathway_rrho2_phf1", celltype)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sel_file <- file.path(out_dir, paste0("selection_", slug, ".tsv"))

KEEP_COLS <- c("geneset", "description", "database",
               "size", "overlap", "odds_ratio", "pval", "FDR", "genes")

# Read one RRHO2 quadrant's GO table; tag it with the direction label and sign.
read_rrho2_contrast <- function(ct, quad, direction, sign) {
  f <- file.path(rrho2_dir, paste0("go_rrho2_", quad, "_", ct, ".tsv"))
  if (!file.exists(f)) { message("  (no GO table for quadrant '", quad, "': ", f, ")"); return(NULL) }
  df <- read.delim(f, sep = "\t", header = TRUE, quote = "",
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df[, KEEP_COLS, drop = FALSE]
  df$direction <- direction
  df$quad      <- quad
  df$sign      <- sign
  df
}

# =====================================================================
# PHASE 1 -- build candidate table + write selection TSV
# =====================================================================
cand_list <- lapply(contrasts, function(cc)
  read_rrho2_contrast(celltype, cc$quad, cc$direction, cc$sign))
cand_list <- Filter(Negate(is.null), cand_list)

if (length(cand_list) == 0)
  stop("No RRHO2 GO tables found for '", celltype, "' in ", rrho2_dir,
       ". Run R/rrho2_decomposition_phf1.r first, and check the `contrasts` quadrants exist.")

cand <- bind_rows(cand_list) %>%
  mutate(
    description      = trimws(description),
    celltype         = celltype,
    signed_oddsratio = sign * odds_ratio,
    plot             = ""
  ) %>%
  arrange(direction, FDR) %>%
  select(plot, celltype, direction, description, database, geneset,
         overlap, odds_ratio, signed_oddsratio, FDR, pval, genes)

message("\n", celltype, " -- candidate pathways by contrast:")
print(as.data.frame(table(cand$direction)))
print(cand[, c("plot", "direction", "description", "overlap", "odds_ratio", "FDR")])

if (!file.exists(sel_file) || reset_selection) {
  write.table(cand, sel_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("\nWrote selection template:\n  ", sel_file,
          "\n-> Code the `plot` column (blank = exclude; 1,2,3,... = category group),",
          " SAVE, then run PHASE 2.\n")
} else {
  # MERGE: keep existing `plot` codes, refresh the stats columns, and APPEND any NEW candidate rows
  # (e.g. a newly added contrast such as ud_top100) as blank -- so you never re-annotate coded rows.
  # Match on (direction, geneset): geneset is the stable GO id. Stale rows (no longer a candidate,
  # e.g. dropped by a min_size change) fall out.
  old <- read.delim(sel_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  pc  <- which(names(old) == "plot")
  old_plot <- if (length(pc) > 1)
    apply(old[, pc, drop = FALSE], 1, function(v) {
      iv <- suppressWarnings(as.integer(trimws(as.character(v)))); iv <- iv[!is.na(iv) & iv > 0]
      if (length(iv)) as.character(iv[1]) else "" })
    else as.character(old$plot)
  old_plot[is.na(old_plot)] <- ""
  old_key  <- paste(old$direction, old$geneset, sep = "\t")
  cand_key <- paste(cand$direction, cand$geneset, sep = "\t")
  cand$plot <- old_plot[match(cand_key, old_key)]
  cand$plot[is.na(cand$plot)] <- ""
  n_new     <- sum(!(cand_key %in% old_key))
  n_dropped <- sum(!(old_key %in% cand_key))
  write.table(cand, sel_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("\nSelection file merged (existing codes preserved):\n  ", sel_file,
          "\n  ", n_new, " new candidate row(s) added blank; ", n_dropped, " stale row(s) dropped.",
          "\n-> Code any new blank rows, SAVE, run PHASE 2. Set reset_selection <- TRUE to fully regenerate.\n")
}

# ===========================  >>> STOP <<<  ==========================
# Edit `sel_file`: `plot` column blank = exclude, integer 1,2,3,... = category group.
# Save, then source everything below. Name groups via `category_labels` (CONFIG).
# =====================================================================

# =====================================================================
# PHASE 2 -- read selection + draw triple
# =====================================================================
sel_all <- read.delim(sel_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE)

# Coalesce duplicate `plot` columns a spreadsheet may introduce (first positive code wins).
plot_cols <- which(names(sel_all) == "plot")
if (length(plot_cols) > 1) {
  coalesced <- apply(sel_all[, plot_cols, drop = FALSE], 1, function(v) {
    iv <- suppressWarnings(as.integer(trimws(as.character(v))))
    iv <- iv[!is.na(iv) & iv > 0]; if (length(iv)) iv[1] else NA_integer_
  })
  sel_all <- sel_all[, -plot_cols, drop = FALSE]; sel_all$plot <- coalesced
  message("Note: merged ", length(plot_cols), " `plot` columns (first positive code per row).")
}

sel_all$cat_id <- suppressWarnings(as.integer(sel_all$plot))
sel <- sel_all[!is.na(sel_all$cat_id) & sel_all$cat_id > 0, , drop = FALSE]

if (nrow(sel) == 0)
  stop("No pathways selected in ", sel_file,
       " -- code the `plot` column with a positive integer on >=1 row, save, re-run PHASE 2.")

# Map each direction label -> quadrant, group (concordant = uu/dd, discordant = ud/du),
# and within-group sign (up-first quadrant = +, down-first = -). Uses the CONFIG `contrasts`
# as the source of truth, so YOUR existing `plot`-column annotations are reused unchanged.
.quad_by_dir <- setNames(vapply(contrasts, `[[`, "", "quad"),
                         vapply(contrasts, `[[`, "", "direction"))
.wsign       <- c(uu = 1, dd = -1, ud = 1, ud_top100 = 1, du = -1)
# Axis wording per quadrant (decoupled from the stored `direction` label so the selection
# file's annotations stay valid). Edit these to taste.
axis_dir     <- c(uu = "up PHF1+/up near", dd = "down PHF1+/down near",
                  ud = "up PHF1+/down near", ud_top100 = "up PHF1+/down near",
                  du = "down PHF1+/up near")
sel$quad  <- .quad_by_dir[sel$direction]
sel$group <- ifelse(sel$quad %in% c("uu", "dd"), "concordant", "discordant")

# Draw ONE group's dotplot (a separate figure). One-sided positive OR when only one direction
# is present in the group (avoids a lone discordant set reading as "down"); signed when both.
draw_group <- function(sub, group_name) {
  quads_present <- unique(sub$quad)
  one_sided     <- length(unique(.wsign[quads_present])) == 1
  sub$signed_oddsratio <- if (one_sided) sub$odds_ratio else .wsign[sub$quad] * sub$odds_ratio
  if (one_sided) {
    x_lab_g <- paste0("Odds ratio  (", axis_dir[quads_present[1]], ")")
  } else {
    pos_q <- quads_present[.wsign[quads_present] > 0][1]
    neg_q <- quads_present[.wsign[quads_present] < 0][1]
    x_lab_g <- paste0("Odds ratio  (", axis_dir[pos_q], " +, ", axis_dir[neg_q], " -)")
  }
  # categories ordered top->bottom by max OR within this group
  cmo  <- tapply(sub$odds_ratio, sub$cat_id, max, na.rm = TRUE)
  cids <- as.integer(names(sort(cmo, decreasing = TRUE)))
  cnm  <- ifelse(as.character(cids) %in% names(category_labels),
                 category_labels[as.character(cids)], paste0("Group ", cids))
  sub$category_name <- cnm[match(sub$cat_id, cids)]
  sub$category <- factor(sub$cat_id, levels = cids, labels = str_wrap(cnm, cat_wrap))
  sub <- sub[order(match(sub$cat_id, cids), sub$signed_oddsratio), , drop = FALSE]
  sub$label <- reorder(str_wrap(sub$description, wrap_width), sub$signed_oddsratio)
  # FDR colour scale
  fmin <- suppressWarnings(min(sub$FDR[sub$FDR > 0], na.rm = TRUE))
  if (!is.finite(fmin) || fmin >= fdr_cap) fmin <- fdr_cap / 10
  flo  <- 10^floor(log10(fmin))
  finn <- 10^seq(floor(log10(flo)) + 1, ceiling(log10(fdr_cap)) - 1)
  finn <- finn[finn > flo & finn < fdr_cap]
  fbrk <- sort(unique(c(flo, finn, fdr_cap)))
  # size breaks (fixed mapping, show only those encompassing this plot)
  li  <- suppressWarnings(max(which(overlap_breaks <= min(sub$overlap))))
  hii <- suppressWarnings(min(which(overlap_breaks >= max(sub$overlap))))
  if (!is.finite(li))  li  <- 1L
  if (!is.finite(hii)) hii <- length(overlap_breaks)
  sbrk <- overlap_breaks[li:hii]

  p <- ggplot(sub, aes(x = signed_oddsratio, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    geom_point(aes(fill = FDR, size = overlap), shape = 21, colour = "black",
               stroke = 0.3, alpha = 0.9) +
    scale_size(name = "Overlap", range = c(2, 6), limits = overlap_limits,
               breaks = sbrk, guide = guide_legend(order = 1)) +
    scale_fill_gradient(low = "navy", high = "gold", name = "FDR", transform = "log10",
                        limits = c(flo, fdr_cap), oob = scales::squish, breaks = fbrk,
                        labels = function(x) formatC(x, format = "fg", digits = 2, drop0trailing = TRUE),
                        guide = guide_colourbar(reverse = TRUE, order = 2,
                          theme = theme(legend.key.height = grid::unit(1.6, "cm"),
                                        legend.key.width  = grid::unit(0.30, "cm")))) +
    facet_grid(rows = vars(category), scales = "free_y", space = "free_y") +
    labs(title = paste0(celltype, " -- ", group_name), x = x_lab_g, y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          strip.background = element_blank(), strip.text = element_blank(),
          panel.spacing.y = grid::unit(4.5, "pt"),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 2, unit = "pt"))
  h_cm  <- max(4, nrow(sub) * row_cm + length(cids) * 0.35 + 2.2)
  gslug <- paste0(slug, "_", group_name)
  ggsave(file.path(out_dir, paste0("plot_", gslug, ".pdf")), p,
         width = plot_w_cm, height = h_cm, units = "cm", device = "pdf")
  sdf <- sub %>% select(celltype, category_code = cat_id, category = category_name,
                        direction, quad, description, database, geneset, overlap,
                        odds_ratio, signed_oddsratio, FDR, pval, genes)
  write.table(sdf, file.path(out_dir, paste0("source_data_", gslug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c(
    paste0("RRHO2 dotplot -- ", celltype, " -- ", group_name),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("x axis: ", x_lab_g, "  (one-sided = ", one_sided, ")"),
    paste0("quadrants present: ", paste(quads_present, collapse = ", ")),
    paste0("n pathways: ", nrow(sub), " in ", length(cids), " categories"),
    "size = overlap; colour = FDR (log10, capped at fdr_cap)."),
    file.path(out_dir, paste0("stats_", gslug, ".txt")))
  message("  ", group_name, ": ", nrow(sub), " pathways -> plot_", gslug, ".pdf")
}

# concordant first; guard each group so an externally-locked output file (e.g. a TSV open
# in a viewer -> "Resource temporarily unavailable") can't abort the other group's plot.
group_order <- intersect(c("concordant", "discordant"), unique(sel$group))
for (g in group_order) {
  tryCatch(draw_group(sel[sel$group == g, , drop = FALSE], g),
           error = function(e) message("  [", g, "] FAILED: ", conditionMessage(e)))
}
message("Done. Separate per-group dotplots written to: ", out_dir)
