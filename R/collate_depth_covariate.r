#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# collate_depth_covariate.r
#
# Figure panels: S1B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# collate_depth_covariate.r
#
# Compares the depth-adjusted DEG (deg_dream_linear_distance_depth_phf1.r,
# ~ ... + depth_s + ... + (1|sample_id)) against the canonical unadjusted parent
# (deg_dream_linear_distance_phf1.r, same model without depth_s), and asks
# whether the tangle-distance coefficients survive conditioning on where in the
# section a cell sits.
#
# This is the celltype-facetted analogue of collate_exclusion_radius.r. Both share
# the same helper definitions so the numbers are computed identically and are
# directly comparable. There is ONE model contrast (+ depth_s) but many
# celltypes, so the sweep dimension - and the facet - is CELLTYPE.
#
# Reads, per celltype:
#   <deg_dir>/de_linear_distance/<ct>/<ct>_dist_to_phf1_um_scaled.tsv         (parent)
#   <deg_dir>/de_linear_distance_depth/<ct>/<ct>_dist_to_phf1_um_scaled.tsv   (depth)
#   <deg_dir>/de_linear_distance_depth/<ct>/<ct>_depth_summary.tsv            (optional
#     context: how much of the distance predictor depth_s shares, how many genes
#     respond to depth at all, and the laminar-ordinal correlation of the axis -
#     joined into the summary so the concordance sits next to it)
#
# WHAT depth_s IS
# ---------------
# depth_s is the standardised position of each cell along its sample's manual
# pia->white-matter compass axis. laminar_rho_pooled (the Spearman correlation
# between the laminar excitatory ordinal L2-3 -> L3-5 -> L5 -> L6-IT -> L6-CT and
# depth_rel) is carried into the summary TSV alongside the numbers.
#
# WHICH COEFFICIENT IS COMPARABLE. logFC is logCPM per SD of transformed distance,
# and dist_sd is recomputed per run from the retained cells. Both models here are
# fitted on the SAME cells (depth_rel is defined for every cell, so adding depth_s
# cannot drop one), so dist_sd should be IDENTICAL and logFC directly comparable.
# That is not assumed - the DEG script asserts it at write time, and it is
# verified again below by recovering dist_sd as logFC / logFC_per_unit from each
# side. The plots and stats nonetheless use logFC_per_unit (logCPM per unit
# log(um)), both because it is the scale-free quantity and because it matches the
# axis of the exclusion-radius concordance figures.
#
# USE THE WITHIN-DONOR OVERLAP, NOT THE POOLED ONE
# ------------------------------------------------
# Both models carry (1|sample_id), so they condition on DONOR-DEMEANED
# predictors. The overlap that governs how much depth_s can move the distance
# coefficient is therefore the within-donor correlation, `depth_r2_within`, not
# the pooled `depth_r2_dist`. Measured on this cohort the pooled figure
# understates the within-donor R^2 by 2x to >1000x, and for Oligo it has the
# OPPOSITE SIGN (+0.052 pooled vs -0.078 within) - a Simpson's reversal across
# donors. The pooled columns are retained for reference only; quoting them as
# "how much of the predictor depth can remove" is wrong.
#
# HOW TO READ AN ATTENUATED SLOPE. depth_s and the distance
# predictor are both cell-level fixed effects and they share variance. Adjusting
# for depth therefore removes part of the predictor, and a shrunken coefficient is
# expected on that ground alone. It is NOT by itself evidence that the gradient
# was a positional confound. The discriminating evidence is the pairing of:
#   (a) the Deming slope here (how much the coefficients moved), against
#   (b) depth_r2_within (how much of the predictor depth could remove at all) and
#       n_genes_sig_for_depth (whether depth_s explains any expression variance,
#       or is an inert column that cost a degree of freedom), and
#   (c) the sign-flip count (a confound would scramble signs; loss of shared
#       variance alone shrinks magnitudes toward zero without flipping them).
# All three are emitted here.
#
# Genes that change sign or lose significance under this adjustment are flagged
# per gene in the source data (sign_flip, sig_class).
#
# Deliverables (plotting triple: PDF + source_data TSV + stats TXT):
#   1. plot_depth_concordance_allgenes_parentsig.pdf -- all genes, facet per
#      celltype, points coloured ONLY by significance in the PARENT model. The
#      figure for "no parent-significant gene changes sign": no vermillion
#      point falls in an off-diagonal quadrant. Threshold crossings are ignored
#      here on purpose so nothing distracts from the sign claim.
#   2. plot_depth_concordance_allgenes_bysig.pdf -- all genes, coloured by
#      significance status in BOTH models, so DEGs gained and lost under the
#      depth model can be told apart from those that hold up.
#   3. plot_depth_concordance_sig.pdf -- parent-significant genes only.
#   4. source_data_depth_concordance*.tsv -- one row per drawn point.
#   5. stats_depth_adjustment_summary.tsv -- per celltype: correlations, Deming
#      and OLS slopes, coefficient movement, sign flips, DEG retention, plus the
#      depth context.
#   6. stats_depth_adjustment.txt -- human-readable log.
#
# No statistics are drawn on the figures; they live in the summary TSV and the log.
#
# Lightweight: reads small result TSVs only. Run directly (no PBS), e.g.
#   Rscript R/collate_depth_covariate.r
#   Rscript R/collate_depth_covariate.r --celltypes Exc-IT-L2-3-CBLN2-HOPX,Oligo

suppressPackageStartupMessages({
  library(argparse); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(purrr); library(stringr)
})

# Harden against dplyr generics being masked when this script is sourced into a
# session that already attached Seurat/AnnotationDbi/etc.
select    <- dplyr::select;    filter    <- dplyr::filter
rename    <- dplyr::rename;    mutate    <- dplyr::mutate
transmute <- dplyr::transmute; summarise <- dplyr::summarise
arrange   <- dplyr::arrange;   group_by  <- dplyr::group_by

## ---------------------------------------------------------------------------
## Arguments
## ---------------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--deg_dir", default = "deg",
  help = "Base DEG dir holding de_linear_distance/ and de_linear_distance_depth/ [default: %(default)s]")
parser$add_argument("--parent_tag", default = "de_linear_distance",
  help = "Subdir of --deg_dir holding the single-random-effect results [default: %(default)s]")
parser$add_argument("--depth_tag", default = "de_linear_distance_depth",
  help = "Subdir of --deg_dir holding the depth-adjusted results [default: %(default)s]")
parser$add_argument("--celltypes", default = "",
  help = "Comma-separated celltypes. Default: every celltype present in BOTH dirs.")
parser$add_argument("--padj", type = "double", default = 0.1,
  help = "Significance threshold, matching the project DEG default [default: %(default)s]")
parser$add_argument("--output_dir", default = "plots/depth_adjustment",
  help = "Output dir [default: %(default)s]")
parser$add_argument("--ncol", type = "integer", default = 4,
  help = "Facet columns [default: %(default)s]")
# Canvas overrides. The default size is derived from the facet grid
# (4.2 cm per column + margin), which is right for the full 13-celltype sheet but
# too narrow for a single-celltype panel: `ncol` is capped at the number of
# facets, so --ncol cannot be used to widen a one-facet figure, and the x-axis
# label and legend then run off the edge. These let a standalone panel be sized
# properly without touching the multi-facet defaults.
parser$add_argument("--width_cm", type = "double", default = 0,
  help = "Override figure width in cm (0 = derive from the facet grid) [default: %(default)s]")
parser$add_argument("--height_cm", type = "double", default = 0,
  help = "Override figure height in cm (0 = derive from the facet grid) [default: %(default)s]")
args <- parser$parse_args()

# Locate palettes.R whether run via Rscript (--file=) or interactively (cwd = root).
.get_script_dir <- function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1])) else "R"
}
.pal_candidates <- c(file.path(.get_script_dir(), "palettes.R"), "R/palettes.R", "palettes.R")
.pal <- .pal_candidates[file.exists(.pal_candidates)][1]
if (is.na(.pal)) stop("Cannot locate palettes.R (looked in: ",
                      paste(.pal_candidates, collapse = ", "), ")")
source(.pal)  # fig_theme

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

parent_dir <- file.path(args$deg_dir, args$parent_tag)
depth_dir    <- file.path(args$deg_dir, args$depth_tag)
if (!dir.exists(parent_dir)) stop("Parent results dir not found: ", parent_dir)
if (!dir.exists(depth_dir))    stop("Depth results dir not found: ", depth_dir,
                                  "\nRun deg_dream_linear_distance_depth_phf1.r first.")

res_file <- function(d, ct) file.path(d, ct, paste0(ct, "_dist_to_phf1_um_scaled.tsv"))

if (nzchar(args$celltypes)) {
  cts <- trimws(strsplit(args$celltypes, ",")[[1]])
} else {
  cts <- sort(intersect(basename(list.dirs(parent_dir, recursive = FALSE)),
                        basename(list.dirs(depth_dir,    recursive = FALSE))))
}

# ATTRITION REPORT. A celltype that aborted on the DEG script's
# MIN_CELLS_PER_DONOR = 20 donor-sufficiency guard has no result TSV in the PARENT
# arm either, which is unrelated to the covariate under test.
# A celltype missing from the parent is expected; one missing only from the depth
# arm means the depth run failed and is worth knowing about, so the two cases are
# reported separately here rather than folded into one silent intersection.
cts_all    <- sort(union(basename(list.dirs(parent_dir, recursive = FALSE)),
                         basename(list.dirs(depth_dir,  recursive = FALSE))))
has_parent <- file.exists(res_file(parent_dir, cts_all))
has_depth  <- file.exists(res_file(depth_dir,  cts_all))
attrition  <- data.frame(celltype = cts_all, parent = has_parent, depth = has_depth,
                         stringsAsFactors = FALSE)
miss_both  <- cts_all[!has_parent & !has_depth]
miss_par   <- cts_all[!has_parent &  has_depth]
miss_dep   <- cts_all[ has_parent & !has_depth]

cts <- cts[file.exists(res_file(parent_dir, cts)) & file.exists(res_file(depth_dir, cts))]
if (length(cts) == 0)
  stop("No celltype has a result TSV in BOTH ", parent_dir, " and ", depth_dir, ".")

cat("Celltypes with both models:", length(cts), "of", length(cts_all), "\n")
cat(" ", paste(cts, collapse = ", "), "\n")
if (length(miss_both))
  cat("  Not fitted in EITHER arm (expected - donor-sufficiency guard):\n    ",
      paste(miss_both, collapse = ", "), "\n", sep = "")
if (length(miss_dep))
  cat("  *** In the parent but NOT the depth arm - the depth run FAILED for these,\n",
      "      check R/logs/de_linear_distance_depth/:\n        ",
      paste(miss_dep, collapse = ", "), "\n", sep = "")
if (length(miss_par))
  cat("  In the depth arm but not the parent (no comparator):\n    ",
      paste(miss_par, collapse = ", "), "\n", sep = "")
cat("\n")

## ---------------------------------------------------------------------------
## Helpers carried over from collate_exclusion_radius.r (same definitions, so the
## two sets of numbers are computed identically and are directly comparable)
## ---------------------------------------------------------------------------
deming_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
  if (!is.finite(sxy) || sxy == 0) return(NA_real_)
  (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
}

cor_stats <- function(d) {
  d <- d %>% filter(is.finite(x), is.finite(y))
  if (nrow(d) < 3)
    return(tibble(n_pairs = nrow(d), pearson_r = NA_real_, pearson_lcl = NA_real_,
                  pearson_ucl = NA_real_, pearson_p = NA_real_, spearman_rho = NA_real_,
                  spearman_p = NA_real_, kendall_tau = NA_real_, kendall_p = NA_real_,
                  deming_slope = NA_real_, ols_slope = NA_real_,
                  median_abs_delta = NA_real_, median_pct_change = NA_real_,
                  max_abs_delta = NA_real_))
  pe <- suppressWarnings(cor.test(d$x, d$y, method = "pearson"))
  sp <- suppressWarnings(cor.test(d$x, d$y, method = "spearman"))
  kd <- suppressWarnings(cor.test(d$x, d$y, method = "kendall"))
  delta <- d$y - d$x
  pct   <- 100 * delta / abs(d$x)
  tibble(n_pairs      = nrow(d),
         pearson_r    = unname(pe$estimate),
         pearson_lcl  = unname(pe$conf.int[1]),
         pearson_ucl  = unname(pe$conf.int[2]),
         pearson_p    = pe$p.value,
         spearman_rho = unname(sp$estimate),
         spearman_p   = sp$p.value,
         kendall_tau  = unname(kd$estimate),
         kendall_p    = kd$p.value,
         deming_slope = deming_slope(d$x, d$y),
         # OLS for contrast: attenuated relative to Deming by the error in x, so a
         # gap between the two is expected and is not evidence of anything.
         ols_slope    = unname(coef(stats::lm(y ~ x, data = d))[2]),
         median_abs_delta  = median(abs(delta), na.rm = TRUE),
         median_pct_change = median(pct[is.finite(pct)], na.rm = TRUE),
         max_abs_delta     = max(abs(delta), na.rm = TRUE))
}

# Colours: Okabe-Ito (CB-safe) per the project palette rule, and the SAME mapping
# as collate_exclusion_radius.r so the depth figure reads like the radius figures.
SIG_CLASS_COLOURS <- c(
  "Significant in both"      = "#D55E00",
  "Parent only (lost)"       = "#0072B2",
  "Depth model only (gained)" = "#009E73",
  "Not significant"          = "grey80"
)
classify_sig <- function(padj_p, padj_d, thr) {
  a <- !is.na(padj_p) & padj_p < thr
  b <- !is.na(padj_d) & padj_d < thr
  factor(ifelse(a & b, "Significant in both",
         ifelse(a & !b, "Parent only (lost)",
         ifelse(!a & b, "Depth model only (gained)", "Not significant"))),
         levels = names(SIG_CLASS_COLOURS))
}
# Two-level split on the PARENT model only. Says nothing about the depth model, so
# the figure shows the fate of the parent-model DEGs without the eye being
# drawn to threshold crossings. Used for the "no sign changes" figure.
PARENT_CLASS_COLOURS <- c(
  "Significant in parent model"     = "#D55E00",
  "Not significant in parent model" = "grey80"
)
classify_parent <- function(padj_p, thr) {
  factor(ifelse(!is.na(padj_p) & padj_p < thr,
                "Significant in parent model", "Not significant in parent model"),
         levels = names(PARENT_CLASS_COLOURS))
}

## ---------------------------------------------------------------------------
## Read and pair
## ---------------------------------------------------------------------------
rd <- function(f) read.delim(f, stringsAsFactors = FALSE)
# dist_sd may be absent from the parent table, so recover it as
# logFC / logFC_per_unit. Both models are fitted on the same cells
# here, so the two values must agree; a mismatch means the cell sets differ and the
# pairing is invalid.
recover_dist_sd <- function(d) {
  r <- d$logFC / d$logFC_per_unit
  r <- r[is.finite(r)]
  if (!length(r)) return(NA_real_)
  stats::median(r)
}

pair_list <- list(); integ <- list()
for (ct in cts) {
  p <- rd(res_file(parent_dir, ct))
  f <- rd(res_file(depth_dir,    ct))
  sd_p <- recover_dist_sd(p); sd_f <- recover_dist_sd(f)
  j <- f %>% select(gene, y = logFC_per_unit, padj_d = padj, logFC_d = logFC) %>%
    inner_join(p %>% select(gene, x = logFC_per_unit, padj_p = padj, logFC_p = logFC),
               by = "gene") %>%
    mutate(celltype = ct,
           sig_class    = classify_sig(padj_p, padj_d, args$padj),
           parent_class = classify_parent(padj_p, args$padj))
  pair_list[[ct]] <- j
  integ[[ct]] <- tibble(
    celltype = ct,
    n_genes_parent = nrow(p), n_genes_depth = nrow(f), n_genes_paired = nrow(j),
    n_parent_only = length(setdiff(p$gene, f$gene)),
    n_depth_only    = length(setdiff(f$gene, p$gene)),
    dist_sd_parent = sd_p, dist_sd_depth = sd_f,
    dist_sd_ratio  = sd_f / sd_p,
    dist_sd_agrees = isTRUE(abs(sd_f / sd_p - 1) < 1e-6))
}
pair_df <- bind_rows(pair_list)
integ_df <- bind_rows(integ)

# Facet order: by descending parent DEG count, so the celltypes that carry the
# claim come first and the empty ones fall to the end.
ord <- pair_df %>% group_by(celltype) %>%
  summarise(n_sig = sum(padj_p < args$padj, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(n_sig), celltype)
pair_df$celltype <- factor(pair_df$celltype, levels = ord$celltype)

cat("=== Integrity: are the two models on the same cells and genes? ===\n")
print(as.data.frame(integ_df), row.names = FALSE)
if (any(!integ_df$dist_sd_agrees, na.rm = TRUE)) {
  cat("\n*** WARNING: recovered dist_sd DIFFERS between parent and depth model for: ",
      paste(integ_df$celltype[!integ_df$dist_sd_agrees], collapse = ", "), "\n")
  cat("*** The two fits are then NOT on the same cell set and logFC is not directly\n")
  cat("*** comparable. logFC_per_unit (used throughout here) still is.\n")
}
cat("\n")

## ---------------------------------------------------------------------------
## The figure
## ---------------------------------------------------------------------------
# Per-celltype square limits. A single shared scale across celltypes would squash
# most panels, because the coefficient range differs by an order of magnitude
# between celltypes; but x and y must share ONE range WITHIN each panel or the
# dashed identity line is not a true 45 degrees. So: free scales, plus a blank
# layer pinning both axes of each panel to the same padded range, plus
# aspect.ratio = 1. (theme(aspect.ratio) rather than coord_fixed, which does not
# combine cleanly with scales = "free".)
square_limits <- function(d) {
  d %>% group_by(celltype) %>%
    summarise(lo = min(c(x, y), na.rm = TRUE), hi = max(c(x, y), na.rm = TRUE),
              .groups = "drop") %>%
    mutate(pad = ifelse(hi - lo > 0, (hi - lo) * 0.04, 0.01),
           lo = lo - pad, hi = hi + pad) %>%
    select(celltype, lo, hi) %>%
    pivot_longer(c(lo, hi), values_to = "v") %>%
    transmute(celltype, x = v, y = v)
}

make_concordance <- function(d, tag, subset_label, colour_mode = c("none", "class", "parentsig")) {
  colour_mode <- match.arg(colour_mode)
  if (nrow(d) == 0) {
    warning("No genes in the '", subset_label, "' subset -- skipping that plot.",
            call. = FALSE)
    return(cor_stats(tibble(x = numeric(0), y = numeric(0)))[0, ] %>%
             mutate(celltype = character(0), gene_subset = character(0)))
  }
  ann <- d %>% group_by(celltype) %>% group_modify(~ cor_stats(.x)) %>% ungroup()
  n_ct <- length(unique(d$celltype))
  ncol <- max(1, min(args$ncol, n_ct)); nrow_ <- ceiling(n_ct / ncol)

  p <- ggplot(d, aes(x, y)) +
    geom_blank(data = square_limits(d)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey35")
  # Grey background genes first, so the coloured ones sit on top.
  if (colour_mode == "class") {
    p <- p +
      geom_point(data = d %>% filter(sig_class == "Not significant"),
                 aes(colour = sig_class), size = 0.4, alpha = 0.35, stroke = 0) +
      geom_point(data = d %>% filter(sig_class != "Not significant"),
                 aes(colour = sig_class), size = 0.6, alpha = 0.8, stroke = 0) +
      scale_colour_manual(values = SIG_CLASS_COLOURS, name = NULL, drop = TRUE)
  } else if (colour_mode == "parentsig") {
    p <- p +
      geom_point(data = d %>% filter(parent_class == "Not significant in parent model"),
                 aes(colour = parent_class), size = 0.4, alpha = 0.35, stroke = 0) +
      geom_point(data = d %>% filter(parent_class == "Significant in parent model"),
                 aes(colour = parent_class), size = 0.6, alpha = 0.8, stroke = 0) +
      scale_colour_manual(values = PARENT_CLASS_COLOURS, name = NULL, drop = TRUE)
  } else {
    p <- p + geom_point(size = 0.5, alpha = 0.5, colour = "grey25", stroke = 0)
  }
  p <- p +
    facet_wrap(~ celltype, ncol = ncol, scales = "free") +
    # (1|sample_id) is in BOTH models, so naming it on the x axis distinguishes
    # nothing; and "(+ depth_s)" restates "depth-adjusted". The only difference
    # between the two axes is the depth term, so that is all the labels say.
    labs(x = "logFC per unit log(distance), parent model",
         y = "logFC per unit log(distance),\ndepth-adjusted model") +
    theme_classic(base_size = 8) + fig_theme +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))
  if (colour_mode != "none")
    p <- p + theme(legend.position = "bottom",
                   legend.box.spacing = grid::unit(2, "pt")) +
      guides(colour = guide_legend(nrow = if (colour_mode == "class") 2 else 1,
                                   override.aes = list(size = 1.3, alpha = 1)))

  w <- if (args$width_cm  > 0) args$width_cm  else min(17, 4.2 * ncol + 1.5)
  h <- if (args$height_cm > 0) args$height_cm else
         min(23, 4.2 * nrow_ + switch(colour_mode, none = 1.0, class = 2.4, parentsig = 1.8))
  ggsave(file.path(args$output_dir, paste0("plot_depth_concordance_", tag, ".pdf")), p,
         width = w, height = h,
         units = "cm", device = "pdf", limitsize = FALSE)

  write.table(d %>% transmute(gene, celltype,
                              logFC_per_unit_parent = x, logFC_per_unit_depth = y,
                              padj_parent = padj_p, padj_depth = padj_d,
                              sig_class, parent_class,
                              sign_flip = sign(x) != sign(y)),
              file.path(args$output_dir,
                        paste0("source_data_depth_concordance_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ann %>% mutate(gene_subset = subset_label)
}

sig_parent <- pair_df$gene[pair_df$padj_p < args$padj]
ann_sig <- make_concordance(pair_df %>% filter(padj_p < args$padj), "sig",
                            sprintf("padj < %s in parent model", args$padj))
ann_all <- make_concordance(pair_df, "allgenes", "all genes")
invisible(make_concordance(pair_df, "allgenes_bysig",
                           "all genes, coloured by significance in both", "class"))
invisible(make_concordance(pair_df, "allgenes_parentsig",
                           "all genes, coloured by parent significance", "parentsig"))

## ---------------------------------------------------------------------------
## Per-celltype summary
## ---------------------------------------------------------------------------
retention <- pair_df %>% group_by(celltype) %>%
  summarise(
    n_tested            = dplyr::n(),
    n_sig_parent        = sum(padj_p < args$padj, na.rm = TRUE),
    n_sig_depth           = sum(padj_d < args$padj, na.rm = TRUE),
    n_sig_both          = sum(padj_p < args$padj & padj_d < args$padj, na.rm = TRUE),
    pct_parentsig_retained = ifelse(n_sig_parent > 0,
                              round(100 * n_sig_both / n_sig_parent, 2), NA_real_),
    jaccard_sig         = ifelse((n_sig_parent + n_sig_depth - n_sig_both) > 0,
                              round(n_sig_both / (n_sig_parent + n_sig_depth - n_sig_both), 4),
                              NA_real_),
    n_sign_flip_all     = sum(sign(x) != sign(y), na.rm = TRUE),
    n_sign_flip_parentsig = sum(padj_p < args$padj & sign(x) != sign(y), na.rm = TRUE),
    median_abs_logFC_parent = median(abs(x), na.rm = TRUE),
    median_abs_logFC_depth    = median(abs(y), na.rm = TRUE),
    attenuation_ratio   = round(median(abs(y), na.rm = TRUE) /
                                median(abs(x), na.rm = TRUE), 4),
    .groups = "drop")

# Depth context, if the depth run wrote it. This is what tells you whether an
# attenuated slope means depth absorbed something real or merely cost a degree of
# freedom. depth_s is a FIXED effect, so there is no ICC and no variance
# partition: the equivalent quantities are how much of the distance predictor
# depth could remove at all (depth_r2_dist) and whether depth explains any
# expression variance (n_sig_depth_0.1). laminar_rho_pooled is carried alongside.
read_depth_ctx <- function(ct) {
  f <- file.path(depth_dir, ct, paste0(ct, "_depth_summary.tsv"))
  if (!file.exists(f)) return(NULL)
  s <- rd(f)
  keep <- intersect(c("celltype", "depth_r_within", "depth_r2_within", "depth_rho_within",
                      "depth_rho_dist", "depth_r_dist", "depth_r2_dist",
                      "depth_rel_sd", "median_abs_t_depth", "median_abs_t_dist",
                      "n_sig_depth_0.1", "n_sig_depth_0.05", "pct_sig_depth_0.1",
                      "laminar_rho_pooled", "laminar_axis_warn"), colnames(s))
  s <- s[, keep, drop = FALSE]
  # n_sig_depth_0.1 would collide with the retention table's n_sig_depth prefix
  # on a partial-match read; rename to something unambiguous in the summary.
  names(s)[names(s) == "n_sig_depth_0.1"]  <- "n_genes_sig_for_depth"
  names(s)[names(s) == "n_sig_depth_0.05"] <- "n_genes_sig_for_depth_05"
  s
}
re_ctx <- bind_rows(lapply(cts, read_depth_ctx))

summ <- ann_all %>% select(-gene_subset) %>%
  left_join(retention, by = "celltype")
if (!is.null(re_ctx) && nrow(re_ctx) > 0)
  summ <- summ %>% left_join(re_ctx, by = "celltype")
summ <- summ %>% left_join(integ_df %>% select(celltype, dist_sd_parent, dist_sd_depth,
                                               dist_sd_agrees), by = "celltype") %>%
  arrange(desc(n_sig_parent))

write.table(summ, file.path(args$output_dir, "stats_depth_adjustment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ann_sig, file.path(args$output_dir, "stats_depth_adjustment_summary_sig.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---------------------------------------------------------------------------
## Human-readable log
## ---------------------------------------------------------------------------
sink(file.path(args$output_dir, "stats_depth_adjustment.txt"))
cat("Cortical-depth-adjusted DEG vs unadjusted parent\n")
cat("================================================\n\n")
cat("Parent model: ~ dist_to_phf1_um_scaled +           nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
cat("Depth model:  ~ dist_to_phf1_um_scaled + depth_s + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n\n")
cat("Parent results: ", parent_dir, "\n")
cat("Depth results:  ", depth_dir, "\n")
cat("padj threshold: ", args$padj, "\n")
cat("Celltypes:      ", length(cts), " (", paste(cts, collapse = ", "), ")\n\n", sep = "")
cat("Axis quantity: logFC_per_unit = logCPM per unit log(um), i.e. logFC / dist_sd.\n")
cat("Both models are fitted on the same cells so dist_sd should be identical; this is\n")
cat("verified below rather than assumed.\n\n")
cat("--- Celltype attrition ---\n")
cat("Which of the ", nrow(attrition), " celltype dirs produced a result TSV in each arm.\n",
    "A celltype absent from BOTH aborted on the DEG script's MIN_CELLS_PER_DONOR = 20\n",
    "donor-sufficiency guard and has nothing to do with depth. One absent from the DEPTH\n",
    "arm only means that run failed and the contrast is incomplete.\n\n", sep = "")
print(attrition, row.names = FALSE)
if (length(miss_dep))
  cat("\n*** INCOMPLETE: in the parent but not the depth arm: ",
      paste(miss_dep, collapse = ", "),
      "\n*** Check R/logs/de_linear_distance_depth/ before quoting this table.\n", sep = "")
cat("\n--- Integrity: same cells, same genes? ---\n")
print(as.data.frame(integ_df), row.names = FALSE)
cat("\n--- Per-celltype concordance and DEG retention ---\n")
print(as.data.frame(summ), row.names = FALSE)
cat("\n--- How to read this ---\n")
cat("Adjusting for a correlated covariate removes the shared component of the\n")
cat("predictor. **depth_r2_within** is how much of the distance predictor depth_s\n")
cat("could remove at all. Use the WITHIN-DONOR column, not the pooled depth_r2_dist:\n")
cat("both models carry (1|sample_id), so they condition on donor-demeaned\n")
cat("predictors, and the pooled figure mixes in between-donor differences the random\n")
cat("intercept has already absorbed. Measured here the pooled value understates the\n")
cat("within-donor R^2 by 2x to >1000x and for some celltypes has the OPPOSITE SIGN.\n")
cat("The pooled columns are retained for reference only.\n\n")
cat("An attenuated coefficient is EXPECTED in proportion to that overlap, and is not\n")
cat("by itself evidence that the gradient was a positional confound. Discriminate\n")
cat("using all three of:\n")
cat("  (a) deming_slope / attenuation_ratio - how far the coefficients moved;\n")
cat("  (b) depth_r2_within and n_genes_sig_for_depth - whether depth_s could remove\n")
cat("      much of the predictor, and whether it explains any expression variance at\n")
cat("      all. If depth explains nothing, the fit IS the parent model plus a wasted\n")
cat("      degree of freedom, and any change is noise;\n")
cat("  (c) n_sign_flip_parentsig - a genuine confound scrambles signs, whereas lost\n")
cat("      shared variance shrinks magnitudes toward zero WITHOUT flipping them.\n")
cat("\n--- The depth axis ---\n")
cat("depth_s is the standardised position along each sample's pia->white-matter axis.\n")
cat("laminar_rho_pooled is the Spearman correlation between the laminar excitatory\n")
cat("ordinal (L2-3 -> L3-5 -> L5 -> L6-IT -> L6-CT) and depth_rel.\n")
cat("Genes that change sign or lose significance under this adjustment are flagged per\n")
cat("gene in the source data (sign_flip, sig_class).\n")
if (!is.null(summ$laminar_axis_warn) && any(as.logical(summ$laminar_axis_warn), na.rm = TRUE)) {
  cat("\nlaminar_axis_warn is set (pooled laminar rho ",
      sprintf("%+.3f", suppressWarnings(mean(summ$laminar_rho_pooled, na.rm = TRUE))),
      ").\n",
      sep = "")
}
if (!is.null(summ$n_genes_sig_for_depth)) {
  inert <- summ$celltype[which(summ$n_genes_sig_for_depth == 0)]
  if (length(inert))
    cat("\n*** depth_s reached padj < ", args$padj, " for ZERO genes in: ",
        paste(inert, collapse = ", "),
        "\n*** For those celltypes depth_s explains no detectable expression variance, so\n",
        "*** the depth model is effectively the parent model plus one extra parameter.\n",
        sep = "")
}
cat("\n--- sessionInfo ---\n")
print(sessionInfo())
sink()

cat("=== Per-celltype summary ===\n")
print(as.data.frame(summ %>% select(any_of(c("celltype", "n_tested", "n_sig_parent",
  "n_sig_depth", "pct_parentsig_retained", "pearson_r", "deming_slope",
  "attenuation_ratio", "n_sign_flip_parentsig", "depth_r2_within",
  "n_genes_sig_for_depth")))), row.names = FALSE)
cat("\n[depth] wrote outputs to", args$output_dir, "\n")
