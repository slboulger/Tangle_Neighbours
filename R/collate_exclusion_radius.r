#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# collate_exclusion_radius.r
#
# Figure panels: S1A
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# collate_exclusion_radius.r
#
# Collates the exclusion-radius sensitivity sweep produced by
# deg_dream_linear_distance_phf1.r --min_distance_um (see run_deg_linear_distance.sh)
# for ONE celltype, and asks whether the tangle-distance coefficients survive
# progressively excluding the cells nearest a PHF1+ neuron.
#
# Reads, per radius R:
#   <deg_dir>/de_linear_distance/<ct>/<ct>_dist_to_phf1_um_scaled.tsv          (R = 0)
#   <deg_dir>/de_linear_distance_min<R>um/<ct>/<ct>_dist_to_phf1_um_scaled.tsv (R > 0)
#
# IMPORTANT — which coefficient is comparable. dist_sd is recomputed from the cells
# retained at each radius, so logFC (logCPM per SD of transformed distance) is NOT
# comparable across radii. Everything here therefore uses the back-divided
# logFC_per_unit / CI.L_per_unit / CI.R_per_unit columns, which are in logCPM per
# unit log(um) and are comparable. dist_sd per radius is reported in the stats log
# so the difference is visible rather than implicit.
#
# The background gene set is fixed pre-radius upstream, so the tested gene set
# should be identical at every radius. That is verified here and reported.
#
# Deliverables (plotting triple: PDF + source_data TSV + stats TXT):
#   1. source_data_exclusion_radius_<ct>.tsv          -- tidy gene x radius
#   2. plot_exclusion_radius_staircase_<ct>.pdf       -- coefficient staircase, facet/gene
#   3. plot_exclusion_radius_concordance_<ct>.pdf     -- radius 0 vs each radius (r0-sig)
#      plot_exclusion_radius_concordance_allgenes_<ct>.pdf  -- same, all genes
#      plot_exclusion_radius_concordance_allgenes_bysig_<ct>.pdf -- all genes, points
#        coloured by significance status (both / radius-0 only / this radius only /
#        neither), so the DEGs stand out against the full cloud and genes gained vs
#        lost at each radius are distinguishable.
#      plot_exclusion_radius_concordance_allgenes_r0sig_<ct>.pdf -- all genes, points
#        coloured ONLY by significance in the reference (radius-0) model. This is the
#        figure for "no radius-0 significant gene changes sign": no vermillion
#        point falls in an off-diagonal quadrant. Threshold crossings are ignored
#        here on purpose, so nothing distracts from the sign claim.
#      All three share one scale across x and y, identical breaks and coord_fixed(1),
#      so the dashed identity line is a true 45 degrees and the panels are directly
#      comparable. No statistics are drawn on the figures -- they are in the summary
#      TSV and the stats log.
#   4. stats_exclusion_radius_<ct>_summary.tsv        -- per-radius summary + correlations
#   5. stats_exclusion_radius_<ct>.txt                -- human-readable log
#
# Lightweight: reads small result TSVs only. Run directly (no PBS), e.g.
#   Rscript R/collate_exclusion_radius.r --deg_dir deg --celltype Exc-IT-L2-3-CBLN2-HOPX

suppressPackageStartupMessages({
  library(argparse); library(dplyr); library(tidyr); library(tibble)
  library(readr); library(ggplot2); library(purrr); library(stringr)
})

# Harden against dplyr generics being masked when this script is sourced into a
# session that already attached Seurat/AnnotationDbi/etc. (a bare select() then
# fails with "unused arguments" -- the classic symptom).
select    <- dplyr::select;    filter    <- dplyr::filter
rename    <- dplyr::rename;    mutate    <- dplyr::mutate
transmute <- dplyr::transmute; summarise <- dplyr::summarise
arrange   <- dplyr::arrange;   group_by  <- dplyr::group_by
slice_head <- dplyr::slice_head

# Locate palettes.R whether run via Rscript (--file=) or interactively (cwd = project root).
.get_script_dir <- function() {
  ca <- commandArgs(FALSE)
  f  <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
.pal_candidates <- c(file.path(.get_script_dir(), "palettes.R"), "R/palettes.R", "palettes.R")
.pal <- .pal_candidates[file.exists(.pal_candidates)][1]
if (is.na(.pal)) stop("Cannot locate palettes.R (looked in: ",
                      paste(.pal_candidates, collapse = ", "), ")")
source(.pal)  # fig_theme

## ---------------------------------------------------------------------------
## Arguments
## ---------------------------------------------------------------------------
parser <- ArgumentParser(
  description = paste("Collate the exclusion-radius sensitivity sweep for one celltype:",
                      "coefficient staircase and radius-vs-radius concordance."))
parser$add_argument("--deg_dir", default = "deg",
  help = "Base deg output dir holding the per-radius result dirs [default: deg]")
parser$add_argument("--celltype", required = TRUE,
  help = "Celltype, exactly as the SCE/result directory names it (may contain spaces)")
parser$add_argument("--radii", default = "0,10,20,30,50",
  help = "Comma-separated exclusion radii in um [default: 0,10,20,30,50]")
parser$add_argument("--genes", default = NULL,
  help = "Optional path to a one-column file of gene symbols for the staircase facets")
parser$add_argument("--output_dir", default = NULL,
  help = "Output dir [default: plots/exclusion_radius_<celltype>]")
parser$add_argument("--dir_prefix", default = "de_linear_distance",
  help = paste("Result dir prefix, so a variant sweep can be collated too",
               "(e.g. de_linear_distance_pct3p5) [default: de_linear_distance]"))
parser$add_argument("--log_dir", default = "R/logs/de_linear_distance",
  help = paste("Log dir used ONLY to recover n_cells for radius 0 when its table",
               "lacks the n_cells column [default: R/logs/de_linear_distance]"))
parser$add_argument("--coef_file_suffix", default = "_dist_to_phf1_um_scaled.tsv",
  help = "Per-celltype main results filename suffix [default: _dist_to_phf1_um_scaled.tsv]")
parser$add_argument("--staircase_top_n", default = 12L, type = "integer",
  help = "Genes in the staircase when --genes is not supplied [default: 12]")
parser$add_argument("--padj", default = 0.1, type = "double",
  help = "padj threshold defining 'significant at radius 0' [default: 0.1]")
args <- parser$parse_args()

ct      <- args$celltype
ct_slug <- gsub("[^A-Za-z0-9]", "_", ct)
out_dir <- if (is.null(args$output_dir)) {
  file.path("plots", paste0("exclusion_radius_", ct_slug))
} else {
  args$output_dir
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

radii <- as.numeric(trimws(strsplit(args$radii, ",")[[1]]))
radii <- sort(unique(radii[is.finite(radii)]))
if (length(radii) < 2) stop("--radii must name at least two radii to compare.")
if (!0 %in% radii) stop("--radii must include 0: it is the reference the sweep is compared against.")

# Mirrors the dir_tag logic in deg_dream_linear_distance_phf1.r exactly: radius 0
# keeps the canonical path, any other radius gets a min<R>um suffix (decimals as p).
radius_dir_tag <- function(r) {
  if (r <= 0) return(args$dir_prefix)
  paste0(args$dir_prefix,
         sprintf("_min%sum", gsub("\\.", "p", format(r, trim = TRUE, scientific = FALSE))))
}
fmt_r <- function(r) format(r, trim = TRUE, scientific = FALSE)

## ---------------------------------------------------------------------------
## Orthogonal (Deming) regression slope
## ---------------------------------------------------------------------------
# Both axes are estimated coefficients carrying comparable error, so an OLS slope
# is attenuated by the error in x and would understate agreement. This is the
# closed-form Deming / major-axis solution under EQUAL error variances in x and y
# -- a reasonable assumption here (same model, same gene, similar cell counts),
# and stated rather than hidden because the slope is only interpretable under it.
# No orthogonal-regression package (deming/lmodel2/smatr) exists in this project,
# so it is implemented inline rather than adding a dependency to dgenv.
deming_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  sxx <- var(x); syy <- var(y); sxy <- cov(x, y)
  if (!is.finite(sxy) || sxy == 0) return(NA_real_)
  (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy)
}

## ---------------------------------------------------------------------------
## Significance class of a gene at a radius, relative to radius 0
## ---------------------------------------------------------------------------
# Colours are Okabe-Ito (CB-safe) per the project palette rule: vermillion for the
# genes that hold up, blue for those lost, bluish-green for those gained. Defined
# before pair_df because that frame calls classify_sig().
SIG_CLASS_COLOURS <- c(
  "Significant at both"        = "#D55E00",
  "Radius 0 only (lost)"       = "#0072B2",
  "This radius only (gained)"  = "#009E73",
  "Not significant"            = "grey80"
)
classify_sig <- function(padj_r0, padj_r, thr) {
  a <- !is.na(padj_r0) & padj_r0 < thr
  b <- !is.na(padj_r)  & padj_r  < thr
  factor(ifelse(a & b, "Significant at both",
         ifelse(a & !b, "Radius 0 only (lost)",
         ifelse(!a & b, "This radius only (gained)", "Not significant"))),
         levels = names(SIG_CLASS_COLOURS))
}

# Simpler two-level split on the REFERENCE model only: was the gene significant at
# radius 0? Deliberately says nothing about significance at the exclusion radius,
# so the figure shows the fate of the radius-0 DEGs without the eye being
# drawn to threshold crossings. Used for the "no sign changes" figure.
R0_CLASS_COLOURS <- c(
  "Significant at radius 0"     = "#D55E00",
  "Not significant at radius 0" = "grey80"
)
classify_r0 <- function(padj_r0, thr) {
  factor(ifelse(!is.na(padj_r0) & padj_r0 < thr,
                "Significant at radius 0", "Not significant at radius 0"),
         levels = names(R0_CLASS_COLOURS))
}

## ---------------------------------------------------------------------------
## Read the per-radius result tables
## ---------------------------------------------------------------------------
# n_cells for radius 0: if the table lacks the n_cells column, recover it from
# "Cells modelled:" in the session summary of that run's log.
recover_n_cells_from_log <- function(celltype) {
  f <- file.path(args$log_dir, paste0(celltype, ".log"))
  if (!file.exists(f)) return(NA_integer_)
  txt <- readLines(f, warn = FALSE)
  ln <- grep("^Cells modelled:", txt, value = TRUE)
  if (length(ln))
    return(suppressWarnings(as.integer(trimws(sub("^Cells modelled:", "", tail(ln, 1))))))
  # Fallback: if the session summary is absent, use the post-cap count. At radius 0
  # nothing is dropped by the radius filter, so the post-cap count IS the modelled
  # count.
  ln <- grep("^Cells after excluding PHF1\\+", txt, value = TRUE)
  if (length(ln))
    return(suppressWarnings(as.integer(trimws(sub(".*:", "", tail(ln, 1))))))
  NA_integer_
}

read_one <- function(r) {
  f <- file.path(args$deg_dir, radius_dir_tag(r), ct, paste0(ct, args$coef_file_suffix))
  if (!file.exists(f)) {
    warning("No result table for radius ", fmt_r(r), " um at ", f,
            " -- dropping this radius. (Expected if the run aborted on donor",
            " sufficiency at this radius.)", call. = FALSE)
    return(NULL)
  }
  d <- suppressWarnings(read_tsv(f, show_col_types = FALSE, progress = FALSE))
  needed <- c("gene", "logFC", "logFC_per_unit", "CI.L_per_unit", "CI.R_per_unit",
              "pval", "padj", "distance_scale")
  miss <- setdiff(needed, names(d))
  if (length(miss)) stop("Result table ", f, " lacks required column(s): ",
                         paste(miss, collapse = ", "))
  # Provenance columns (n_cells, dist_sd) may be absent from a table; handle that.
  n_cells_col <- if ("n_cells" %in% names(d)) as.integer(d$n_cells[1]) else NA_integer_
  if (is.na(n_cells_col)) n_cells_col <- recover_n_cells_from_log(ct)
  d %>%
    transmute(
      gene, logFC, logFC_per_unit, CI.L_per_unit, CI.R_per_unit, pval, padj,
      distance_scale,
      min_distance_um = r,
      n_cells         = n_cells_col,
      dist_sd         = if ("dist_sd" %in% names(d)) dist_sd else NA_real_,
      file            = f
    )
}

cat("Celltype:", ct, "| radii:", paste(fmt_r(radii), collapse = ", "), "um\n")
res_list <- lapply(radii, read_one)
names(res_list) <- fmt_r(radii)
res_list <- res_list[!vapply(res_list, is.null, logical(1))]
if (!length(res_list)) stop("No result tables found under ", args$deg_dir, " for ", ct, ".")

dat <- bind_rows(res_list)
radii_found <- sort(unique(dat$min_distance_um))
if (!0 %in% radii_found)
  stop("The radius-0 reference table is missing; nothing to compare against.")
cat("Read", length(radii_found), "radii:", paste(fmt_r(radii_found), collapse = ", "), "um\n")

# Guard: per-unit coefficients are only comparable in a common unit.
scales_seen <- unique(dat$distance_scale)
if (length(scales_seen) > 1)
  stop("Mixed distance_scale across radii (", paste(scales_seen, collapse = ", "),
       "); logFC_per_unit are not in the same units and must not be compared.")

# Check the upstream pre-radius gene background actually held.
gene_sets  <- split(dat$gene, dat$min_distance_um)
ref_genes  <- gene_sets[[as.character(0)]]
bkg_report <- tibble(
  min_distance_um = as.numeric(names(gene_sets)),
  n_genes         = lengths(gene_sets),
  n_shared_with_r0 = vapply(gene_sets, function(g) length(intersect(g, ref_genes)), integer(1)),
  n_only_here      = vapply(gene_sets, function(g) length(setdiff(g, ref_genes)), integer(1)),
  n_only_r0        = vapply(gene_sets, function(g) length(setdiff(ref_genes, g)), integer(1))
) %>% arrange(min_distance_um)
bkg_identical <- all(bkg_report$n_only_here == 0 & bkg_report$n_only_r0 == 0)
if (!bkg_identical)
  warning("The tested gene set is NOT identical across radii -- coefficients are ",
          "not strictly comparable. See the gene-set check in the stats log.",
          call. = FALSE)

## ---------------------------------------------------------------------------
## (1) Tidy gene x radius source data
## ---------------------------------------------------------------------------
tidy_out <- dat %>%
  select(gene, min_distance_um, logFC_per_unit, CI.L_per_unit, CI.R_per_unit,
         pval, padj, n_cells) %>%
  arrange(gene, min_distance_um)
write.table(tidy_out, file.path(out_dir, paste0("source_data_exclusion_radius_", ct_slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

r0 <- dat %>% filter(min_distance_um == 0)
sig_r0 <- r0 %>% filter(!is.na(padj), padj < args$padj) %>% pull(gene)
cat("Genes significant at radius 0 (padj <", args$padj, "):", length(sig_r0), "\n")

## ---------------------------------------------------------------------------
## (2) Coefficient staircase
## ---------------------------------------------------------------------------
if (!is.null(args$genes)) {
  if (!file.exists(args$genes)) stop("--genes file not found: ", args$genes)
  gsel <- trimws(readLines(args$genes, warn = FALSE))
  gsel <- gsel[nzchar(gsel)]
  if (length(gsel) && tolower(gsel[1]) %in% c("gene", "genes", "symbol"))
    gsel <- gsel[-1]                      # tolerate a header line
  absent <- setdiff(gsel, unique(dat$gene))
  if (length(absent))
    warning("--genes: ", length(absent), " not in the results and dropped: ",
            paste(head(absent, 10), collapse = ", "), call. = FALSE)
  stair_genes <- intersect(gsel, unique(dat$gene))
  gene_source <- sprintf("--genes file (%s)", args$genes)
} else {
  stair_genes <- r0 %>% arrange(padj, pval) %>%
    slice_head(n = args$staircase_top_n) %>% pull(gene)
  gene_source <- sprintf("top %d by padj at radius 0", args$staircase_top_n)
}
if (!length(stair_genes)) stop("No genes selected for the staircase plot.")

stair_df <- dat %>%
  filter(gene %in% stair_genes) %>%
  mutate(gene = factor(gene, levels = stair_genes),
         sig  = !is.na(padj) & padj < args$padj)

p_stair <- ggplot(stair_df, aes(min_distance_um, logFC_per_unit)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey60") +
  geom_line(linewidth = 0.3, colour = "grey40") +
  geom_pointrange(aes(ymin = CI.L_per_unit, ymax = CI.R_per_unit, colour = sig),
                  linewidth = 0.3, size = 0.25) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey55"),
                      labels = c(`TRUE` = paste0("padj < ", args$padj), `FALSE` = "ns"),
                      name = NULL) +
  scale_x_continuous(breaks = radii_found) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  labs(x = expression("Exclusion radius (" * mu * "m)"),
       y = "logFC per unit log(distance)") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5),
        legend.position = "bottom",
        legend.box.spacing = grid::unit(2, "pt"),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt"))

n_row_stair <- ceiling(length(stair_genes) / 4)
ggsave(file.path(out_dir, paste0("plot_exclusion_radius_staircase_", ct_slug, ".pdf")),
       p_stair, width = 17, height = 2.6 * n_row_stair + 1.6, units = "cm", device = "pdf")
write.table(stair_df %>% select(gene, min_distance_um, logFC_per_unit,
                                CI.L_per_unit, CI.R_per_unit, pval, padj, sig),
            file.path(out_dir, paste0("source_data_exclusion_radius_staircase_", ct_slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---------------------------------------------------------------------------
## (3) Panel-wide concordance: radius 0 vs each non-zero radius
## ---------------------------------------------------------------------------
nz_radii <- setdiff(radii_found, 0)
if (!length(nz_radii)) stop("Only radius 0 was read; there is nothing to compare it with.")

pair_df <- dat %>%
  filter(min_distance_um %in% nz_radii) %>%
  select(gene, min_distance_um, y = logFC_per_unit, padj_r = padj) %>%
  inner_join(r0 %>% select(gene, x = logFC_per_unit, padj_r0 = padj), by = "gene") %>%
  mutate(radius_label = factor(paste0(fmt_r(min_distance_um), " um"),
                               levels = paste0(fmt_r(nz_radii), " um")),
         sig_class = classify_sig(padj_r0, padj_r, args$padj),
         r0_class  = classify_r0(padj_r0, args$padj))

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
  # Per-gene movement of the coefficient. The correlations say the ranking is
  # preserved; these say how far individual coefficients actually shifted, which
  # is the quantity a spillover artefact would inflate.
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

make_concordance <- function(d, tag, subset_label,
                             colour_mode = c("none", "class", "r0sig")) {
  colour_mode <- match.arg(colour_mode)
  # An empty subset (e.g. nothing significant at radius 0) must not abort the run:
  # return typed-but-empty stats so the summary join still works, and skip the plot.
  if (nrow(d) == 0) {
    warning("No genes in the '", subset_label, "' subset -- skipping that concordance ",
            "plot. The all-genes version is still written.", call. = FALSE)
    return(cor_stats(tibble(x = numeric(0), y = numeric(0)))[0, ] %>%
             mutate(min_distance_um = numeric(0), gene_subset = character(0)))
  }
  # Correlations are still computed (they feed the summary TSVs) but are NOT drawn
  # on the figure: the numbers live in stats_*_summary.tsv and the stats log.
  ann <- d %>% group_by(min_distance_um, radius_label) %>% group_modify(~ cor_stats(.x)) %>%
    ungroup()
  # One shared scale for both axes, with identical breaks, and coord_fixed so the
  # identity line is a true 45 degrees. Without this the two axes get independent
  # ranges and a point on y = x can sit visibly off the diagonal.
  lim <- range(c(d$x, d$y), na.rm = TRUE, finite = TRUE)
  pad <- if (diff(lim) > 0) diff(lim) * 0.04 else 0.01
  lim <- c(lim[1] - pad, lim[2] + pad)
  brk <- pretty(lim, n = 5)
  brk <- brk[brk >= lim[1] & brk <= lim[2]]
  p <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey35")
  # Draw the grey background genes first so the coloured ones sit on top of them.
  if (colour_mode == "class") {
    p <- p +
      geom_point(data = d %>% filter(sig_class == "Not significant"),
                 aes(colour = sig_class), size = 0.4, alpha = 0.35, stroke = 0) +
      geom_point(data = d %>% filter(sig_class != "Not significant"),
                 aes(colour = sig_class), size = 0.6, alpha = 0.8, stroke = 0) +
      scale_colour_manual(values = SIG_CLASS_COLOURS, name = NULL, drop = TRUE)
  } else if (colour_mode == "r0sig") {
    p <- p +
      geom_point(data = d %>% filter(r0_class == "Not significant at radius 0"),
                 aes(colour = r0_class), size = 0.4, alpha = 0.35, stroke = 0) +
      geom_point(data = d %>% filter(r0_class == "Significant at radius 0"),
                 aes(colour = r0_class), size = 0.6, alpha = 0.8, stroke = 0) +
      scale_colour_manual(values = R0_CLASS_COLOURS, name = NULL, drop = TRUE)
  } else {
    p <- p + geom_point(size = 0.5, alpha = 0.5, colour = "grey25", stroke = 0)
  }
  p <- p +
    facet_wrap(~ radius_label, nrow = 1) +
    scale_x_continuous(breaks = brk) +
    scale_y_continuous(breaks = brk) +
    coord_fixed(ratio = 1, xlim = lim, ylim = lim) +
    labs(x = "logFC per unit log(distance), radius 0",
         y = "logFC per unit log(distance),\nexclusion radius") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))
  if (colour_mode != "none")
    p <- p + theme(legend.position = "bottom",
                   legend.box.spacing = grid::unit(2, "pt")) +
      guides(colour = guide_legend(nrow = if (colour_mode == "class") 2 else 1,
                                   override.aes = list(size = 1.3, alpha = 1)))
  ggsave(file.path(out_dir, paste0("plot_exclusion_radius_concordance_", tag, ct_slug, ".pdf")),
         p, width = min(17, 4.2 * length(nz_radii) + 1.5),
         height = switch(colour_mode, none = 6.4, class = 7.8, r0sig = 7.4),
         units = "cm", device = "pdf")
  write.table(d %>% transmute(gene, min_distance_um, logFC_per_unit_r0 = x,
                              logFC_per_unit_radius = y, padj_r0, padj_radius = padj_r,
                              sig_class, r0_class,
                              sign_flip = sign(x) != sign(y)),
              file.path(out_dir, paste0("source_data_exclusion_radius_concordance_",
                                        tag, ct_slug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ann %>% mutate(gene_subset = subset_label)
}

ann_sig <- make_concordance(pair_df %>% filter(gene %in% sig_r0), "",
                            sprintf("padj < %s at radius 0", args$padj))
ann_all <- make_concordance(pair_df, "allgenes_", "all genes")
# Third variant: every tested gene, coloured by significance status at BOTH radii,
# so gains and losses can be told apart from the genes that hold up.
invisible(make_concordance(pair_df, "allgenes_bysig_", "all genes, coloured by significance",
                           colour_mode = "class"))
# Fourth variant: every tested gene, coloured ONLY by significance in the reference
# (radius-0) model. This is the figure for the claim "none of the radius-0
# significant genes changes sign": every vermillion point stays in the same
# x/y quadrant, and nothing about threshold crossings distracts from that.
invisible(make_concordance(pair_df, "allgenes_r0sig_",
                           "all genes, coloured by radius-0 significance",
                           colour_mode = "r0sig"))

## ---------------------------------------------------------------------------
## (4) Per-radius summary TSV
## ---------------------------------------------------------------------------
# Effect-size column: median |logFC_per_unit| is taken over the radius-0
# significant set at EVERY radius, so it tracks the same genes rather than a set
# that shrinks as significance is lost.
per_radius <- dat %>%
  group_by(min_distance_um) %>%
  summarise(
    n_cells        = first(n_cells),
    dist_sd        = first(dist_sd),
    n_genes_tested = dplyr::n(),
    n_sig_005      = sum(!is.na(padj) & padj < 0.05),
    n_sig_01       = sum(!is.na(padj) & padj < 0.1),
    median_abs_logFC_per_unit_r0sig =
      median(abs(logFC_per_unit[gene %in% sig_r0]), na.rm = TRUE),
    .groups = "drop"
  )

# Direction agreement, retained significance and DEG-set overlap vs radius 0.
# Sign flips are counted among the radius-0 significant genes: a spillover-driven
# gradient losing its source should attenuate or reverse, so a flip count near 0
# with retention high is the "no artefact" signature.
# Named by radius, and explicitly covering every radius read so a radius with
# zero significant genes still contributes a row rather than vanishing.
sig_at <- lapply(radii_found, function(r) {
  dat$gene[dat$min_distance_um == r & !is.na(dat$padj) & dat$padj < args$padj]
})
names(sig_at) <- as.character(radii_found)

concordance_tbl <- pair_df %>%
  filter(gene %in% sig_r0, is.finite(x), is.finite(y)) %>%
  group_by(min_distance_um) %>%
  summarise(
    n_r0sig            = dplyr::n(),
    n_same_sign        = sum(sign(x) == sign(y)),
    n_sign_flip        = sum(sign(x) != sign(y)),
    pct_sign_flip      = round(100 * sum(sign(x) != sign(y)) / dplyr::n(), 3),
    n_r0sig_still_sig  = sum(!is.na(padj_r) & padj_r < args$padj),
    pct_r0sig_retained = round(100 * sum(!is.na(padj_r) & padj_r < args$padj) / dplyr::n(), 2),
    .groups = "drop"
  )

# Gained/lost/Jaccard on the DEG sets themselves (not restricted to r0-significant).
setops_tbl <- tibble(min_distance_um = as.numeric(names(sig_at))) %>%
  mutate(
    n_sig_gained_vs_r0 = vapply(sig_at, function(g) length(setdiff(g, sig_r0)), integer(1)),
    n_sig_lost_vs_r0   = vapply(sig_at, function(g) length(setdiff(sig_r0, g)), integer(1)),
    jaccard_sig_vs_r0  = vapply(sig_at, function(g) {
      u <- length(union(g, sig_r0)); if (u == 0) NA_real_ else round(length(intersect(g, sig_r0)) / u, 4)
    }, numeric(1))
  )

cor_cols <- c("n_pairs", "pearson_r", "pearson_lcl", "pearson_ucl", "pearson_p",
              "spearman_rho", "spearman_p", "kendall_tau", "kendall_p",
              "deming_slope", "ols_slope", "median_abs_delta", "median_pct_change",
              "max_abs_delta")

build_summary <- function(ann) {
  per_radius %>%
    left_join(ann %>% select(min_distance_um, all_of(cor_cols)), by = "min_distance_um") %>%
    left_join(concordance_tbl, by = "min_distance_um") %>%
    left_join(setops_tbl, by = "min_distance_um") %>%
    mutate(celltype = ct,
           # Attenuation of the typical effect size relative to radius 0. A
           # spillover artefact carried by the excluded cells would push this < 1.
           median_abs_logFC_ratio_vs_r0 = round(
             median_abs_logFC_per_unit_r0sig /
               median_abs_logFC_per_unit_r0sig[min_distance_um == 0], 4)) %>%
    select(celltype, everything()) %>%
    arrange(min_distance_um)
}

summary_tbl <- build_summary(ann_sig)
write.table(summary_tbl,
            file.path(out_dir, paste0("stats_exclusion_radius_", ct_slug, "_summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_all <- build_summary(ann_all)
write.table(summary_all,
            file.path(out_dir, paste0("stats_exclusion_radius_", ct_slug, "_summary_allgenes.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---------------------------------------------------------------------------
## (5) Stats log
## ---------------------------------------------------------------------------
flips <- concordance_tbl

sink(file.path(out_dir, paste0("stats_exclusion_radius_", ct_slug, ".txt")))
cat("EXCLUSION-RADIUS SENSITIVITY -- linear tangle-distance DEG\n")
cat("==========================================================\n\n")
cat("Celltype:        ", ct, "\n")
cat("Date:            ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("deg_dir:         ", args$deg_dir, "\n")
cat("dir_prefix:      ", args$dir_prefix, "\n")
cat("Radii requested: ", paste(fmt_r(radii), collapse = ", "), "um\n")
cat("Radii read:      ", paste(fmt_r(radii_found), collapse = ", "), "um\n")
cat("padj threshold:  ", args$padj, "\n")
cat("distance_scale:  ", scales_seen, "\n")
cat("Staircase genes: ", gene_source, "\n\n")

cat("--- Source files ---\n")
for (f in unique(dat$file)) cat("  ", f, "\n")
cat("\n")

cat("--- Comparability of the coefficient ---\n")
cat("dist_sd is recomputed from the cells retained at each radius, so logFC (per SD\n")
cat("of transformed distance) is NOT comparable across radii. All numbers below use\n")
cat("logFC_per_unit = logFC / dist_sd (logCPM per unit log(um)), which is.\n")
cat("dist_sd per radius (NA = dist_sd column absent from the table):\n")
print(as.data.frame(per_radius %>% select(min_distance_um, n_cells, dist_sd)), row.names = FALSE)
cat("\n")

cat("--- Tested gene set across radii ---\n")
cat("The gene background is computed pre-radius upstream, so these should be identical.\n")
print(as.data.frame(bkg_report), row.names = FALSE)
cat("Identical across all radii read: ", bkg_identical, "\n")
if (!bkg_identical)
  cat("WARNING: gene sets differ -- the cross-radius comparison is confounded by a\n",
      "moving gene set. Check that the upstream background is computed on the full\n",
      "PHF1-negative set BEFORE the radius filter.\n", sep = "")
if (any(is.na(per_radius$n_cells)))
  cat("\nNOTE: n_cells unavailable for at least one radius (column absent from the table\n",
      "and no 'Cells modelled:' line was found under ", args$log_dir, ").\n", sep = "")
cat("\n")

cat("--- Per-radius summary (correlations vs radius 0, padj <", args$padj, "set) ---\n")
print(as.data.frame(summary_tbl), row.names = FALSE)
cat("\n--- Per-radius summary (correlations vs radius 0, ALL genes) ---\n")
print(as.data.frame(summary_all), row.names = FALSE)
cat("\n")

cat("--- Effect size: Pearson r with 95% CI, radius-0 significant genes ---\n")
cat("(r is the effect size for cross-radius agreement; the Deming slope is the\n")
cat(" orthogonal-regression slope, 1 = coefficients unchanged by the exclusion.)\n")
for (i in seq_len(nrow(ann_sig))) {
  a <- ann_sig[i, ]
  cat(sprintf("  %5s um: r = %.3f [%.3f, %.3f], p = %.3g | rho = %.3f, p = %.3g | Deming slope = %.3f | n = %d\n",
              fmt_r(a$min_distance_um), a$pearson_r, a$pearson_lcl, a$pearson_ucl,
              a$pearson_p, a$spearman_rho, a$spearman_p, a$deming_slope, a$n_pairs))
}
cat("\n")

cat("--- Sign flips and retained significance among radius-0 significant genes ---\n")
print(as.data.frame(flips), row.names = FALSE)
if (nrow(flips) && all(flips$n_sign_flip == 0))
  cat("\nNo gene reverses sign at any radius: the direction of every radius-0\n",
      "significant effect is preserved once the innermost cells are excluded.\n", sep = "")
cat("\n--- DEG set overlap vs radius 0 (padj <", args$padj, ") ---\n")
print(as.data.frame(setops_tbl %>% left_join(per_radius %>% select(min_distance_um, n_sig_01),
                                             by = "min_distance_um")), row.names = FALSE)
cat("\n--- Interpretation for a transcript-spillover artefact ---\n")
cat("Spillover from a PHF1+ neuron into neighbouring segmentation masks is a very\n")
cat("short-range effect, so excluding cells within 10-50 um should remove it. If the\n")
cat("distance gradient were driven by spillover, these coefficients would ATTENUATE\n")
cat("(Deming slope < 1, median_abs_logFC_ratio_vs_r0 < 1) and DEGs would be lost.\n")
cat("Stability across radii is therefore evidence AGAINST a spillover explanation.\n")
cat("\nRadius-0 significant genes (padj <", args$padj, "):", length(sig_r0), "\n")
cat("\n--- sessionInfo ---\n"); print(sessionInfo())
sink()

cat("\n[exclusion_radius] wrote outputs to", out_dir, "\n")
print(as.data.frame(summary_tbl %>% select(min_distance_um, n_cells, n_genes_tested,
                                           n_sig_005, n_sig_01, pearson_r, deming_slope)),
      row.names = FALSE)
