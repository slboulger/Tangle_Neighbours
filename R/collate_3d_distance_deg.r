#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# collate_3d_distance_deg.r
#
# Figure panels: S1C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# collate_3d_distance_deg.r
#
# Contrast of the out-of-plane (3D-corrected distance) DEG against the primary
# unadjusted model, in the same shape as plots/depth_adjustment_CBLN2/.
#
# The stats helpers (deming_slope, cor_stats, the significance classes and the
# Okabe-Ito mapping) are carried over VERBATIM from R/collate_depth_covariate.r,
# which took them from R/collate_exclusion_radius.r, so all three sets of
# numbers are computed identically and are directly comparable.
#
# ============================================================================
# THREE THINGS THAT DIFFER FROM THE DEPTH CONTRAST
# ============================================================================
#
# 1. dist_sd DELIBERATELY DISAGREES BETWEEN THE ARMS. The depth contrast fits
#    both models on the same cells, so it asserts dist_sd is identical. Here the
#    augmented distances are shorter and more cells fall inside the 1000 um cap
#    (24,807 against 24,672), so dist_sd moves (0.8167 against 0.7715). logFC is
#    per SD of log-distance and is therefore NOT comparable across the arms;
#    logFC_per_unit = logFC / dist_sd IS, and is the only axis used here. The
#    disagreement is reported, not asserted away.
#
# 2. THE VARIANT ARM IS R REPLICATE DRAWS, NOT ONE FIT. Simulated anchors are
#    drawn stochastically, so each config has R fits that are pooled with the
#    standard within/between decomposition (Qbar, T = Ubar + (1+1/R) B,
#    Barnard-Rubin df) before the contrast. Ubar comes from the MODERATED SE, and
#    the pooling treats the draw distribution as a modelling assumption (a
#    variance-inflation rule) rather than as a posterior predictive.
#
# 3. FOR THE IN-PLANE HIDDEN-ANCHOR ARMS, THE ATTENUATION IS READ AGAINST THE
#    DILUTION CONSTANT. Hidden anchors are a position-weighted random draw, not
#    an identification of which cells are mislabelled, so augmentation attenuates
#    the slope BY CONSTRUCTION. The measured constant is
#    slope_aug/slope_obs = 0.906 [0.896, 0.915] (R/test_tangle_3d_utils.R).
#    So this script reports
#
#        attenuation_ratio            what was observed
#        attenuation_vs_dilution      observed / 0.906 (in-plane arms only)
#
#    A value near 1 means the arm did nothing the PROCEDURE does not do on its
#    own; a value well below 1 indicates sensitivity to the anchor set.
#
# Usage:
#   Rscript R/collate_3d_distance_deg.r [--padj 0.1]

suppressPackageStartupMessages({
  library(argparse); library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})

.this_dir <- local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else "R"
})
source(file.path(.this_dir, "palettes.R"))

parser <- ArgumentParser()
parser$add_argument("--celltype",   default = "Exc-IT-L2-3-CBLN2-HOPX")
parser$add_argument("--parent_dir", default = "deg/de_linear_distance")
parser$add_argument("--deg_root",   default = "deg")
parser$add_argument("--manifest",   default = "dist3d/MANIFEST.tsv")
parser$add_argument("--dilution_tsv", default = "results/tangle_3d/armB_dilution_constant.tsv")
parser$add_argument("--padj",       type = "double", default = 0.1)
parser$add_argument("--output_dir", default = "plots/tangle_3d_adjustment_CBLN2")
args <- parser$parse_args()

CT  <- args$celltype
THR <- args$padj
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Helpers -- verbatim from collate_depth_covariate.r
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
         ols_slope    = unname(coef(stats::lm(y ~ x, data = d))[2]),
         median_abs_delta  = median(abs(delta), na.rm = TRUE),
         median_pct_change = median(pct[is.finite(pct)], na.rm = TRUE),
         max_abs_delta     = max(abs(delta), na.rm = TRUE))
}

SIG_CLASS_COLOURS <- c(
  "Significant in both"    = "#D55E00",
  "Parent only (lost)"     = "#0072B2",
  "3D model only (gained)" = "#009E73",
  "Not significant"        = "grey80"
)
classify_sig <- function(padj_p, padj_d, thr) {
  a <- !is.na(padj_p) & padj_p < thr
  b <- !is.na(padj_d) & padj_d < thr
  factor(ifelse(a & b, "Significant in both",
         ifelse(a & !b, "Parent only (lost)",
         ifelse(!a & b, "3D model only (gained)", "Not significant"))),
         levels = names(SIG_CLASS_COLOURS))
}
PARENT_CLASS_COLOURS <- c(
  "Significant in parent model"     = "#D55E00",
  "Not significant in parent model" = "grey80"
)
classify_parent <- function(padj_p, thr) {
  factor(ifelse(!is.na(padj_p) & padj_p < thr,
                "Significant in parent model", "Not significant in parent model"),
         levels = names(PARENT_CLASS_COLOURS))
}

rd <- function(f) read.delim(f, stringsAsFactors = FALSE)
res_file <- function(d, ct) file.path(d, ct, paste0(ct, "_dist_to_phf1_um_scaled.tsv"))

## ---------------------------------------------------------------------------
## Pool R replicate draws
##
## Qbar = mean(logFC_r); T = Ubar + (1 + 1/R) B; Barnard-Rubin df. Ubar is the
## mean WITHIN-draw variance, taken from the moderated SE recovered from the
## reported CI (CI.R - CI.L)/(2 * qt(0.975, df)) -- limma writes the CI at the
## moderated df, so this recovers exactly the SE that produced it.
##
## For R = 1 this reduces to the single fit and B = 0, which is what the
## identity arm needs.
## ---------------------------------------------------------------------------
pool_draws <- function(tabs) {
  R <- length(tabs)
  key <- Reduce(intersect, lapply(tabs, function(t) t$gene))
  g   <- lapply(tabs, function(t) t[match(key, t$gene), , drop = FALSE])

  # per-unit scale throughout: dist_sd differs between draws as well as arms
  Q <- vapply(g, function(t) t$logFC_per_unit, numeric(length(key)))
  # SE = logFC / t. EXACT, and it makes no assumption about the reference
  # distribution. Recovering it from the CI as (CI.R - CI.L)/(2*qnorm(0.975))
  # would be wrong: limma writes the interval at the MODERATED t df, not the
  # normal quantile, so that route inflates every SE.
  U <- vapply(g, function(t) {
    se_scaled <- abs(t$logFC / t$t)
    (se_scaled / t$dist_sd)^2
  }, numeric(length(key)))
  if (is.null(dim(Q))) { Q <- matrix(Q, nrow = 1); U <- matrix(U, nrow = 1) }

  # R = 1 is the identity arm: there is no between-draw variance to add and no
  # pooling to do, so pass the fitted result through untouched rather than
  # re-deriving p from a recovered SE. Re-deriving cannot improve on the
  # original and can only introduce discrepancies.
  if (R == 1L) {
    t1 <- g[[1]]
    return(tibble(gene = key, logFC_per_unit = t1$logFC_per_unit,
                  se_per_unit = abs(t1$logFC / t1$t) / t1$dist_sd,
                  t = t1$t, pval = t1$pval, padj = t1$padj, R = 1L,
                  frac_between = 0, dist_sd = t1$dist_sd[1]))
  }

  Qbar <- rowMeans(Q)
  Ubar <- rowMeans(U)
  B    <- if (R > 1) apply(Q, 1, var) else rep(0, length(key))
  Tt   <- Ubar + (1 + 1 / R) * B
  df   <- if (R > 1) (R - 1) * (1 + Ubar / ((1 + 1 / R) * B))^2 else Inf
  df[!is.finite(df)] <- 1e6
  tstat <- Qbar / sqrt(Tt)
  pval  <- 2 * stats::pt(-abs(tstat), df = df)

  tibble(gene = key,
         logFC_per_unit = Qbar,
         se_per_unit = sqrt(Tt),
         t = tstat, pval = pval,
         padj = stats::p.adjust(pval, "BH"),
         R = R,
         frac_between = ifelse(Tt > 0, (1 + 1 / R) * B / Tt, 0),
         dist_sd = mean(vapply(g, function(t) t$dist_sd[1], numeric(1))))
}

## ---------------------------------------------------------------------------
## Load
## ---------------------------------------------------------------------------
pf <- res_file(args$parent_dir, CT)
if (!file.exists(pf)) stop("Parent result not found: ", pf)
parent <- rd(pf)

man <- rd(args$manifest)
dil <- { tb <- rd(args$dilution_tsv); setNames(tb$value, tb$quantity) }
DILUTION <- unname(dil[["dilution_slope_aug_over_obs"]])

configs <- sub("^aug_", "", man$config)
found <- list(); missing <- character(0)
for (i in seq_along(configs)) {
  reps <- sprintf("rep_%02d", seq_len(man$R[i]))
  tabs <- list()
  for (rp in reps) {
    d <- file.path(args$deg_root, sprintf("de_linear_distance_3d_%s_%s", configs[i], rp))
    f <- res_file(d, CT)
    if (file.exists(f)) tabs[[rp]] <- rd(f) else missing <- c(missing, basename(d))
  }
  if (length(tabs)) found[[man$config[i]]] <- list(tabs = tabs, man = man[i, ])
}

if (!length(found))
  stop("No 3D DEG results found under ", args$deg_root, "/de_linear_distance_3d_*.\n",
       "Run the array first:  qsub R/run_deg_linear_distance_3d.sh\n",
       "(index 1 is the identity check -- confirm it passes before the rest.)")

## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------
# Formatting carried over VERBATIM from collate_depth_covariate.r's
# make_concordance(), so the 3D figure is visually interchangeable with the
# depth and exclusion-radius sensitivity figures: theme_classic + fig_theme,
# aspect.ratio 1, free square limits via geom_blank, dashed identity line,
# grey background genes drawn first, cm sizing.
square_limits <- function(d) {
  r <- range(c(d$x, d$y), na.rm = TRUE)
  pad <- ifelse(diff(r) > 0, diff(r) * 0.04, 0.01)
  data.frame(x = c(r[1] - pad, r[2] + pad), y = c(r[1] - pad, r[2] + pad))
}

make_concordance <- function(d, tag, colour_mode = c("none", "class", "parentsig")) {
  colour_mode <- match.arg(colour_mode)
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
    labs(x = "logFC per unit log(distance), parent model",
         y = "logFC per unit log(distance),\n3D-corrected model") +
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

  # 8.78 x 9.17 cm is the page size of the depth figure this one sits beside
  # (plots/depth_adjustment_CBLN2/plot_depth_concordance_allgenes_parentsig.pdf,
  # measured from its MediaBox). Its formula, 4.2*ncol + 1.5, assumes facets; at
  # one celltype that gives 5.7 cm, which clips both the x axis label and the
  # legend. Matching the rendered page keeps the two figures interchangeable.
  ggsave(file.path(args$output_dir, paste0("plot_tangle3d_concordance_", tag, ".pdf")),
         p, width = 8.78,
         height = 7.4 + switch(colour_mode, none = 0.0, class = 1.8, parentsig = 1.2),
         units = "cm", device = "pdf", limitsize = FALSE)

  write.table(d %>% transmute(gene, config, logFC_per_unit_parent = x,
                              logFC_per_unit_3d = y, padj_parent, padj_3d,
                              sig_class, parent_class, sign_flip),
              file.path(args$output_dir,
                        paste0("source_data_tangle3d_concordance_", tag, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(p)
}

## ---------------------------------------------------------------------------
## Contrast per config
## ---------------------------------------------------------------------------
all_d <- list(); summ <- list()
for (nm in names(found)) {
  pooled <- pool_draws(found[[nm]]$tabs)
  mi     <- found[[nm]]$man

  d <- parent %>%
    transmute(gene, x = logFC_per_unit, padj_parent = padj,
              dist_sd_parent = dist_sd) %>%
    inner_join(pooled %>% transmute(gene, y = logFC_per_unit, padj_3d = padj,
                                    dist_sd_3d = dist_sd, frac_between, R),
               by = "gene") %>%
    mutate(config = nm,
           sig_class    = classify_sig(padj_parent, padj_3d, THR),
           parent_class = classify_parent(padj_parent, THR),
           sign_flip    = sign(x) != sign(y))
  all_d[[nm]] <- d

  # PER-REPLICATE significance, before pooling. This is the diagnostic that
  # decides how to read a drop in n_sig_3d: if single draws reproduce the parent
  # count but the POOLED result does not, the loss is variance inflation from
  # draw-to-draw disagreement, NOT the distance correction weakening the field.
  n_sig_reps <- vapply(found[[nm]]$tabs,
                       function(t) sum(t$padj < THR, na.rm = TRUE), numeric(1))

  cs <- cor_stats(d)
  att <- median(abs(d$y), na.rm = TRUE) / median(abs(d$x), na.rm = TRUE)
  summ[[nm]] <- bind_cols(
    tibble(config = nm, tag = mi$tag, p_detect = mi$p_detect, weight = mi$weight,
           R_requested = mi$R, R_found = length(found[[nm]]$tabs)),
    cs,
    tibble(
      n_sig_parent = sum(d$padj_parent < THR, na.rm = TRUE),
      n_sig_3d     = sum(d$padj_3d < THR, na.rm = TRUE),
      n_sig_both   = sum(d$padj_parent < THR & d$padj_3d < THR, na.rm = TRUE),
      pct_parentsig_retained = round(100 * sum(d$padj_parent < THR & d$padj_3d < THR,
                                               na.rm = TRUE) /
                                     max(1, sum(d$padj_parent < THR, na.rm = TRUE)), 2),
      jaccard_sig = round(sum(d$padj_parent < THR & d$padj_3d < THR, na.rm = TRUE) /
                          max(1, sum(d$padj_parent < THR | d$padj_3d < THR, na.rm = TRUE)), 4),
      n_sign_flip_all       = sum(d$sign_flip, na.rm = TRUE),
      n_sign_flip_parentsig = sum(d$sign_flip & d$padj_parent < THR, na.rm = TRUE),
      median_abs_logFC_parent = median(abs(d$x), na.rm = TRUE),
      median_abs_logFC_3d     = median(abs(d$y), na.rm = TRUE),
      attenuation_ratio       = round(att, 4),
      # attenuation_vs_dilution (below): near 1 => the arm did nothing that
      # dilution does not do on its own.
      dilution_constant       = round(DILUTION, 4),
      # NA where no hidden anchors were added: with H = 0 there is no dilution
      # to correct for, and dividing by the constant would invent an apparent
      # 10% strengthening out of an arm that is by construction identical to
      # the parent.
      # NA unless the config is an IN-PLANE augmented-anchor arm. The 0.906
      # constant was calibrated on that procedure's own anchor-set dilution and
      # has no meaning for the unseen-plane arms, which add tangles outside the
      # section rather than relabelling cells inside it. It is also NA for the
      # identity arm, where nothing is added at all.
      attenuation_vs_dilution = if (grepl("^aug_p0", nm)) round(att / DILUTION, 4)
                                else NA_real_,
      median_frac_between     = round(median(d$frac_between, na.rm = TRUE), 4),
      n_sig_rep_median = median(n_sig_reps),
      n_sig_rep_min    = min(n_sig_reps),
      n_sig_rep_max    = max(n_sig_reps),
      dist_sd_parent = d$dist_sd_parent[1],
      dist_sd_3d     = d$dist_sd_3d[1],
      dist_sd_ratio  = round(d$dist_sd_3d[1] / d$dist_sd_parent[1], 4)))
}

summ <- bind_rows(summ)
write.table(summ, file.path(args$output_dir, "stats_tangle3d_adjustment_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Figures for the headline config.
hl <- summ$config[summ$tag == "headline" & summ$weight == "intensity"]
hl <- if (length(hl)) hl[1] else names(all_d)[1]
dh <- all_d[[hl]]
make_concordance(dh, "allgenes", "none")
make_concordance(dh, "allgenes_bysig", "class")
make_concordance(dh, "allgenes_parentsig", "parentsig")
make_concordance(dh %>% filter(padj_parent < THR), "sig", "class")

## ---------------------------------------------------------------------------
## Stats log
## ---------------------------------------------------------------------------
sp <- file.path(args$output_dir, "stats_tangle3d_adjustment.txt")
con <- file(sp, open = "wt"); sink(con, split = TRUE)

cat("Out-of-plane (3D) corrected-distance DEG vs the primary model\n")
cat("=============================================================\n\n")
cat("Primary model: ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
cat("3D arm:        the SAME model, with dist_to_phf1_um replaced by a 3D-corrected\n")
cat("               distance (one table per config; see dist3d/MANIFEST.tsv).\n\n")
cat("Parent results:", args$parent_dir, "\n")
cat("padj threshold:", THR, "| celltype:", CT, "\n")
if (length(missing))
  cat("\nMISSING result dirs (", length(missing), "):\n  ",
      paste(head(missing, 20), collapse = "\n  "), "\n", sep = "")

cat("\nAxis quantity: logFC_per_unit = logFC / dist_sd, i.e. logFC per unit log(distance).\n")
cat("dist_sd DIFFERS between the arms by construction -- the augmented distances\n")
cat("are shorter so more cells fall inside the 1000 um cap and the SD of\n")
cat("log(distance) moves. logFC (per SD) is therefore NOT comparable across arms;\n")
cat("logFC_per_unit is. Observed ratios:\n")
print(as.data.frame(summ[, c("config", "dist_sd_parent", "dist_sd_3d", "dist_sd_ratio")]),
      row.names = FALSE, digits = 5)

cat("\n--- Replicate pooling ---\n")
cat("Each config is R stochastic draws pooled as Qbar, T = Ubar + (1+1/R)B, with\n")
cat("Barnard-Rubin df. Ubar comes from the MODERATED SE; the pooling treats the\n")
cat("draw distribution as a modelling assumption (a variance-inflation rule), not\n")
cat("as a posterior predictive given the data.\n")
cat("median_frac_between is the share of total variance coming from draw-to-draw\n")
cat("disagreement; a large value means the result depends on which cells were drawn.\n")

cat("\n--- Summary ---\n")
print(as.data.frame(summ[, c("config", "tag", "p_detect", "weight", "R_found",
                             "pearson_r", "spearman_rho", "kendall_tau",
                             "deming_slope", "attenuation_ratio",
                             "attenuation_vs_dilution", "n_sig_parent", "n_sig_3d",
                             "pct_parentsig_retained", "n_sign_flip_parentsig",
                             "median_frac_between", "n_sig_rep_median",
                             "n_sig_rep_min", "n_sig_rep_max")]),
      row.names = FALSE, digits = 4)

cat("\n=== POOLED vs PER-REPLICATE SIGNIFICANCE ===\n")
cat("n_sig_3d is the POOLED count; n_sig_rep_median is a typical SINGLE draw.\n")
cat("Where a single draw reproduces the parent count and the pooled one does not,\n")
cat("the shortfall is VARIANCE INFLATION from draw-to-draw disagreement, not the\n")
cat("distance correction weakening the field. Cross-check it against median_frac_between\n")
cat("(the share of total variance that is between-draw) and against deming_slope\n")
cat("and n_sign_flip_parentsig, which say whether the coefficients moved at all.\n")
cat("A config with few replicates and a high frac_between is the least stable:\n")
cat("one outlying draw dominates B and the pooled count collapses.\n")

cat("\n=== attenuation_vs_dilution (in-plane hidden-anchor arms only) ===\n")
cat(sprintf("Arm B dilution constant: %.3f [%.3f, %.3f]\n", DILUTION,
            dil[["dilution_ci_lo"]], dil[["dilution_ci_hi"]]))
cat("Hidden anchors are a POSITION-WEIGHTED RANDOM DRAW, not an identification of\n")
cat("which cells are mislabelled, so the in-plane arm attenuates the distance\n")
cat("slope BY CONSTRUCTION.\n\n")
cat("attenuation_ratio is what was observed. attenuation_vs_dilution divides it by\n")
cat("the constant (NA for the unseen-plane and identity arms):\n")
cat("  ~1.0  the arm did nothing that the procedure does not do on its own.\n")
cat("  <<1.0 sensitive to the anchor set.\n")
cat("  >>1.0 the augmented distances strengthen the gradient.\n")

cat("\n--- Identity check ---\n")
idc <- summ[summ$tag == "identity", , drop = FALSE]
if (nrow(idc)) {
  cat(sprintf("max |dlogFC_per_unit| = %.3g | sign flips = %d | sig calls differing = %d\n",
              idc$max_abs_delta[1], idc$n_sign_flip_all[1],
              abs(idc$n_sig_parent[1] - idc$n_sig_3d[1])))
  if (idc$max_abs_delta[1] > 1e-8)
    cat("FAILED -- the fork is not a pure distance swap. Nothing else here is\n",
        "interpretable until this is resolved.\n")
  else
    cat("PASSED -- the fork changes nothing but the distance.\n")
} else {
  cat("NOT PRESENT. Run array index 1 first; it verifies that the fork is a pure\n")
  cat("distance swap.\n")
}

cat("\n"); print(sessionInfo())
sink(); close(con)

cat("\nWritten to", args$output_dir, "\n")
