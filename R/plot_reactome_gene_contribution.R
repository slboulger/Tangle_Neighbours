#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_reactome_gene_contribution.R
#
# Figure panels: S2B, S2C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_reactome_gene_contribution.R
#
# Which genes carry a module score's shift? Decomposes module scores into their per-gene
# contributions, for ONE celltype, for BOTH contrasts the main script reports:
#
#   distance  -- the slope against distance to the nearest PHF1+ neuron (PHF1-negative
#                cells, capped, canonical log-distance transform)
#   PHF1+/-   -- the PHF1+ minus PHF1- difference (ALL cells, no distance cap)
#
# --module selects what to decompose: 'significant' (the default -- every module the
# log-distance model called significant for this celltype, read from the main run's
# stats_reactome_lmm_coeffs.tsv), 'all', or an explicit comma-separated list.
# Default celltype: Exc-IT-L2-3-CBLN2-HOPX.
#
# GENES ARE FITTED ONCE, over the UNION of the selected modules. A gene's effect depends
# only on the celltype and the contrast, never on which module it is attributed to, so
# fitting per module would refit shared genes and could let one gene carry two different
# numbers in two tables. BH is still applied WITHIN each module, so a shared gene's
# ADJUSTED p legitimately differs between them -- the families differ, the estimate does not.
#
# WHY NOT JUST USE THE PER-GENE DEG TABLES. deg/de_linear_distance/ is a different
# pipeline: raw counts -> CPM (no TMM) -> voomWithDreamWeights -> dream, so its logFC is
# in logCPM per SD of log-distance, carrying voom precision weights. The module score
# lives in SCT data space and is fitted with an unweighted lmer. Coefficients from the two
# do not compose, so a "contribution" built from DEG logFCs would not sum to the slope
# actually drawn in the module-score figure. Everything here is therefore refitted in the
# module score's own space, with the SAME model as
# R/plot_reactome_stress_death_vs_phf1_distance.R.
#
# THE DECOMPOSITION. AddModuleScore is a difference of two means over SCT data,
#
#     score = mean(geneset genes) - mean(control genes)
#
# so, writing k for the number of scored geneset genes,
#
#     slope(score)         = slope(geneset_mean) - slope(control_mean)
#     slope(geneset_mean)  = (1/k) * SUM_g slope_g        contribution_g = slope_g / k
#
# The control set never has to be reverse-engineered out of Seurat's binning: it is
# recovered by subtraction, control_mean = geneset_mean - score, which is exact.
#
# The second identity needs the LMM slope of a mean to equal the mean of the per-gene LMM
# slopes. That is not free -- REML re-estimates the variance components for every response
# -- so it is CHECKED rather than assumed: the stats log reports the reconstruction
# residual for both identities. In simulation with heterogeneous per-gene slopes and noise
# (SD 0.2-1.2) the discrepancy was 0.07% of the slope.
#
# SCORES ARE BIT-IDENTICAL TO THE MAIN SCRIPT. AddModuleScore draws its control genes with
# a sequence of sample() calls, one per gene per set, so the control genes a set receives
# depend on how many sets preceded it. Scoring this module ALONE would therefore give a
# different score from the figure it is meant to explain. All modules are scored in one
# call here, exactly as the main script does, and the target module is then extracted.
# The module configuration is read out of the main script for the same reason -- so the
# two can never disagree about which pathways exist or what they are called.
#
# WHAT THIS IS. A descriptive decomposition of a fitted slope: the genes are the module's
# members, so each contribution describes how much of the module effect that member
# carries; rank by contribution and read the CI. The control-gene term is reported
# explicitly so the share of the effect carried by expression-matched control genes is
# visible alongside the geneset effect.
#
# OUTPUT (under --output_dir), per module <mod> x celltype <ct>. The distance files carry
# no contrast tag; the PHF1+/- ones are tagged "_phf1":
#   plot_genecontrib[_phf1]_<mod>_<ct>.pdf                -- per-gene effect + 95% CI, ranked
#   source_data_genecontrib[_phf1]_<mod>_<ct>.tsv         -- exact drawn rows + contributions
#   plot_genecontrib[_phf1]_decomposition_<mod>_<ct>.pdf  -- module / geneset / control effects
#   source_data_genecontrib[_phf1]_decomposition_<mod>_<ct>.tsv
# Written once for the whole run:
#   stats_genecontrib_summary_<ct>.tsv                    -- one row per module x contrast:
#       module/geneset/control effect, control share, n significant genes, top gene
#   stats_genecontrib_<ct>.txt                            -- ONE log: every module, both contrasts
#
# Use run_reactome_gene_contribution.sh

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(argparse)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")                    # fig_theme, forest_theme, neuron_order
source("R/reactome_module_genesets.R")    # REACTOME_GMT_DEFAULT
source("R/reactome_filtered_modules.R")   # reactome_filtered_modules(), CATALOGUE_RULES

set.seed(42)

# Sibling figure this one decomposes. Its module config comes from
# reactome_filtered_modules(); --module significant reads that script's coefficient table.
MAIN_SCRIPT <- "R/plot_reactome_stress_death_vs_phf1_distance.R"

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--seu", default = "seu_PHF1.rds", help = "Path to seu_PHF1.rds")
parser$add_argument("--gmt", default = REACTOME_GMT_DEFAULT, help = "Cached Enrichr Reactome GMT")
parser$add_argument("--celltype", default = "Exc-IT-L2-3-CBLN2-HOPX",
  help = "Celltype to decompose [default: Exc-IT-L2-3-CBLN2-HOPX]")
parser$add_argument("--module", default = "significant",
  help = paste("Which modules to decompose: 'significant' = every module whose",
               "log-distance model was significant for this celltype (read from",
               "--lmm_coeffs); 'all' = every scored module; or a comma-separated list of",
               "keys from the main script's REACTOME_PRIMARY [default: significant]"))
parser$add_argument("--lmm_coeffs",
  default = "plots/reactome_stress_death_vs_phf1_distance_1000um/stats_reactome_lmm_coeffs.tsv",
  help = paste("Main-script coefficient table, used to decide which modules are",
               "significant. Must come from a run at the SAME cap and seed."))
parser$add_argument("--output_dir", default = "plots/reactome_gene_contribution_1000um",
  help = "Output directory")
parser$add_argument("--seed", type = "integer", default = 42,
  help = "Seed; MUST match the main script's for the scores to agree [default: 42]")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "Restrict modelled cells to dist_to_phf1_um <= this (um); 0 = no cap [default: 1000]")
parser$add_argument("--top_n", type = "integer", default = 0,
  help = "Draw only the N largest |contribution| genes; 0 = draw all [default: 0]")
parser$add_argument("--overwrite", default = "yes", choices = c("yes", "no"),
  help = "Clear --output_dir before writing [default: yes]")
args <- parser$parse_args()

MAX_DIST <- if (args$max_dist_um > 0) args$max_dist_um else Inf

if (args$overwrite == "yes" && dir.exists(args$output_dir)) {
  old <- list.files(args$output_dir,
                    pattern = "^(plot|source_data|stats)_.*\\.(pdf|tsv|txt)$", full.names = TRUE)
  if (length(old)) { cat("Clearing", length(old), "existing file(s)\n"); unlink(old) }
}
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Constants -- MUST match the main script or the scores diverge.
MIN_GENES_PRESENT <- 10
CTRL              <- 100
NBIN              <- 24
SIG_ALPHA         <- 0.05

# Okabe-Ito up/down pairs: vermillion = tangle-associated, blue = the opposite.
# One pair per contrast, because the sign means different things:
#   distance -- the slope is per +1 SD of LOG DISTANCE, so NEGATIVE = higher near PHF1+
#   phf1     -- the estimate is PHF1+ minus PHF1-, so POSITIVE = higher in PHF1+
DIR_COLOURS      <- c(`higher near PHF1+` = "#D55E00", `higher far from PHF1+` = "#0072B2")
DIR_COLOURS_PHF1 <- c(`higher in PHF1+`   = "#D55E00", `lower in PHF1+`         = "#0072B2")

# The two contrasts this script decomposes.
TERMS <- list(
  dist = list(key = "dist", predictor = "dist_scaled", paired = FALSE,
              label   = "distance gradient",
              file    = "",                       # no suffix: the original file names
              xlab    = "log1p SCT expression per SD of log distance",
              declab  = "log1p SCT expression per SD of log distance",
              colours = DIR_COLOURS,
              dir_fun = function(est) ifelse(est < 0, names(DIR_COLOURS)[1], names(DIR_COLOURS)[2])),
  phf1 = list(key = "phf1", predictor = "phf1_pos", paired = TRUE,
              label   = "PHF1+ vs PHF1- contrast",
              file    = "_phf1",
              xlab    = "log1p SCT expression, PHF1+ minus PHF1-",
              declab  = "log1p SCT expression, PHF1+ minus PHF1-",
              colours = DIR_COLOURS_PHF1,
              dir_fun = function(est) ifelse(est > 0, names(DIR_COLOURS_PHF1)[1], names(DIR_COLOURS_PHF1)[2]))
)

##  ............................................................................
##  Helpers (copied from the main script unless noted)                      ####

pick_col <- function(meta, candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0)
    stop(sprintf("None of the %s columns (%s) found in seu meta.data",
                 what, paste(candidates, collapse = ", ")))
  hit[1]
}

prep_sets <- function(sets, panel, min_genes, collection, group) {
  present <- lapply(sets, function(g) intersect(g, panel))
  n_total <- vapply(sets, length, integer(1))
  n_pres  <- vapply(present, length, integer(1))
  scored  <- n_pres >= min_genes
  overlap <- tibble(
    collection = collection, group = group, set = names(sets),
    n_total = as.integer(n_total), n_present = as.integer(n_pres),
    frac_present = round(n_pres / pmax(n_total, 1L), 4), scored = scored)
  list(present = present[scored], overlap = overlap)
}

score_group <- function(seu_ct, present_sets, prefix, seed) {
  set.seed(seed)
  seu_ct <- tryCatch(
    AddModuleScore(seu_ct, features = present_sets, name = prefix,
                   assay = "SCT", ctrl = CTRL, nbin = NBIN, seed = seed),
    error = function(e) {
      cat("    AddModuleScore failed at nbin=", NBIN, " (", conditionMessage(e),
          ") - retrying at nbin=15\n", sep = "")
      set.seed(seed)
      AddModuleScore(seu_ct, features = present_sets, name = prefix,
                     assay = "SCT", ctrl = CTRL, nbin = 15, seed = seed)
    })
  cols <- paste0(prefix, seq_along(present_sets))
  stopifnot("AddModuleScore column count != number of scored sets" =
              length(cols) == length(present_sets) && all(cols %in% colnames(seu_ct@meta.data)))
  sc <- seu_ct@meta.data[, cols, drop = FALSE]
  colnames(sc) <- names(present_sets)
  sc
}

# One LMM, one response, one predictor of interest. `predictor` is "dist_scaled" for the
# distance gradient or "phf1_pos" for the PHF1+/- contrast; the covariates and the random
# effect are identical either way, and identical to the main script's two models. The
# coefficient is matched by prefix because a logical predictor becomes 'phf1_posTRUE'.
fit_term <- function(y, d, predictor) {
  sub <- d; sub$score <- y
  fml <- stats::as.formula(sprintf(
    "score ~ %s + nUMI_log + percent_neg + Sex + Age_s + PMI_s + (1 | sample_id)", predictor))
  m <- tryCatch(lmerTest::lmer(fml, data = sub, REML = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m))
  rn <- grep(paste0("^", predictor), rownames(cf), value = TRUE)
  if (!length(rn)) return(NULL)
  est <- cf[rn[1], "Estimate"]; se <- cf[rn[1], "Std. Error"]
  vc  <- as.data.frame(lme4::VarCorr(m))
  list(estimate = est, se = se, ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
       t = cf[rn[1], "t value"], df = cf[rn[1], "df"], p = cf[rn[1], "Pr(>|t|)"],
       sd_total = sqrt(sum(vc$vcov[is.na(vc$var2)])), singular = isSingular(m))
}

# Donor-paired difference for one gene (PHF1+ minus PHF1-). Returns dz and the count of
# donors moving with the mean.
donor_paired <- function(y, d) {
  pd <- data.frame(sample_id = d$sample_id, phf1_pos = d$phf1_pos, y = y) %>%
    group_by(sample_id) %>%
    summarise(np = sum(phf1_pos), nn = sum(!phf1_pos),
              diff = mean(y[phf1_pos]) - mean(y[!phf1_pos]), .groups = "drop") %>%
    filter(np >= 1, nn >= 1, is.finite(diff))
  if (!nrow(pd)) return(list(dz = NA_real_, n_donors = 0L, n_same_dir = NA_integer_))
  list(dz = if (nrow(pd) >= 2 && stats::sd(pd$diff) > 0)
              mean(pd$diff) / stats::sd(pd$diff) else NA_real_,
       n_donors = nrow(pd),
       n_same_dir = sum(sign(pd$diff) == sign(mean(pd$diff))))
}

# Fit every gene ONCE per contrast, over the UNION of the selected modules' genes.
#
# A gene's effect depends only on the celltype and the contrast, never on which module it
# is being attributed to, so fitting per module would refit shared genes and -- worse --
# could let the same gene carry two different numbers in two tables. Fitting the union
# once makes UBB's slope identical everywhere it appears, by construction.
fit_gene_cache <- function(d, E, genes, predictor, paired = FALSE) {
  cache <- list()
  for (i in seq_along(genes)) {
    g <- genes[i]
    r <- fit_term(E[, g], d, predictor)
    if (is.null(r)) { cat("  [warn] LMM failed for ", g, " (", predictor, ")\n", sep = ""); next }
    dp <- if (paired) donor_paired(E[, g], d)
          else list(dz = NA_real_, n_donors = NA_integer_, n_same_dir = NA_integer_)
    cache[[g]] <- tibble(
      gene = g, estimate = r$estimate, se = r$se, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
      t = r$t, df = r$df, p_LMM = r$p, cohens_d = r$estimate / r$sd_total,
      sd_total = r$sd_total, dz_donor_paired = dp$dz,
      n_donors_paired = dp$n_donors, n_donors_same_dir = dp$n_same_dir,
      singular = r$singular, mean_expr = mean(E[, g]))
    if (i %% 25 == 0) cat("    ", i, "/", length(genes), "\n", sep = "")
  }
  if (!length(cache)) stop("every per-gene LMM failed for predictor '", predictor, "'")
  cache
}

# Assemble ONE module's decomposition from the shared gene cache plus its own three
# aggregate fits. BH is applied WITHIN the module, so a gene shared between modules can
# have different adjusted p-values in each -- that is correct, the families differ.
decompose_module <- function(d, gene_cache, module_genes,
                             module_score, geneset_mean, control_mean, predictor) {
  K <- length(module_genes)
  agg <- list(module_score = fit_term(module_score, d, predictor),
              geneset_mean = fit_term(geneset_mean, d, predictor),
              control_mean = fit_term(control_mean, d, predictor))
  if (any(vapply(agg, is.null, logical(1))))
    stop("an aggregate LMM failed to fit for predictor '", predictor, "'")

  have <- module_genes[module_genes %in% names(gene_cache)]
  if (!length(have)) stop("no gene fits available for this module")
  gt <- bind_rows(gene_cache[have])
  gt$lmm_padj        <- p.adjust(gt$p_LMM, method = "BH")   # BH within THIS module
  gt$significant     <- gt$lmm_padj < SIG_ALPHA
  gt$contribution    <- gt$estimate / K                     # additive share of the geneset effect
  gt$contribution_lo <- gt$ci_lo / K
  gt$contribution_hi <- gt$ci_hi / K
  gt$pct_of_geneset  <- 100 * gt$contribution / agg$geneset_mean$estimate
  gt <- gt[order(gt$estimate), ]

  sum_contrib <- sum(gt$contribution)
  list(agg = agg, gt = gt, k = K, n_fitted = length(have), sum_contrib = sum_contrib,
       rec_geneset = sum_contrib - agg$geneset_mean$estimate,
       rec_module  = (agg$geneset_mean$estimate - agg$control_mean$estimate) -
                      agg$module_score$estimate)
}

# Which modules to decompose: 'significant' (from the main run's coefficient table),
# 'all', or an explicit comma-separated list.
select_modules <- function(spec, scored, ct, lmm_path) {
  if (identical(spec, "all")) return(scored)
  if (identical(spec, "significant")) {
    if (!file.exists(lmm_path))
      stop("--module significant needs the main script's coefficient table, not found:\n  ",
           lmm_path, "\n  Run plot_reactome_stress_death_vs_phf1_distance.R first, or pass",
           " an explicit module list.")
    tb <- utils::read.delim(lmm_path, stringsAsFactors = FALSE)
    need <- c("celltype", "module", "significant")
    if (!all(need %in% names(tb)))
      stop("coefficient table is missing ", paste(setdiff(need, names(tb)), collapse = ", "))
    sub <- tb[tb$celltype == ct, ]
    if (!nrow(sub))
      stop("no rows for celltype '", ct, "' in ", lmm_path,
           "\n  (celltypes present: ", paste(unique(tb$celltype), collapse = ", "), ")")
    sel <- sub$module[sub$significant %in% c(TRUE, "TRUE")]
    sel <- scored[scored %in% sel]          # keep REACTOME_PRIMARY order, drop unscored
    if (!length(sel))
      stop("no module was significant for ", ct, " in ", lmm_path,
           " -- nothing to decompose.")
    return(sel)
  }
  req <- trimws(strsplit(spec, ",")[[1]]); req <- req[nzchar(req)]
  bad <- setdiff(req, scored)
  if (length(bad))
    stop("not scored (or not a module key): ", paste(bad, collapse = ", "),
         "\n  Available: ", paste(scored, collapse = ", "))
  scored[scored %in% req]
}

##  ............................................................................
##  Configuration + gene sets                                               ####
# The module list is DERIVED from CATALOGUE_RULES by the same function the main
# script calls, so the two cannot disagree about the pathway list.
#
# Assigned AFTER the Seurat object loads, because the rules are defined relative to
# the CosMx panel and the panel is rownames(seu[["SCT"]]).

CT <- args$celltype

##  ............................................................................
##  Load, score (all modules, one call), extract the target                 ####
cat("\nLoading Seurat object:", args$seu, "\n")
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "SCT"
stopifnot(
  "celltype column missing"        = "celltype"        %in% colnames(seu@meta.data),
  "PHF1 column missing"            = "PHF1"            %in% colnames(seu@meta.data),
  "dist_to_phf1_um column missing" = "dist_to_phf1_um" %in% colnames(seu@meta.data),
  "percent.neg column missing"     = "percent.neg"     %in% colnames(seu@meta.data),
  "nCount_RNA column missing"      = "nCount_RNA"      %in% colnames(seu@meta.data),
  "sample_id column missing"       = "sample_id"       %in% colnames(seu@meta.data)
)
if (!CT %in% unique(as.character(seu$celltype)))
  stop("Celltype not present in the object: ", CT)
sex_col <- pick_col(seu@meta.data, c("Sex"), "Sex")
age_col <- pick_col(seu@meta.data, c("Age"), "Age")
pmi_col <- pick_col(seu@meta.data, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")

panel_genes <- rownames(seu[["SCT"]])

# Derive the modules from CATALOGUE_RULES -- same call, same rules, same genes as the
# main script. GENES COME FROM THE CATALOGUE, not from reactome_genesets() + intersect():
# the catalogue applies HGNC alias resolution (PARK2 -> PRKN, IL8 -> CXCL8, H2AFX ->
# H2AX, UFD1L -> UFD1), so a fresh intersect would drop PRKN from PINK1-PRKN Mediated
# Mitophagy.
cat("\nDeriving Reactome modules from CATALOGUE_RULES...\n")
mods <- reactome_filtered_modules(panel_genes, gmt_path = args$gmt)
REACTOME_PRIMARY <- mods$REACTOME_PRIMARY
MODULE_LABELS    <- mods$MODULE_LABELS
MODULE_GROUPS    <- mods$MODULE_GROUPS
rx <- list(genes = mods$genes)
cat("Modules:", length(REACTOME_PRIMARY), "\n")

prep <- prep_sets(rx$genes, panel_genes, MIN_GENES_PRESENT, "reactome", "reactome_primary")
SCORED <- names(REACTOME_PRIMARY)[names(REACTOME_PRIMARY) %in% names(prep$present)]

MODS <- select_modules(args$module, SCORED, CT, args$lmm_coeffs)
cat("Scored modules:", length(SCORED), "| decomposing", length(MODS), "\n")
for (m in MODS)
  cat(sprintf("  %-16s %-14s k = %3d  %s\n", m, unname(REACTOME_PRIMARY[m]),
              length(prep$present[[m]]), unname(MODULE_LABELS[m])))

GENE_UNION <- unique(unlist(prep$present[MODS], use.names = FALSE))
cat("Union of genes to fit:", length(GENE_UNION), "(vs",
    sum(lengths(prep$present[MODS])), "if fitted per module)\n")

seu_ct <- subset(seu, celltype == CT)
DefaultAssay(seu_ct) <- "SCT"
md <- seu_ct@meta.data
cat(CT, "cells:", nrow(md), "\n")

# All SCORED modules in ONE call, matching the main script's RNG sequence exactly --
# not just the ones being decomposed, or the control genes (and hence the scores) differ.
scores <- score_group(seu_ct, prep$present[SCORED], "RXM_", args$seed)

expr <- as.matrix(SeuratObject::GetAssayData(seu_ct, assay = "SCT", layer = "data")[GENE_UNION, , drop = FALSE])
# Per-module aggregates. control_mean is exact: AddModuleScore subtracts it, so
# control_mean = geneset_mean - score, with no need to replicate Seurat's binning.
AGG_VEC <- lapply(setNames(MODS, MODS), function(m) {
  gm <- colMeans(expr[prep$present[[m]], , drop = FALSE])
  list(geneset_mean = gm, module_score = scores[[m]], control_mean = gm - scores[[m]])
})

id_cands <- c("cell_ID", "cell_id", "CellID")
id_col   <- id_cands[id_cands %in% colnames(md)]
base <- data.frame(
  cell_id         = if (length(id_col)) as.character(md[[id_col[1]]]) else rownames(md),
  phf1_pos        = as.logical(as.character(md$PHF1)),
  dist_to_phf1_um = as.numeric(md$dist_to_phf1_um),
  percent_neg     = as.numeric(md$percent.neg),
  nUMI_log        = log2(as.numeric(md$nCount_RNA) + 1),
  sample_id       = as.character(md$sample_id),
  Sex             = as.character(md[[sex_col]]),
  Age             = as.numeric(md[[age_col]]),
  PMI             = as.numeric(md[[pmi_col]]),
  stringsAsFactors = FALSE)
expr_t <- t(expr)                                  # cells x genes, aligned with base
rm(seu, seu_ct); invisible(gc())


##  ............................................................................
##  Model frames: one per contrast                                          ####
# DISTANCE frame: PHF1-negative cells only, capped, canonical log-distance transform.
covar_ok <- !is.na(base$percent_neg) & !is.na(base$nUMI_log) &
            !is.na(base$Age) & !is.na(base$PMI) & !is.na(base$Sex)

keep_d <- covar_ok & !is.na(base$dist_to_phf1_um) & !base$phf1_pos
rows_d <- which(keep_d)          # row indices INTO `base`, so the per-module aggregate
                                 # vectors (which span all cells) can be subset to match
d_dist <- base[keep_d, ]; E_dist <- expr_t[keep_d, , drop = FALSE]
d_dist$Age_s <- as.numeric(scale(d_dist$Age)); d_dist$PMI_s <- as.numeric(scale(d_dist$PMI))
if (is.finite(MAX_DIST)) {
  sel <- d_dist$dist_to_phf1_um <= MAX_DIST
  d_dist <- d_dist[sel, ]; E_dist <- E_dist[sel, , drop = FALSE]; rows_d <- rows_d[sel]
}
if (any(!is.finite(d_dist$dist_to_phf1_um) | d_dist$dist_to_phf1_um <= 0))
  stop("non-positive or non-finite dist_to_phf1_um among modelled cells")
# Canonical transform: natural log, scaled by SD only, NOT centred.
dist_sd <- stats::sd(log(d_dist$dist_to_phf1_um))
d_dist$dist_scaled <- log(d_dist$dist_to_phf1_um) / dist_sd
d_dist$Sex <- droplevels(factor(d_dist$Sex)); d_dist$sample_id <- factor(d_dist$sample_id)

# PHF1+/- frame: ALL cells, NO distance cap and NO distance filter. PHF1+ cells have no
# dist_to_phf1_um at all, so filtering on distance would delete the very group being
# contrasted; and capping would select the PHF1-negative comparison group on proximity to
# PHF1+ cells, biasing the contrast toward the null. Matches the main script exactly.
rows_p <- which(covar_ok)
d_phf1 <- base[covar_ok, ]; E_phf1 <- expr_t[covar_ok, , drop = FALSE]
d_phf1$Age_s <- as.numeric(scale(d_phf1$Age)); d_phf1$PMI_s <- as.numeric(scale(d_phf1$PMI))
d_phf1$Sex <- droplevels(factor(d_phf1$Sex)); d_phf1$sample_id <- factor(d_phf1$sample_id)

cat(sprintf("\nDistance frame : %d cells, %d donors, dist_sd = %.6f\n",
            nrow(d_dist), dplyr::n_distinct(d_dist$sample_id), dist_sd))
cat(sprintf("PHF1+/- frame  : %d cells (%d PHF1+, %d PHF1-), %d donors\n",
            nrow(d_phf1), sum(d_phf1$phf1_pos), sum(!d_phf1$phf1_pos),
            dplyr::n_distinct(d_phf1$sample_id)))
if (sum(d_phf1$phf1_pos) == 0)
  stop("No PHF1+ cells in the contrast frame for ", CT,
       " -- the PHF1+/- decomposition cannot be fitted.")

FRAMES <- list(dist = list(d = d_dist, E = E_dist, rows = rows_d),
               phf1 = list(d = d_phf1, E = E_phf1, rows = rows_p))
stopifnot("frame row index out of step with its data" =
            length(rows_d) == nrow(d_dist) && length(rows_p) == nrow(d_phf1))

##  ............................................................................
##  Fit every gene once per contrast, then assemble each module                ####
pct <- function(x, ref) if (is.finite(ref) && ref != 0) 100 * abs(x) / abs(ref) else NA_real_

GENE_CACHE <- list()
for (tk in names(TERMS)) {
  tm <- TERMS[[tk]]
  cat("\n=== Fitting", length(GENE_UNION), "genes for the", tm$label, "===\n")
  GENE_CACHE[[tk]] <- fit_gene_cache(FRAMES[[tk]]$d, FRAMES[[tk]]$E, GENE_UNION,
                                     tm$predictor, paired = tm$paired)
}

RES <- list()                                   # RES[[module]][[contrast]]
for (m in MODS) {
  RES[[m]] <- list()
  for (tk in names(TERMS)) {
    tm <- TERMS[[tk]]
    av <- AGG_VEC[[m]]; sel <- FRAMES[[tk]]$rows
    r <- decompose_module(FRAMES[[tk]]$d, GENE_CACHE[[tk]], prep$present[[m]],
                          av$module_score[sel], av$geneset_mean[sel], av$control_mean[sel],
                          tm$predictor)
    r$gt$direction <- tm$dir_fun(r$gt$estimate)
    RES[[m]][[tk]] <- r
  }
  cat(sprintf("  %-16s dist: module %+.6f control-share %5.1f%%  |  phf1: module %+.6f control-share %5.1f%%\n",
              m, RES[[m]]$dist$agg$module_score$estimate,
              pct(RES[[m]]$dist$agg$control_mean$estimate, RES[[m]]$dist$agg$geneset_mean$estimate),
              RES[[m]]$phf1$agg$module_score$estimate,
              pct(RES[[m]]$phf1$agg$control_mean$estimate, RES[[m]]$phf1$agg$geneset_mean$estimate)))
}

##  ............................................................................
##  Figures + source data, per module x contrast                            ####
dec_table <- function(r) tibble(
  term = factor(c("Module score", "Geneset mean", "Control mean"),
                levels = c("Control mean", "Geneset mean", "Module score")),
  estimate = c(r$agg$module_score$estimate, r$agg$geneset_mean$estimate, r$agg$control_mean$estimate),
  ci_lo    = c(r$agg$module_score$ci_lo,    r$agg$geneset_mean$ci_lo,    r$agg$control_mean$ci_lo),
  ci_hi    = c(r$agg$module_score$ci_hi,    r$agg$geneset_mean$ci_hi,    r$agg$control_mean$ci_hi),
  p_LMM    = c(r$agg$module_score$p,        r$agg$geneset_mean$p,        r$agg$control_mean$p))

summary_rows <- list()
for (m in MODS) {
  MOD_LAB <- unname(MODULE_LABELS[m])
  for (tk in names(TERMS)) {
    tm <- TERMS[[tk]]; r <- RES[[m]][[tk]]; gt <- r$gt

    draw <- if (args$top_n > 0 && args$top_n < nrow(gt)) {
      o <- order(abs(gt$contribution), decreasing = TRUE)[seq_len(args$top_n)]
      gt[sort(o), ]
    } else gt
    draw$gene <- factor(draw$gene, levels = draw$gene[order(draw$estimate)])

    # Row count for the figure height must be taken BEFORE padding: the filler rows are
    # undrawn, so counting them would add blank space to every panel.
    n_drawn <- nrow(draw)

    # Complete the "Gene FDR" shape legend. A module whose genes are all n.s. would
    # otherwise render the "padj<0.05" key with a label and no point. drop = FALSE does
    # not fix it; see pad_levels() in R/palettes.R. `significant` is logical here, so the
    # levels are the character forms the scale's limits already use.
    draw <- pad_levels(draw, "significant", c("TRUE", "FALSE"))

    p_genes <- ggplot(draw, aes(x = estimate, y = gene, colour = direction)) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40", linetype = "dashed") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.4) +
      geom_point(aes(shape = significant), size = 1.3, stroke = 0.5) +
      scale_colour_manual(values = tm$colours, name = NULL, drop = FALSE) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Gene FDR",
                         labels = c(`TRUE` = "padj<0.05", `FALSE` = "ns"),
                         limits = c("TRUE", "FALSE"), drop = FALSE) +
      labs(x = tm$xlab, y = NULL) +
      forest_theme + theme(legend.position = "bottom", legend.box = "vertical",
                           legend.spacing.y = grid::unit(1, "pt"))
    ggsave(file.path(args$output_dir,
           sprintf("plot_genecontrib%s_%s_%s.pdf", tm$file, m, CT)), p_genes,
           width = 10, height = 2.4 + 0.30 * n_drawn, units = "cm", device = "pdf",
           limitsize = FALSE)

    dec <- dec_table(r)
    p_dec <- ggplot(dec, aes(estimate, term)) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40", linetype = "dashed") +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.5, colour = "grey20") +
      geom_point(size = 1.8, colour = "grey20", stroke = 0) +
      labs(x = tm$declab, y = NULL) + forest_theme
    ggsave(file.path(args$output_dir,
           sprintf("plot_genecontrib%s_decomposition_%s_%s.pdf", tm$file, m, CT)), p_dec,
           width = 8.0, height = 3.2, units = "cm", device = "pdf")

    write.table(gt %>% mutate(contrast = tm$key, module = m, module_label = MOD_LAB,
                              celltype = CT, reactome_id = unname(REACTOME_PRIMARY[m]),
                              k_genes = r$k, drawn = gene %in% as.character(draw$gene)) %>%
                  dplyr::select(contrast, module, module_label, reactome_id, celltype, k_genes,
                                gene, estimate, se, ci_lo, ci_hi, contribution, contribution_lo,
                                contribution_hi, pct_of_geneset, cohens_d, dz_donor_paired,
                                n_donors_paired, n_donors_same_dir, mean_expr, t, df, p_LMM,
                                lmm_padj, significant, direction, singular, drawn),
                file.path(args$output_dir,
                          sprintf("source_data_genecontrib%s_%s_%s.tsv", tm$file, m, CT)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(dec %>% mutate(contrast = tm$key, module = m, celltype = CT) %>%
                  dplyr::select(contrast, module, celltype, term, estimate, ci_lo, ci_hi, p_LMM),
                file.path(args$output_dir,
                  sprintf("source_data_genecontrib%s_decomposition_%s_%s.tsv", tm$file, m, CT)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    top <- gt[which.max(abs(gt$contribution)), ]
    summary_rows[[paste(m, tk)]] <- tibble(
      celltype = CT, module = m, module_label = MOD_LAB, contrast = tm$key,
      k_genes = r$k, n_genes_fitted = r$n_fitted,
      module_effect  = r$agg$module_score$estimate,
      geneset_effect = r$agg$geneset_mean$estimate,
      control_effect = r$agg$control_mean$estimate,
      control_share_pct = pct(r$agg$control_mean$estimate, r$agg$geneset_mean$estimate),
      sum_contributions = r$sum_contrib,
      resid_geneset = r$rec_geneset, resid_module = r$rec_module,
      n_genes_sig = sum(gt$significant %in% TRUE),
      top_gene = top$gene, top_contribution = top$contribution,
      top_pct_of_geneset = top$pct_of_geneset)
  }
}

SUMMARY <- bind_rows(summary_rows)
write.table(SUMMARY, file.path(args$output_dir,
            sprintf("stats_genecontrib_summary_%s.tsv", CT)),
            sep = "\t", quote = FALSE, row.names = FALSE)

##  ............................................................................
##  Stats log (all modules, both contrasts)                                 ####
sink(file.path(args$output_dir, sprintf("stats_genecontrib_%s.txt", CT)))
cat("Per-gene decomposition of module scores, BOTH contrasts\n")
cat("Celltype:", CT, "| modules:", length(MODS), "(", args$module, ")\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("MODULES DECOMPOSED\n")
print(as.data.frame(tibble(module = MODS,
                           reactome_id = unname(REACTOME_PRIMARY[MODS]),
                           label = unname(MODULE_LABELS[MODS]),
                           k_genes = lengths(prep$present[MODS]))))
if (identical(args$module, "significant"))
  cat("\nSelected as significant in the log-distance model for this celltype, read from\n  ",
      args$lmm_coeffs, "\n", sep = "")
cat("\nGenes were fitted ONCE over the union of these modules (", length(GENE_UNION),
    " unique genes,\nvs ", sum(lengths(prep$present[MODS])),
    " if refitted per module), so a gene shared between modules carries the\nSAME effect in each.",
    " BH is applied WITHIN each module, so its adjusted p can differ.\n\n", sep = "")

cat("MODELS (one per response; identical to the main script's two models):\n")
cat("  distance : y ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
cat("  PHF1+/-  : y ~ phf1_pos               + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
cat("  Distance transform (canonical): log(dist)/sd(log(dist)); dist_sd =", round(dist_sd, 6), "\n")
cat("  Distance frame: PHF1-NEGATIVE cells, cap",
    if (is.finite(MAX_DIST)) sprintf("%g um", MAX_DIST) else "none", "--",
    nrow(d_dist), "cells,", dplyr::n_distinct(d_dist$sample_id), "donors.\n")
cat("  PHF1+/- frame : ALL cells, NO distance cap and no distance filter --",
    nrow(d_phf1), "cells (", sum(d_phf1$phf1_pos), "PHF1+,", sum(!d_phf1$phf1_pos), "PHF1-),",
    dplyr::n_distinct(d_phf1$sample_id), "donors.\n")
cat("    PHF1+ cells have no dist_to_phf1_um, so a distance filter would delete the group\n")
cat("    being contrasted; a cap would select the PHF1-negative comparison group on\n")
cat("    proximity to PHF1+ cells and bias the contrast toward the null.\n")
cat("  AddModuleScore: SCT, ctrl =", CTRL, ", nbin =", NBIN, ", seed =", args$seed, ";",
    length(SCORED), "modules scored in ONE call so the scores match the main figures.\n\n")

cat("DECOMPOSITION\n")
cat("  score = mean(geneset genes) - mean(control genes), so for EITHER contrast\n")
cat("    effect(score)        = effect(geneset_mean) - effect(control_mean)\n")
cat("    effect(geneset_mean) = (1/k) * SUM_g effect_g,  contribution_g = effect_g / k\n\n")

cat("CROSS-MODULE SUMMARY (also in stats_genecontrib_summary_", CT, ".tsv)\n", sep = "")
print(as.data.frame(SUMMARY %>% dplyr::transmute(
  module, contrast, k = k_genes,
  module_eff = signif(module_effect, 3), geneset_eff = signif(geneset_effect, 3),
  control_eff = signif(control_effect, 3), ctrl_share = round(control_share_pct, 1),
  n_sig = n_genes_sig, top_gene, top_pct = round(top_pct_of_geneset, 1))))
cat("\nctrl_share = |control effect| as a percentage of |geneset effect|: the part of the shift\n")
cat("shared with the expression-matched control genes.\n")

for (m in MODS) {
  cat("\n\n==========================================================================\n")
  cat(m, " -- ", unname(MODULE_LABELS[m]), " (", unname(REACTOME_PRIMARY[m]), ")\n", sep = "")
  cat("==========================================================================\n")
  for (tk in names(TERMS)) {
    tm <- TERMS[[tk]]; r <- RES[[m]][[tk]]
    cat("\n--- ", toupper(tm$label), " ---\n", sep = "")
    if (tk == "dist") cat("SIGN: per +1 SD of log distance, so NEGATIVE = higher NEAR PHF1+.\n\n")
    else              cat("SIGN: PHF1+ minus PHF1-, so POSITIVE = higher IN PHF1+ cells.\n\n")
    print(as.data.frame(dec_table(r)))
    cat("\n  k =", r$k, "genes |", r$n_fitted, "fitted | sum of contributions",
        format(r$sum_contrib, digits = 6), "\n")
    cat("  RECONSTRUCTION RESIDUALS:\n")
    cat("    sum(contributions) - effect(geneset_mean)  =", format(r$rec_geneset, digits = 3),
        sprintf("(%.2f%% of the geneset effect)\n", pct(r$rec_geneset, r$agg$geneset_mean$estimate)))
    cat("    [geneset - control] - effect(module_score) =", format(r$rec_module, digits = 3),
        sprintf("(%.2f%% of the module effect)\n", pct(r$rec_module, r$agg$module_score$estimate)))
    cat("  CONTROL-GENE SHARE:",
        sprintf("%.1f%%", pct(r$agg$control_mean$estimate, r$agg$geneset_mean$estimate)),
        "of the geneset effect.\n\n")
    tab <- r$gt %>% dplyr::transmute(
      gene, estimate = round(estimate, 5),
      ci95 = sprintf("[%.5f, %.5f]", ci_lo, ci_hi),
      contribution = signif(contribution, 3), pct_of_geneset = round(pct_of_geneset, 1),
      cohens_d = round(cohens_d, 4), mean_expr = round(mean_expr, 3),
      padj = signif(lmm_padj, 3), significant)
    if (tm$paired) tab <- tab %>% mutate(dz = round(r$gt$dz_donor_paired, 3),
                                         same_dir = r$gt$n_donors_same_dir,
                                         n_donors = r$gt$n_donors_paired)
    print(as.data.frame(tab))
  }
}

cat("\n\n==========================================================================\n")
cat("=== UNITS ===\n")
cat("Every response here -- each gene, the geneset mean, the control mean and the module\n")
cat("score -- is on ONE scale: the SCT 'data' layer, which Seurat defines as\n")
cat("    log1p(corrected counts) = log_e(corrected UMI + 1)\n")
cat("where 'corrected' means SCT's counts reverse-transformed to a fixed sequencing depth.\n")
cat("AddModuleScore reads that same layer (its slot default is 'data'), so gene-level and\n")
cat("module-level numbers are directly comparable and the contributions really do add up.\n\n")
cat("So an estimate of, say, +0.10 means the adjusted mean of log_e(corrected UMI + 1) is\n")
cat("0.10 higher in PHF1+ cells (or per +1 SD of log distance, for the distance contrast).\n\n")
cat("DO NOT read exp(estimate) as a fold change. That identity needs log1p ~ log, i.e.\n")
cat("counts comfortably above 1. Much of a CosMx panel sits at 0-3 counts per cell, where\n")
cat("log1p is markedly compressive and exp(estimate) UNDERSTATES the multiplicative change.\n")
cat("The mean_expr column is in these same units and tells you which regime a gene is in:\n")
cat("the approximation is tolerable for high mean_expr, misleading for low.\n")
cat("If a fold change is needed for a specific gene, get it from the DEG pipeline\n")
cat("(deg/de_linear_distance/, logCPM via voom+dream), which is built for that -- but note\n")
cat("its coefficients are on a different scale and CANNOT be summed into these module\n")
cat("contributions, which is the whole reason this script refits in SCT space.\n\n")

cat("=== Effect sizes and how to read them ===\n")
cat("estimate is unstandardised, in the log1p SCT units above, with its 95% CI. cohens_d =\n")
cat("estimate / sqrt(donor variance + residual). contribution = estimate/k, in the SAME\n")
cat("units as the module-score effect, and the contributions SUM to it.\n")
cat("pct_of_geneset can exceed 100 or be negative: genes opposing the module direction\n")
cat("subtract from the total, so the shares do not partition a positive whole.\n")
cat("For the PHF1+/- contrast, dz is the donor-paired effect (mean per-donor difference /\n")
cat("SD of those differences) and same_dir counts donors moving with the mean.\n\n")

cat("Both reconstruction residuals are non-zero only because REML re-estimates the variance\n")
cat("components for every response, so the LMM is not exactly linear in y. Judge them on the\n")
cat("ABSOLUTE residual first: the module score is a DIFFERENCE of two larger effects, so\n")
cat("whenever it is small the second percentage is inflated by that cancellation and can read\n")
cat("alarmingly high while the absolute error is negligible (in simulation: 6.6e-05 absolute\n")
cat("= 2.1% of a module effect of -0.003, against 0.25% on the geneset identity). A residual\n")
cat("that is a material fraction of the GENESET effect is the one worth worrying about.\n\n")

cat("NOTE: the genes are the module's members, so the decomposition is descriptive; rank by\n")
cat("  contribution and read the CI.\n")
cat("NOTE (shared genes): a gene shared across several modules can appear as a top\n")
cat("  contributor in each of them; the summary table's top_gene column shows this across\n")
cat("  modules.\n")

cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

cat("\nDone.", length(MODS), "module(s) x", length(TERMS), "contrast(s).\n")
print(as.data.frame(SUMMARY %>% dplyr::transmute(
  module, contrast, module_eff = signif(module_effect, 3),
  ctrl_share = round(control_share_pct, 1), n_sig = n_genes_sig, top_gene)))
cat("Outputs written to:", args$output_dir, "\n")
