#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# deg_dream_linear_pseudotime_phf1.r
#
# Produces: Table S4
# Genes vs pseudotime along the CBLN2 trajectory; upstream of pathway_enrichment_cbln2_trajectory.R
# (Table S5) and the 2G dotplot.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# deg_dream_linear_pseudotime_phf1.r
#
# Step 5 of R/run_traj_pipeline.sh: genes that change along the CBLN2 tangle-associated
# trajectory fitted by R/cbln2_slingshot_trajectory.R. The counterpart of Foiani et al. 2026's
# "Correlation analysis" (per protein, Pearson vs pseudotime, BH, FDR < 0.05), done with the
# project's canonical cell-level mixed model so it is consistent with every other DEG here.
#
# Model (per gene) -- docs/MODELS.md Set 1, identical to R/deg_dream_linear_distance_phf1.r with
# the distance predictor swapped for pseudotime:
#   ~ pt_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1 | sample_id)
# voomWithDreamWeights on CPM (no TMM) + dream + eBayes, BH. pt_scaled = pseudotime / SD
# (divided, not centred -- the distance convention). logFC = logCPM change per SD of
# pseudotime: POSITIVE = rises toward the tangle-rich end. logFC_per_unit is per unit of
# slingshot pseudotime.
#
# Differences from the distance script, all structural:
#   * ALL cells of the celltype, PHF1+ included (Foiani correlates across all samples; the
#     tangle-bearing cells are the far end of the axis). The distance script models
#     PHF1-negative cells only because distance is undefined for PHF1+ cells.
#   * No distance cap / exclusion radius (no distance here).
#   * Input SCE is celltype_sce/<ct>_sce.qs (holds the PHF1+ cells), joined to pseudotime by
#     cell_id from plots/cbln2_slingshot_trajectory/source_data_cbln2_traj_percell.tsv.
# Gene filter as canonical: detected in >= 5% of modelled cells.
#
# OUTPUT  deg/de_linear_pseudotime/<celltype>/
#   <celltype>_pt_scaled.tsv   gene, logFC, logFC_per_unit, logFC_trajectory, CI.L, CI.R, CI.L_per_unit,
#                              CI.R_per_unit, CI.L_trajectory, CI.R_trajectory, pval, padj, AveExpr, t, B,
#                              contrast, model, n_cells, pt_sd, pt_range, fit_space, k, root
#     logFC            log2 CPM change per 1 SD of pseudotime (the model coefficient)
#     logFC_per_unit   per 1 pseudotime unit (= logFC / pt_sd)
#     logFC_trajectory across the whole trajectory, root to tip (= logFC_per_unit x pt_range); the
#                      scale the enrichment |log2FC| > 0.25 filter uses. A linear slope multiplied
#                      out, not an observed start-to-end difference.
#   Same column names as the distance DEG tables, so the enrichR script reads it.
#
# RUN:
#   Rscript R/deg_dream_linear_pseudotime_phf1.r --ncores 8                      # CBLN2
#   Rscript R/deg_dream_linear_pseudotime_phf1.r --celltype Exc-IT-L3-5-CHGA-IL1RAPL2 --ncores 8
# or as step 5 of R/run_traj_pipeline.sh.

suppressPackageStartupMessages({
  library(tidyverse); library(argparse); library(SingleCellExperiment); library(SummarizedExperiment)
  library(Matrix); library(edgeR); library(limma); library(variancePartition); library(BiocParallel)
})
.muffleReformulasMigrationWarning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("findbars|nobars", conditionMessage(w)) && grepl("reformulas", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/cbln2_embedding_utils.R")   # ct_sce_path(), ct_percell_file() (no packages needed to source)

parser <- ArgumentParser()
parser$add_argument("--celltype", default = CBLN2_CT,
                    help = "celltype; sets the --sce and --pseudotime defaults (CBLN2 paths unchanged)")
parser$add_argument("--sce", default = NULL, help = "default: celltype_sce/<celltype>_sce.qs")
parser$add_argument("--pseudotime", default = NULL,
                    help = "default: plots/<tag>_slingshot_trajectory/source_data_<tag>_traj_percell.tsv")
parser$add_argument("--confounding_vars", default = "Sex,Age,PMI")
parser$add_argument("--output_dir", default = "deg")
parser$add_argument("--ncores", type = "integer", default = 4L)
parser$add_argument("--min_pct_cells", type = "double", default = 5)
args <- parser$parse_args()
conf_vars <- trimws(strsplit(args$confounding_vars, ",")[[1]])
if (is.null(args$sce))        args$sce        <- ct_sce_path(args$celltype)
if (is.null(args$pseudotime)) args$pseudotime <- ct_percell_file(args$celltype)

celltype <- gsub("_sce\\.qs$", "", basename(args$sce))
outdir <- file.path(args$output_dir, "de_linear_pseudotime", celltype)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("Celltype:", celltype, "\nOutput dir:", outdir, "\n")

## Load + join pseudotime --------------------------------------------------------------------
sce <- qs::qread(args$sce)
rownames(sce) <- rowData(sce)$gene
pt_tab <- read.delim(args$pseudotime, colClasses = c(cell_id = "character", sample_id = "character"))
stopifnot(!anyDuplicated(pt_tab$cell_id), all(c("pt", "fit_space", "k", "root") %in% names(pt_tab)))
m <- match(colnames(sce), pt_tab$cell_id)
if (anyNA(m)) stop(sum(is.na(m)), " cells in the SCE have no pseudotime -- re-run R/cbln2_slingshot_trajectory.R on this SCE.")
if (anyNA(pt_tab$pt[m])) stop("NA pseudotime for ", sum(is.na(pt_tab$pt[m])), " cells.")
stopifnot(identical(as.character(sce$sample_id), pt_tab$sample_id[m]))
sce$pt <- pt_tab$pt[m]
traj <- unique(pt_tab[, c("fit_space", "k", "root")]); stopifnot(nrow(traj) == 1)
cat(sprintf("Trajectory: space %s, k %d, root %s. Cells %d (PHF1+ %d), donors %d\n",
            traj$fit_space, traj$k, traj$root, ncol(sce), sum(as.logical(as.character(sce$PHF1))),
            length(unique(sce$sample_id))))
cat("Cells per sample:\n"); print(table(sce$sample_id))

## Gene filter (canonical: detected in >= 5% of modelled cells) ------------------------------
expr_cells <- rowSums(counts(sce) > 0)
bkg_genes  <- rownames(sce)[expr_cells >= (args$min_pct_cells / 100) * ncol(sce)]
cat(sprintf("Genes: %d total, %d detected in >= %s%% of %d cells\n", nrow(sce), length(bkg_genes),
            format(args$min_pct_cells), ncol(sce)))

## Covariates (identical to the distance script) ----------------------------------------------
cd <- as.data.frame(colData(sce))
cd$nUMI_log    <- log2(colSums(counts(sce)) + 1)
cd$percent_neg <- as.numeric(cd$percent.neg)
cd$Age <- as.numeric(scale(as.numeric(cd$Age)))
cd$PMI <- as.numeric(scale(as.numeric(cd$PMI)))
cd$sample_id <- droplevels(as.factor(cd$sample_id))
pt_sd <- sd(cd$pt)
pt_range <- diff(range(cd$pt))       # whole trajectory, for logFC_trajectory
cd$pt_scaled <- cd$pt / pt_sd
cat(sprintf("Pseudotime: range %.2f-%.2f, SD %.4f (used for scaling)\n", min(cd$pt), max(cd$pt), pt_sd))

form_str <- sprintf("~ pt_scaled + nUMI_log + percent_neg + %s + (1|sample_id)", paste(conf_vars, collapse = " + "))
form <- as.formula(form_str)
missing <- setdiff(setdiff(all.vars(form), "sample_id"), colnames(cd))
if (length(missing)) stop("Missing colData columns: ", paste(missing, collapse = ", "))
cat("Dream formula:", form_str, "\n")

cat("\nDiagnostic | Pearson r between pseudotime and QC/confounders:\n")
for (v in intersect(c("nUMI_log", "percent_neg", conf_vars), colnames(cd))) if (is.numeric(cd[[v]])) {
  r <- cor(cd$pt, cd[[v]], use = "complete.obs")
  cat(sprintf("  %-12s r = %+.3f%s\n", v, r, ifelse(abs(r) > 0.3, "  <-- WARNING: collinear", "")))
}

## voom + dream ---------------------------------------------------------------------------------
dge   <- DGEList(counts = counts(sce)[bkg_genes, ])
param <- MulticoreParam(args$ncores, progressbar = TRUE)
vobj  <- .muffleReformulasMigrationWarning(voomWithDreamWeights(counts = dge, formula = form, data = cd, BPPARAM = param))
fit   <- .muffleReformulasMigrationWarning(dream(exprObj = vobj, formula = form, data = cd, BPPARAM = param))
fit   <- eBayes(fit)

res <- topTable(fit, coef = "pt_scaled", number = nrow(fit), adjust.method = "BH", confint = TRUE) %>%
  mutate(gene = rownames(.), contrast = "pt_scaled", model = form_str,
         logFC_per_unit = logFC / pt_sd, CI.L_per_unit = CI.L / pt_sd, CI.R_per_unit = CI.R / pt_sd,
         logFC_trajectory = logFC_per_unit * pt_range, CI.L_trajectory = CI.L_per_unit * pt_range,
         CI.R_trajectory = CI.R_per_unit * pt_range,
         n_cells = ncol(sce), pt_sd = pt_sd, pt_range = pt_range, fit_space = traj$fit_space, k = traj$k, root = traj$root) %>%
  dplyr::rename(pval = P.Value, padj = adj.P.Val) %>%
  dplyr::select(gene, logFC, logFC_per_unit, logFC_trajectory, CI.L, CI.R, CI.L_per_unit, CI.R_per_unit,
                CI.L_trajectory, CI.R_trajectory, pval, padj,
                AveExpr, t, B, contrast, model, n_cells, pt_sd, pt_range, fit_space, k, root)
out_file <- file.path(outdir, paste0(celltype, "_pt_scaled.tsv"))
write.table(res, out_file, col.names = TRUE, row.names = FALSE, sep = "\t")
cat("\nWritten:", out_file, "\n")

## Summary ----------------------------------------------------------------------------------------
for (thr in c(0.05, 0.1)) {
  s <- res %>% filter(padj < thr)
  cat(sprintf("padj < %.2f: %d genes (%d rise toward the tangle-rich end, %d fall)%s\n", thr, nrow(s),
              sum(s$logFC > 0), sum(s$logFC < 0),
              if (thr == 0.05) "  [Foiani threshold]" else "  [project DEG threshold, feeds enrichment]"))
}
cat("Share of tested genes at padj < 0.1:", sprintf("%.0f%%", 100 * mean(res$padj < 0.1)), "\n")
cat(sprintf("rho(logFC, AveExpr) over all genes: %.3f\n",
            cor(res$logFC, res$AveExpr, method = "spearman")))
cat("\nTop 20 rising / falling (padj < 0.1, by t):\n")
print(as.data.frame(res %>% filter(padj < 0.1) %>% arrange(desc(t)) %>% head(20) %>% select(gene, logFC, CI.L, CI.R, padj)), digits = 3)
print(as.data.frame(res %>% filter(padj < 0.1) %>% arrange(t) %>% head(20) %>% select(gene, logFC, CI.L, CI.R, padj)), digits = 3)
cat("\n--- sessionInfo ---\n"); print(sessionInfo())
