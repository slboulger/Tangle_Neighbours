#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# prepare_seg_geneset_reference.R
#
# Figure panels: S5A-C - upstream step
# Builds the stably expressed gene reference (seg_control/Reference_SEG_panel.qs)
# used as the negative-control module.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# prepare_seg_geneset_reference.R
#
# ONE-OFF, THREE STAGES. Derives PROJECT-SPECIFIC stably expressed gene (SEG) sets
# from our own snRNA-seq reference, one per CosMx celltype, with
# scMerge::scSEGIndex - the same method Lin et al. 2019 used for the published human
# SEG list, but fitted on this tissue and these cell types.
#
#   Lin Y, Ghazanfar S, Strbenac D, et al. Evaluating stably expressed genes in
#   single cells. GigaScience 2019;8(9):giz106.  (scSEGIndex)
#
# WHY THREE STAGES (memory)
#   A single process is memory-bound for two structural reasons, neither to do with
#   the fitting itself:
#     (1) BiocParallel MulticoreParam FORKS the whole R process. Forking while the
#         parent held the 124k-cell reference plus subset copies plus a new logcounts
#         assay multiplied a multi-GB heap by the worker count.
#     (2) Chained subsetting (sce[genes, ] then sce[, cells]) allocates a full copy
#         per step while the previous copy is still reachable.
#   So: stage `extract` runs ONCE with NO forking and does gene+cell subsetting in a
#   SINGLE indexing operation; stage `fit` runs in a FRESH process per celltype and
#   only ever loads that celltype's small matrix; stage `finalise` just assembles.
#   Default --workers 1 => SerialParam => no forking anywhere. Raise it only if the
#   measured per-gene cost makes serial fitting impractical, and know the risk.
#
# REFERENCE
#   .../MAP_analysis/IGFQ002014_Feb26/sce_ref.qs - 25454 genes x 124251 cells,
#   36 donors, counts assay only, ENSEMBL rownames, symbols in rowData$gene.
#   Built by that directory's R/create_ref_sce.R: colData$subcluster = neuron
#   subclusters from harmony_integration_Neuron/neuron_sce3.qs, glia labelled directly.
#
# LABEL MAPPING  reference subcluster names are NOT the CosMx names. The chain is in
#   the project's own scripts and is reproduced in CT_MAP below:
#     R/AUCell_Neuron.R:29-30   subcluster -> broad_subcluster (case_when)
#     R/rename_celltypes.R      broad label -> final CosMx `celltype`
#   e.g. Exc-L2-3-CUX2-CALB1 + Exc-L2-CUX2-PDGFD + Exc-L3-PCP4-CALB1 +
#        Exc-L2-RELN-BMPR1B  ->  Exc-IT-L2-3  ->  Exc-IT-L2-3-CBLN2-HOPX
#
# GENE UNIVERSE  scSEGIndex's segIdx is built from RANKS over the genes supplied
#   (x1 = rank(rho)/(n+1), x2 = 1 - rank(sigma)/(n+1), x3 = 1 - rank(zero)/(n+1)) and
#   mu.scaled is min-max scaled over the same set, so the universe changes every
#   score. Default universe is `panel`: the reference restricted to the CosMx panel
#   BEFORE fitting, so genes are ranked among those CosMx can actually measure.
#   `genomewide` (all reference genes, intersected with the panel afterwards) is
#   supported but NOT run by default.
#
# SELECTION RULE  applied in stage `finalise`, NOT in `fit`, so the rule can be
#   re-cut from the cached fits in seconds without refitting anything.
#     default            : segIdx > 80th percentile          (--feature_gates no)
#     optional, stricter : ... AND x1, x2, x3 each > 60th pct (--feature_gates yes)
#   scSEGIndex returns only the absolute parameters, so x1/x2/x3 are recomputed with
#   its exact internal formulas and asserted to reproduce the returned segIdx.
#
#   WHY THE PER-FEATURE GATES ARE OFF BY DEFAULT (measured on the panel universe):
#   the three-gate conjunction collapses for Micro. Selected counts were
#     Exc-IT-L2-3 155 | Exc-IT-L3-5 95 | Astro 100 | Micro 15
#   because Micro's stability features are anti-correlated where the other celltypes'
#   agree - Spearman x2(sigma)~x3(omega*) is +0.71/+0.88/+0.86 for Astro/L2-3/L3-5 but
#   -0.47 for Micro, so its sigma and omega* gates are near mutually exclusive. A
#   15-gene control module is far noisier than the others (AddModuleScore variance
#   scales ~1/n_genes). segIdx alone gives 305/337/284/237 - comparable and adequately
#   powered. Note also that given the three Q60 gates the segIdx cut was nearly
#   redundant (it removed only 6/29/20/0 further genes), i.e. the conjunction was
#   effectively "all three > Q60" rather than Lin's percentile on the index.
#   The x1/x2/x3 values and their pass flags are still REPORTED per gene either way.
#
# ABSOLUTE STABILITY FEATURES  reported un-ranked in seg_index_<universe>_<CT>.tsv:
#     lambda   = rho         mixing proportion of the normal ("expressed") component
#     sigma                  SD of the normal component (scMerge returns the SD)
#     sigma_sq = sigma^2     variance of the normal component
#     omega                  raw dropout rate (fraction of zero cells)
#     omega_st = zero        = omega * mu.scaled, the mean-scaled dropout ("omega*")
#     mu, mu_scaled          normal-component mean and its min-max scaled value
#     F_donor, p_donor       the F-STATISTIC REPLACEMENT, below
#
# F-STATISTIC REPLACEMENT  scSEGIndex's optional 4th feature is a between-CELLTYPE
#   ANOVA F, which penalises genes differing between cell types. It is not applicable
#   here because each set is fitted WITHIN one celltype (cell_type = NULL). The
#   informative analogue for a control that must be flat across our 9 CosMx donors is
#   a between-DONOR ANOVA F within that celltype (36 reference donors):
#       aov(logcounts_gene ~ donor)
#   High F_donor = donor-dependent = NOT stable. It is REPORTED as a diagnostic and is
#   deliberately NOT in the selection rule (the rule is Lin's three-feature version,
#   matching cell_type = NULL).
#
# OUTPUT
#   <out_dir>/Reference_SEG_<universe>.qs      list(<CosMx celltype> = <symbols>)
#   <out_dir>/seg_index_<universe>_<CT>.tsv    full per-gene table
#   <out_dir>/Reference_SEG_provenance.txt
#   <work_dir>/mat_<universe>_<CT>.qs          intermediate (safe to delete after)
#
# RUN (locally; requires scMerge).
#   Rscript R/prepare_seg_geneset_reference.R --stage extract  --work_dir <local dir>
#   Rscript R/prepare_seg_geneset_reference.R --stage fit      --work_dir <local dir> --celltype Micro
#   ... one fit process per celltype ...
#   Rscript R/prepare_seg_geneset_reference.R --stage finalise --work_dir <local dir>
# Use a LOCAL work_dir (not the RDS mount): the intermediates are a few hundred MB.

suppressPackageStartupMessages({
  library(qs)
  library(argparse)
})

PROJECT_ROOTS <- c(
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>/phf1_v2"
)
.root <- PROJECT_ROOTS[dir.exists(PROJECT_ROOTS)]
if (length(.root) == 0) stop("No project root found.", call. = FALSE)
setwd(.root[1])

REF_PATHS <- c(
  "<RDS_ROOT>/MAP_analysis/IGFQ002014_Feb26/sce_ref.qs",
  "<RDS_ROOT>/MAP_analysis/IGFQ002014_Feb26/sce_ref.qs"
)

parser <- ArgumentParser()
parser$add_argument("--stage", required = TRUE,
  help = "extract | fit | finalise")
parser$add_argument("--work_dir", required = TRUE,
  help = "LOCAL scratch dir for intermediates (a few hundred MB); not the RDS mount")
parser$add_argument("--celltype", default = "",
  help = "stage fit: which CosMx celltype to fit (one per process)")
parser$add_argument("--sce_ref", default = "",
  help = "Reference SCE .qs; default resolves the HPC path then the /Volumes mount")
parser$add_argument("--panel_sce",
  default = "celltype_sce_neighbours/Exc-ET-L5-SPON1-FGD4_sce_neighbours.qs",
  help = "Any per-celltype CosMx SCE; its rownames define the panel (cheap to read)")
parser$add_argument("--out_dir", default = "seg_control",
  help = "Output dir, relative to the project root [default: seg_control]")
parser$add_argument("--universe", default = "panel",
  help = "panel | genomewide [default: panel]")
parser$add_argument("--seg_q",  type = "double", default = 0.80,
  help = "segIdx percentile cut [default: 0.80]")
parser$add_argument("--feat_q", type = "double", default = 0.60,
  help = "Per-feature percentile cut on x1/x2/x3, used only with --feature_gates yes [default: 0.60]")
parser$add_argument("--feature_gates", default = "no",
  help = paste("yes|no - also require x1,x2,x3 each > --feat_q. Default no: the",
               "conjunction collapses for Micro (n=15). See the header. [default: no]"))
parser$add_argument("--zero_max", type = "double", default = 0.8,
  help = "Dropout filter; genes with a higher zero fraction are dropped [default: 0.8]")
parser$add_argument("--workers", type = "integer", default = 1,
  help = "1 = SerialParam (NO forking, recommended). >1 forks - see header [default: 1]")
parser$add_argument("--pilot_genes", type = "integer", default = 0,
  help = "stage fit: fit only the first N genes, to measure s/gene (0 = all)")
parser$add_argument("--seed", type = "integer", default = 42, help = "Seed")
args <- parser$parse_args()

set.seed(args$seed)
stopifnot("--stage must be extract, fit or finalise" =
            args$stage %in% c("extract", "fit", "finalise"))
stopifnot("--universe must be panel or genomewide" =
            args$universe %in% c("panel", "genomewide"))
stopifnot("--feature_gates must be yes or no" = args$feature_gates %in% c("yes", "no"))
USE_GATES <- identical(args$feature_gates, "yes")
dir.create(args$work_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(args$work_dir)) stop("could not create work_dir: ", args$work_dir, call. = FALSE)

# CosMx celltype <- reference `subcluster` components.
# Provenance: R/AUCell_Neuron.R:29-30 then R/rename_celltypes.R.
CT_MAP <- list(
  "Exc-IT-L2-3-CBLN2-HOPX"    = c("Exc-L2-3-CUX2-CALB1", "Exc-L2-CUX2-PDGFD",
                                  "Exc-L3-PCP4-CALB1",   "Exc-L2-RELN-BMPR1B"),
  "Exc-IT-L3-5-CHGA-IL1RAPL2" = c("Exc-L5-RORB-TLL1",    "Exc-L5-RORB-TPBG"),
  "Micro"                     = "Micro",
  "Astro"                     = "Astro"
)
ct_tag <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)
mat_path <- function(ct) file.path(args$work_dir,
                                   sprintf("mat_%s_%s.qs", args$universe, ct_tag(ct)))
res_path <- function(ct) file.path(args$work_dir,
                                   sprintf("res_%s_%s.qs", args$universe, ct_tag(ct)))

report_mem <- function(tag) {
  gc(verbose = FALSE)
  cat(sprintf("    [mem] %s: R heap ~%.2f GB\n", tag,
              sum(gc(verbose = FALSE)[, "used"] * c(8, 8)) / 1024^3))
}

##  ==========================================================================
##  STAGE extract - one process, NO forking, single indexing operation      ####
if (args$stage == "extract") {
  suppressPackageStartupMessages({
    library(SingleCellExperiment); library(scater); library(Matrix)
  })
  ref_path <- if (nzchar(args$sce_ref)) args$sce_ref else {
    h <- REF_PATHS[file.exists(REF_PATHS)]
    if (length(h) == 0) stop("sce_ref.qs not found at:\n  ",
                             paste(REF_PATHS, collapse = "\n  "), call. = FALSE)
    h[1]
  }
  cat("STAGE extract\n  reference:", ref_path, "\n  universe :", args$universe, "\n")

  panel_genes <- rownames(qs::qread(args$panel_sce))
  cat("  CosMx panel genes:", length(panel_genes), "\n")

  cat("  loading reference (large; single read, no forking)...\n")
  sce <- qs::qread(ref_path, nthreads = 4)
  cat("   ", nrow(sce), "genes x", ncol(sce), "cells\n")
  stopifnot("counts assay missing"      = "counts" %in% assayNames(sce),
            "subcluster column missing" = "subcluster" %in% colnames(colData(sce)),
            "manifest (donor) missing"  = "manifest"   %in% colnames(colData(sce)))
  report_mem("after read")

  # Library sizes from the FULL transcriptome, BEFORE any gene subsetting, so size
  # factors are not distorted by restricting to the panel.
  lib_size <- Matrix::colSums(counts(sce))
  cat("  library sizes: median", format(stats::median(lib_size), big.mark = ","), "\n")

  # ENSEMBL -> symbol; keep one ENSEMBL id per symbol, preferring highest total_counts.
  rd <- as.data.frame(rowData(sce))
  stopifnot("rowData$gene (symbols) missing" = "gene" %in% colnames(rd))
  sym <- as.character(rd$gene)
  tc  <- if ("total_counts" %in% colnames(rd)) as.numeric(rd$total_counts) else rep(0, nrow(rd))
  ok  <- !is.na(sym) & nzchar(sym)
  n_dup <- sum(duplicated(sym[ok]))
  ord   <- order(ifelse(ok, 0L, 1L), -tc)
  gsel  <- ord[!duplicated(sym[ord]) & ok[ord]]
  cat(sprintf("  symbols: %d/%d genes, %d duplicate symbols collapsed, %d retained\n",
              sum(ok), nrow(sce), n_dup, length(gsel)))

  gene_sym <- sym[gsel]
  if (args$universe == "panel") {
    keep_g <- gsel[gene_sym %in% panel_genes]
    cat("  panel genes present in reference:", length(keep_g), "of", length(panel_genes), "\n")
    stopifnot("too few panel genes in the reference" = length(keep_g) >= 1000)
  } else {
    keep_g <- gsel
  }

  sub_lab <- as.character(sce$subcluster)
  ct_of   <- rep(NA_character_, length(sub_lab))
  for (ct in names(CT_MAP)) {
    miss <- setdiff(CT_MAP[[ct]], unique(sub_lab))
    if (length(miss) > 0)
      stop(sprintf("reference subcluster(s) missing for %s: %s", ct,
                   paste(miss, collapse = ", ")), call. = FALSE)
    ct_of[sub_lab %in% CT_MAP[[ct]]] <- ct
  }
  keep_c <- which(!is.na(ct_of))
  cat("  cells per CosMx celltype:\n"); print(table(ct_of[keep_c]))

  # SINGLE indexing operation: genes and cells together, so only one copy is made.
  sce <- sce[keep_g, keep_c]
  rownames(sce) <- sym[keep_g]
  sce$cosmx_celltype <- ct_of[keep_c]
  sce$donor <- factor(as.character(sce$manifest))
  sizeFactors(sce) <- lib_size[keep_c] / mean(lib_size[keep_c])
  rm(lib_size, ct_of, sub_lab, rd, sym, tc, ord, gsel, gene_sym); invisible(gc())
  report_mem("after subset")

  cat("  normalising (logNormCounts on the precomputed full-transcriptome factors)\n")
  sce <- scater::logNormCounts(sce)
  report_mem("after logNormCounts")

  meta <- list(universe = args$universe, panel_genes = panel_genes,
               zero_max = args$zero_max, n_donors = nlevels(sce$donor),
               ref_path = ref_path, ct_map = CT_MAP,
               scmerge_version = as.character(utils::packageVersion("scMerge")))
  qs::qsave(meta, file.path(args$work_dir, sprintf("meta_%s.qs", args$universe)))

  for (ct in names(CT_MAP)) {
    idx <- which(sce$cosmx_celltype == ct)
    lc  <- logcounts(sce[, idx])
    z   <- Matrix::rowMeans(lc == 0)
    surv <- z <= args$zero_max
    cat(sprintf("  %-28s %6d cells | %5d genes -> %5d survive >%.0f%%-zero filter\n",
                ct, length(idx), length(z), sum(surv), 100 * args$zero_max))
    if (sum(surv) < 100)
      stop(sprintf("%s: only %d genes survive; scSEGIndex needs >=100", ct, sum(surv)),
           call. = FALSE)
    qs::qsave(list(celltype = ct,
                   lc = lc[surv, , drop = FALSE],          # sparse; densified in `fit`
                   donor = droplevels(sce$donor[idx]),
                   zero_all = as.numeric(z[surv])),
              mat_path(ct))
    rm(lc, z, surv); invisible(gc())
  }
  cat("\nSTAGE extract done. Intermediates in", args$work_dir, "\n")
}

##  ==========================================================================
##  STAGE fit - fresh process per celltype, loads only that small matrix    ####
if (args$stage == "fit") {
  suppressPackageStartupMessages({
    library(scMerge); library(BiocParallel); library(Matrix)
  })
  stopifnot("--celltype required for stage fit" = nzchar(args$celltype))
  ct <- args$celltype
  stopifnot("unknown --celltype" = ct %in% names(CT_MAP))
  if (!file.exists(mat_path(ct)))
    stop("intermediate not found: ", mat_path(ct), "\nRun --stage extract first.",
         call. = FALSE)

  cat("STAGE fit:", ct, "| universe:", args$universe, "\n")
  d  <- qs::qread(mat_path(ct))
  lc <- as.matrix(d$lc)                       # one dense copy, bounded and known
  donor <- d$donor
  zero_all <- setNames(d$zero_all, rownames(lc))
  rm(d); invisible(gc())
  if (args$pilot_genes > 0 && nrow(lc) > args$pilot_genes) {
    lc <- lc[seq_len(args$pilot_genes), , drop = FALSE]
    zero_all <- zero_all[rownames(lc)]
    cat("  PILOT: first", nrow(lc), "genes only\n")
  }
  cat(sprintf("  matrix: %d genes x %d cells (dense %.0f MB), donors: %d\n",
              nrow(lc), ncol(lc), as.numeric(object.size(lc)) / 1024^2, nlevels(donor)))
  report_mem("after load")

  BP <- if (args$workers <= 1) BiocParallel::SerialParam(progressbar = FALSE) else
    BiocParallel::MulticoreParam(workers = args$workers, progressbar = FALSE)
  cat("  BPPARAM:", class(BP)[1], if (args$workers > 1)
    sprintf("(%d forked workers - see header warning)", args$workers) else "(no forking)", "\n")

  t0  <- proc.time()[3]
  res <- scMerge::scSEGIndex(exprs_mat = lc, cell_type = NULL, BPPARAM = BP)
  el  <- proc.time()[3] - t0
  cat(sprintf("  scSEGIndex: %d genes in %.1f min (%.3f s/gene)\n",
              nrow(res), el / 60, el / max(nrow(res), 1)))
  report_mem("after fit")

  res <- res[!is.na(res$segIdx), , drop = FALSE]

  # Recompute scSEGIndex's three rank features with its exact internal formulas and
  # assert they reproduce the returned segIdx (guards a scMerge version change).
  n  <- nrow(res)
  x1 <- rank(res$rho) / (n + 1)                 # higher rho   = more stable
  x2 <- 1 - rank(res$sigma) / (n + 1)           # lower sigma  = more stable
  x3 <- 1 - rank(res$zero) / (n + 1)            # lower omega* = more stable
  chk <- max(abs(rowMeans(cbind(x1, x2, x3)) - res$segIdx), na.rm = TRUE)
  if (!is.finite(chk) || chk > 1e-8)
    stop(sprintf(paste0("recomputed x1/x2/x3 do not reproduce segIdx (max|diff|=%.3g); ",
                        "scSEGIndex internals may have changed - re-check the rank ",
                        "formulas before trusting the selection rule"), chk), call. = FALSE)
  cat(sprintf("  x1/x2/x3 reproduce segIdx (max|diff| = %.2g)\n", chk))

  # Between-donor ANOVA F within this celltype: the F-statistic replacement.
  # Closed-form one-way ANOVA, vectorised over genes.
  k <- nlevels(donor); nc <- ncol(lc)
  if (k >= 2) {
    ind   <- stats::model.matrix(~ donor - 1)
    n_g   <- colSums(ind)
    sums  <- lc[rownames(res), , drop = FALSE] %*% ind
    gmean <- sums / rep(n_g, each = nrow(sums))
    grand <- rowSums(lc[rownames(res), , drop = FALSE]) / nc
    SSB   <- rowSums(sweep(gmean, 1, grand, "-")^2 * rep(n_g, each = nrow(sums)))
    SST   <- rowSums(lc[rownames(res), , drop = FALSE]^2) - nc * grand^2
    SSW   <- pmax(SST - SSB, 0)
    Fv    <- (SSB / (k - 1)) / (SSW / (nc - k))
    res$F_donor <- as.numeric(Fv)
    res$p_donor <- stats::pf(res$F_donor, k - 1, nc - k, lower.tail = FALSE)
  } else {
    res$F_donor <- NA_real_; res$p_donor <- NA_real_
  }

  # NOTE: no selection here. `fit` only produces the per-gene feature table; the rule
  # is applied in `finalise` so it can be re-cut from these cached results in seconds.
  out <- data.frame(
    gene = rownames(res), celltype = ct, universe = args$universe,
    segIdx   = res$segIdx,
    lambda   = res$rho,                     # mixing proportion, normal component
    sigma    = res$sigma,                   # SD of the normal component
    sigma_sq = res$sigma^2,                 # variance
    mu       = res$mu,
    mu_scaled = res$mu.scaled,
    omega    = zero_all[rownames(res)],     # raw dropout rate
    omega_st = res$zero,                    # = omega * mu.scaled
    F_donor  = res$F_donor, p_donor = res$p_donor,
    x1_rho   = x1, x2_sigma = x2, x3_zero = x3,
    stringsAsFactors = FALSE, row.names = NULL)
  out <- out[order(-out$segIdx), ]

  qs::qsave(list(celltype = ct, table = out, n_surv = nrow(out),
                 seconds = el, pilot = args$pilot_genes),
            res_path(ct))
  cat(sprintf("  wrote %s (%d genes; selection deferred to --stage finalise)\n",
              res_path(ct), nrow(out)))
}

##  ==========================================================================
##  STAGE finalise - assemble the geneset .qs, TSVs and provenance          ####
if (args$stage == "finalise") {
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
  .probe <- file.path(args$out_dir, ".write_probe")
  if (!file.create(.probe, showWarnings = FALSE))
    stop("out_dir is not writable: ", args$out_dir, call. = FALSE)
  unlink(.probe)

  meta_p <- file.path(args$work_dir, sprintf("meta_%s.qs", args$universe))
  if (!file.exists(meta_p)) stop("meta not found: ", meta_p, call. = FALSE)
  meta <- qs::qread(meta_p)

  missing <- names(CT_MAP)[!file.exists(vapply(names(CT_MAP), res_path, ""))]
  if (length(missing) > 0)
    stop("no fit result for: ", paste(missing, collapse = ", "),
         "\nRun --stage fit for each celltype first.", call. = FALSE)

  RULE_LINES <- if (USE_GATES)
    sprintf("Rule      : segIdx > Q%.2f AND x1(rho), x2(sigma), x3(omega*) each > Q%.2f",
            args$seg_q, args$feat_q)
  else c(
    sprintf("Rule      : segIdx > Q%.2f  (Lin's percentile on the index). The optional", args$seg_q),
    "            per-feature Q60 gates are OFF: they collapse Micro to n=15 because its",
    "            sigma and omega* features are anti-correlated (Spearman x2~x3 = -0.47,",
    "            vs +0.71/+0.88/+0.86 in Astro/L2-3/L3-5), and given those gates the",
    "            segIdx cut was nearly redundant. Pass flags are still reported per gene.")

  cat(sprintf("Rule: segIdx > Q%.2f%s\n", args$seg_q,
              if (USE_GATES) sprintf(" AND x1,x2,x3 each > Q%.2f", args$feat_q) else
                " (per-feature gates OFF - see script header)"))

  sets <- list(); info <- list()
  for (ct in names(CT_MAP)) {
    r <- qs::qread(res_path(ct))
    if (r$pilot > 0)
      stop(sprintf("%s was fitted in PILOT mode (%d genes) - rerun the full fit",
                   ct, r$pilot), call. = FALSE)

    # ---- apply the selection rule here (cheap; no refitting) ----
    d <- r$table
    seg_cut <- as.numeric(stats::quantile(d$segIdx, args$seg_q, na.rm = TRUE))
    c1 <- as.numeric(stats::quantile(d$x1_rho,   args$feat_q, na.rm = TRUE))
    c2 <- as.numeric(stats::quantile(d$x2_sigma, args$feat_q, na.rm = TRUE))
    c3 <- as.numeric(stats::quantile(d$x3_zero,  args$feat_q, na.rm = TRUE))
    d$pass_segIdx <- d$segIdx   > seg_cut
    d$pass_x1     <- d$x1_rho   > c1
    d$pass_x2     <- d$x2_sigma > c2
    d$pass_x3     <- d$x3_zero  > c3
    d$selected <- if (USE_GATES)
      d$pass_segIdx & d$pass_x1 & d$pass_x2 & d$pass_x3 else d$pass_segIdx
    r$table <- d
    r$cuts  <- c(segIdx = seg_cut, x1 = c1, x2 = c2, x3 = c3)
    cat(sprintf("  %-28s segIdx>%.4f -> %4d selected of %4d",
                ct, seg_cut, sum(d$selected), nrow(d)))
    if (!USE_GATES)
      cat(sprintf("   [would be %d with the 3 gates]",
                  sum(d$pass_segIdx & d$pass_x1 & d$pass_x2 & d$pass_x3)))
    cat("\n")

    g <- sort(d$gene[d$selected])
    if (args$universe == "genomewide") {
      n_pre <- length(g); g <- sort(intersect(g, meta$panel_genes))
      cat(sprintf("  %s: genome-wide selection intersected with panel: %d -> %d\n",
                  ct, n_pre, length(g)))
    }
    sets[[ct]] <- g; info[[ct]] <- r
    write.table(r$table, file.path(args$out_dir,
                sprintf("seg_index_%s_%s.tsv", args$universe, ct_tag(ct))),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  out_qs <- file.path(args$out_dir, sprintf("Reference_SEG_%s.qs", args$universe))
  qs::qsave(sets, out_qs)
  chk <- qs::qread(out_qs)
  stopifnot("round-trip failed" = identical(names(chk), names(sets)),
            "round-trip failed: not character" = all(vapply(chk, is.character, logical(1))))
  cat("\nWrote:", out_qs, "\n")
  cat("Set sizes:", paste(sprintf("%s=%d", names(sets), vapply(sets, length, 1L)),
                          collapse = "  "), "\n")

  # pairwise overlap between celltype sets - how celltype-specific are they?
  cts <- names(sets)
  cat("\nPairwise Jaccard between celltype SEG sets:\n")
  jac <- character(0)
  for (i in seq_along(cts)) for (j in seq_along(cts)) if (i < j) {
    a <- sets[[cts[i]]]; b <- sets[[cts[j]]]
    j2 <- if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA
    line <- sprintf("  %-28s vs %-28s shared=%4d J=%.3f", cts[i], cts[j],
                    length(intersect(a, b)), j2)
    cat(line, "\n"); jac <- c(jac, line)
  }
  core <- Reduce(intersect, sets)
  cat("  genes selected in ALL", length(cts), "celltypes:", length(core), "\n")

  prov <- c(
    "Reference-derived SEG sets - provenance",
    "",
    sprintf("Generated : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Script    : R/prepare_seg_geneset_reference.R (stages extract -> fit -> finalise)",
    sprintf("Reference : %s", meta$ref_path),
    sprintf("scMerge   : %s", meta$scmerge_version),
    sprintf("R         : %s", R.version.string),
    sprintf("Universe  : %s", args$universe),
    "",
    "Method    : scMerge::scSEGIndex (Lin et al. 2019 GigaScience 8(9):giz106), fitted",
    "            per CosMx celltype with cell_type = NULL on logNormCounts. Size factors",
    "            are full-transcriptome library sizes computed BEFORE gene subsetting.",
    sprintf("            Dropout filter: genes with >%.0f%% zero cells dropped by scSEGIndex.",
            100 * meta$zero_max),
    RULE_LINES,
    "            x1 = rank(rho)/(n+1); x2 = 1 - rank(sigma)/(n+1); x3 = 1 - rank(zero)/(n+1)",
    "",
    "Absolute stability features in seg_index_<universe>_<CT>.tsv:",
    "  lambda   = rho        mixing proportion of the normal (expressed) component",
    "  sigma / sigma_sq      SD / variance of the normal component",
    "  omega / omega_st      raw dropout rate / mean-scaled dropout (omega * mu_scaled)",
    "  mu / mu_scaled        normal-component mean / min-max scaled mean",
    "  F_donor / p_donor     between-donor ANOVA F within the celltype: the replacement",
    "                        for scSEGIndex's between-celltype F, which is not applicable",
    "                        within a single celltype. REPORTED ONLY - not in the rule.",
    sprintf("                        Reference donors: %d", meta$n_donors),
    "",
    "Label mapping (reference subcluster -> CosMx celltype);",
    "provenance R/AUCell_Neuron.R:29-30 then R/rename_celltypes.R :")
  for (ct in names(CT_MAP))
    prov <- c(prov, sprintf("  %-28s <- %s", ct, paste(CT_MAP[[ct]], collapse = " + ")))
  prov <- c(prov, "", "Set sizes, cuts and fit cost:")
  for (ct in names(info)) {
    r <- info[[ct]]
    prov <- c(prov, sprintf("  %-28s n=%4d of %4d surviving | segIdx>%.4f | %.1f min",
                            ct, length(sets[[ct]]), r$n_surv,
                            r$cuts[["segIdx"]], r$seconds / 60))
  }
  prov <- c(prov, "", "Pairwise Jaccard between celltype sets:", jac,
            sprintf("  selected in all %d celltypes: %d genes", length(cts), length(core)),
            "",
            "NOTE: the `genomewide` sensitivity universe is not run by default. Rerun",
            "with --universe genomewide to add it.",
            "",
            sprintf("Consumed by R/modulescore_seg_control_vs_phf1_distance.R (--seg_qs %s)",
                    out_qs))
  writeLines(prov, file.path(args$out_dir, "Reference_SEG_provenance.txt"))
  cat("Wrote:", file.path(args$out_dir, "Reference_SEG_provenance.txt"), "\n")
  cat("\nsessionInfo():\n"); print(sessionInfo())
}

cat("\nDone (stage:", args$stage, ").\n")
