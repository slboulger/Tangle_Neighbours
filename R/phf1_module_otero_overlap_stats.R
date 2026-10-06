#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# phf1_module_otero_overlap_stats.R
#
# Produces: Table S9
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# phf1_module_otero_overlap_stats.R
#
# STATS ONLY (no figure). Tests whether the PHF1 module of Figure 4B / Table S9
# overlaps the Otero-Garcia layer 2-3 AT8-up signature beyond expression-matched
# expectation. Otero UP only (no DOWN contrast).
#
# SETS -- exactly the two scored in Figure 4B (run_phf1_module_distance_modelp.sh):
#   PHF1 module : phf1_markers/<CT>/<CT>_phf1_geneset.txt (141 symbols)
#   Otero UP    : otero_signatures/otero_at8_signatures.rds $otero_L23_up   (Ex1 u Ex2)
#   Otero DOWN  : $otero_L23_down -- NOT tested; only excluded from the draw pool, so the
#                 null is drawn from genes in neither Otero set
# Secondary row: the phf1cat N53 module cut from the pseudobulk categorical DEG
# (genesets/cosmx_deg_modules/). phf1_markers/ is the FindMarkers (Wilcoxon) output that
# defines the PHF1 module; the secondary row is reported so the result can be checked
# not to depend on which module is used.
#
# UNIVERSE: the probes tested in the canonical categorical DEG (deg/de_cat/de_PHF1/<CT>/),
# because matching needs a per-gene expression level (AveExpr) on the same cells.
# Symbols are mapped to probes with expand_probe_names() (CosMx packs genes, "FKBP1A/C").
#
# TESTS:
#   PRIMARY  -- overlap vs B random sets of the Otero set's size, drawn from panel
#               genes in neither Otero set, matched on AveExpr in 20 quantile bins.
#               fold = observed / mean matched overlap; interval = observed over the
#               97.5/2.5% quantiles of the draws; emp_p two-sided.
#   CMH      -- common OR across the same AveExpr strata.
#   Fisher   -- unadjusted, ONLY to size the expression inflation.
#
# OUTPUT, under plots/phf1_module_otero_overlap/:
#   stats_phf1_module_otero_overlap.txt  -- log, shared genes, sessionInfo
#   stats_phf1_module_otero_overlap.tsv  -- one row per module
#
# Runs from <PROJECT_ROOT>/phf1_v2 in a few minutes:
#   Rscript R/phf1_module_otero_overlap_stats.R

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

HPC_DIR <- "<PROJECT_ROOT>/phf1_v2"
if (dir.exists(HPC_DIR)) setwd(HPC_DIR)
source("R/cosmx_deg_module_utils.R")  # expand_probe_names()

set.seed(42)

CT        <- "Exc-IT-L2-3-CBLN2-HOPX"
DEG_FILE  <- file.path("deg/de_cat/de_PHF1", CT, paste0(CT, "_PHF1TRUEVsPHF1FALSE.tsv"))
PHF1_SET  <- file.path("phf1_markers", CT, paste0(CT, "_phf1_geneset.txt"))
CAT_SET   <- file.path("genesets/cosmx_deg_modules", CT, paste0(CT, "_phf1cat_N53_geneset.txt"))
OTERO_RDS <- "otero_signatures/otero_at8_signatures.rds"
N_BINS    <- 20
B         <- 99999
OUT_DIR   <- "plots/phf1_module_otero_overlap"
NAME      <- "phf1_module_otero_overlap"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_set <- function(f) { g <- trimws(readLines(f)); unique(g[nzchar(g)]) }

deg <- read.delim(DEG_FILE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(deg$gene))
U   <- deg$gene
bin <- cut(rank(deg$AveExpr, ties.method = "first"), breaks = N_BINS, labels = FALSE)
probe_members <- setNames(lapply(U, expand_probe_names), U)
# A probe matches if a member symbol OR the packed probe name itself is listed
# (phf1_markers/ carries "MZT2A/B" verbatim; Otero carries single symbols).
on_panel <- function(sym) U %in% sym | vapply(probe_members, function(m) any(m %in% sym), logical(1))

otero <- readRDS(OTERO_RDS)
stopifnot(all(c("otero_L23_up", "otero_L23_down") %in% names(otero)))
modules <- list("PHF1 module (Fig. 4B, phf1_markers)" = read_set(PHF1_SET),
                "phf1cat N53 (canonical DEG module)"   = read_set(CAT_SET))
o_sets  <- list("Otero L2-3 UP" = otero$otero_L23_up, "Otero L2-3 DOWN" = otero$otero_L23_down)
o_flag  <- lapply(o_sets, on_panel)
if (any(o_flag[[1]] & o_flag[[2]])) stop("a probe is in both Otero UP and DOWN")
pool <- which(!(o_flag[[1]] | o_flag[[2]]))
pool_by_bin <- split(pool, factor(bin[pool], levels = 1:N_BINS))

matched_null <- function(a, b) {
  need <- as.integer(table(factor(bin[b], levels = 1:N_BINS)))
  if (any(need > lengths(pool_by_bin))) stop("matching infeasible; lower N_BINS")
  vapply(seq_len(B), function(i) {
    d <- unlist(mapply(function(p, n) if (n) p[sample.int(length(p), n)] else integer(0),
                       pool_by_bin, need, SIMPLIFY = FALSE))
    sum(a[d])
  }, numeric(1))
}
cmh <- function(a, b) {
  arr <- table(factor(a, c(FALSE, TRUE)), factor(b, c(FALSE, TRUE)), bin)
  keep <- apply(arr, 3, function(m) all(rowSums(m) > 0) && all(colSums(m) > 0))
  arr <- arr[, , keep, drop = FALSE]
  if (dim(arr)[3] < 2) return(c(NA, NA, NA, NA, dim(arr)[3]))
  mt <- mantelhaen.test(arr, exact = FALSE)
  c(unname(mt$estimate), mt$conf.int, mt$p.value, dim(arr)[3])
}

sink(file.path(OUT_DIR, paste0("stats_", NAME, ".txt")), split = TRUE)
cat("=== Overlap of the Figure 4B PHF1 module with the Otero-Garcia L2-3 AT8 signature ===\n")
cat("Run:", format(Sys.time()), "\n\n")
cat("Universe : probes tested in", DEG_FILE, "(n =", length(U), ")\n")
cat("Otero    :", OTERO_RDS, "\n")
for (k in names(o_sets))
  cat(sprintf("  %-16s %d symbols, %d on the universe\n", k, length(o_sets[[k]]), sum(o_flag[[k]])))
cat("Matching :", N_BINS, "AveExpr quantile bins; B =", B,
    "draws from the", length(pool), "probes in neither Otero set\n\n")

rows <- list()
for (m in names(modules)) {
  sym <- modules[[m]]
  a <- on_panel(sym)
  lost <- setdiff(sym, c(U[a], unlist(probe_members[a])))
  cat("------------------------------------------------------------------\n")
  cat(m, ":", length(sym), "symbols ->", sum(a), "probes on the universe",
      if (length(lost)) paste0("(not on universe: ", paste(lost, collapse = ", "), ")") else "", "\n")
  for (k in "Otero L2-3 UP") {
    b <- o_flag[[k]]
    k_obs <- sum(a & b)
    nul <- matched_null(a, which(b)); mu <- mean(nul)
    ft <- fisher.test(table(factor(a, c(FALSE, TRUE)), factor(b, c(FALSE, TRUE))))
    cm <- cmh(a, b)
    r <- tibble(module = m, otero_set = k, n_module = sum(a), n_otero = sum(b),
                n_universe = length(U), overlap = k_obs,
                pct_of_module = 100 * k_obs / sum(a),
                expected_random = sum(a) * sum(b) / length(U),
                expected_matched = mu,
                fold_matched = k_obs / mu,
                fold_lcl = k_obs / quantile(nul, 0.975, names = FALSE),
                fold_ucl = k_obs / quantile(nul, 0.025, names = FALSE),
                emp_p = (1 + sum(abs(nul - mu) >= abs(k_obs - mu))) / (B + 1),
                n_draws_ge_obs = sum(nul >= k_obs), B = B,
                cmh_or = cm[1], cmh_lcl = cm[2], cmh_ucl = cm[3], cmh_p = cm[4],
                cmh_strata = cm[5],
                fisher_or = unname(ft$estimate), fisher_lcl = ft$conf.int[1],
                fisher_ucl = ft$conf.int[2], fisher_p = ft$p.value)
    rows[[length(rows) + 1]] <- r
    cat(sprintf("\n  vs %s: shared %d (%.0f%% of module); expected %.1f random, %.1f matched\n",
                k, k_obs, r$pct_of_module, r$expected_random, mu))
    cat(sprintf("    PRIMARY  fold over matched = %.2f [%.2f, %.2f]; emp p = %.2g (%d of %d draws >= observed)\n",
                r$fold_matched, r$fold_lcl, r$fold_ucl, r$emp_p, r$n_draws_ge_obs, B))
    cat(sprintf("    CMH      OR = %.2f [%.2f, %.2f]; p = %.2g (%d strata)\n",
                r$cmh_or, r$cmh_lcl, r$cmh_ucl, r$cmh_p, as.integer(r$cmh_strata)))
    cat(sprintf("    Fisher   OR = %.2f [%.2f, %.2f]; p = %.2g  (unadjusted; sizes the expression inflation)\n",
                r$fisher_or, r$fisher_lcl, r$fisher_ucl, r$fisher_p))
    cat("    shared genes:\n")
    cat(strwrap(paste(sort(U[a & b]), collapse = ", "), width = 96, prefix = "      "), sep = "\n")
  }
  cat("\n")
}
res <- bind_rows(rows)
write.table(res, file.path(OUT_DIR, paste0("stats_", NAME, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("------------------------------------------------------------------\n")
cat("Reading notes\n")
cat(" - Quote fold_matched with its interval. Otero DOWN is not tested (pool exclusion only).\n")
cat(" - The fold interval is the spread of expression-matched random sets (gene-set-choice\n")
cat("   variability); it does not carry donor uncertainty in the module itself.\n")
cat(" - The PHF1 module and Otero UP share the genes listed above.\n")
cat(" - Universe = the DEG-tested probes (needed for AveExpr); Figure 4B scores on the full\n")
cat("   panel, so per-set panel counts here can differ slightly from that figure's log.\n\n")
print(sessionInfo())
sink()
