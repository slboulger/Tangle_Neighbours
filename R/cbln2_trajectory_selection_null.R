#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cbln2_trajectory_selection_null.R
#
# Figure panels: 2B
# Within-donor label permutation with the full trajectory selection repeated on every
# permutation; gives the empirical p for tangle-bearing neurons lying later along pseudotime.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# cbln2_trajectory_selection_null.R
#
# P-value for "tangle % increases along the trajectory", calibrated for how the trajectory
# was chosen. R/cbln2_trajectory_optimise.R picks the root (fewest-PHF1+ cluster) and the
# settings (space x k) that best follow tangle %, so a test on the chosen trajectory with the
# same labels is optimistic. This script re-runs the WHOLE selection on PHF1 labels shuffled
# within donor (each donor keeps its tangle count) and asks how often a random labelling
# does as well as the real one.
#
# HOW
#   * Every slingshot fit is label-free given (space, k, root). The search grid is read from
#     step 1's export (subcluster_phf1/<tag>_traj_selected_params.rds); every k-means clustering
#     is re-derived (seed 100) and slingshot is fitted from EVERY eligible root (>= min_clus_cells),
#     because a shuffled labelling can make a different cluster the fewest-PHF1+ one. Fits are
#     cached in subcluster_phf1/<tag>_traj_allroot_fits.qs.
#   * Check: with the REAL labels the selection must reproduce step 1's choice and objective
#     exactly (stops otherwise).
#   * Each permutation: shuffle PHF1 within donor -> per clustering, root = fewest-PHF1+ cluster ->
#     objective (donor-mean AUC, donors >= min_pos_donor PHF1+) for every space x k -> keep the
#     best (ties -> pooled AUC), exactly as step 1.
#   * Two statistics, both from the permutation-chosen trajectory:
#       auc    donor-mean AUC of PHF1+ later in pseudotime (step 1's objective)
#       slope  donor-mean slope of % PHF1+ across the 10 equal-width pseudotime bins
#              (percentage points per bin) -- the quantity in plot_<tag>_traj_phf1_by_bin_mean.pdf,
#              i.e. "tangle % increases"
#     p = (1 + #{null >= observed}) / (1 + B), one-sided. With B = 1000 the smallest p is 0.001.
#
# OUTPUT (plots/<tag>_trajectory_optimise/)
#   plot_<tag>_trajopt_null.pdf        null distributions with the observed values
#   source_data_<tag>_trajopt_null.tsv one row per permutation (+ observed)
#   stats_<tag>_trajopt_null.tsv       observed, null mean / sd / 95th pct, p
#   stats_<tag>_trajopt_null.txt       log (+ sessionInfo)
#   R/plot_cbln2_traj_phf1_by_bin_mean.R prints this p on its figure when the tsv exists.
#
# RUN after steps 1-3:
#   Rscript R/cbln2_trajectory_selection_null.R --celltype <ct> --n_cores <cores> [--n_perm 1000]
# First run fits ~90 trajectories (harmony space is the slow part); later runs reuse them.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(argparse)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/cbln2_embedding_utils.R")
source("R/cbln2_trajectory_utils.R")

parser <- ArgumentParser(description = "Selection-adjusted permutation p for tangle % along the trajectory")
parser$add_argument("--celltype", default = CBLN2_CT)
parser$add_argument("--n_perm", type = "integer", default = 1000L)
parser$add_argument("--n_bins", type = "integer", default = 10L)
parser$add_argument("--n_cores", type = "integer", default = max(1L, min(8L, parallel::detectCores() - 1L)))
parser$add_argument("--seed", type = "integer", default = 1L)
parser$add_argument("--refit", action = "store_true", default = FALSE)
args <- parser$parse_args()

CT  <- args$celltype
TAG <- ct_tag(CT)
OUT_DIR     <- file.path("plots", paste0(TAG, "_trajectory_optimise"))
PFX         <- paste0(TAG, "_trajopt_null")
PARAMS_FILE <- file.path(EMB_CACHE_DIR, paste0(TAG, "_traj_selected_params.rds"))
GRID_FILE   <- file.path(OUT_DIR, paste0("stats_", TAG, "_trajopt_grid.tsv"))
FITS_CACHE  <- file.path(EMB_CACHE_DIR, paste0(TAG, "_traj_allroot_fits.qs"))
PHF1_COLOUR <- "#D55E00"
if (!file.exists(PARAMS_FILE)) stop("Run R/cbln2_trajectory_optimise.R --celltype '", CT, "' first (", PARAMS_FILE, ").")
sel <- readRDS(PARAMS_FILE)
stopifnot(identical(sel$celltype %||% CBLN2_CT, CT))
SPACES   <- sel$grid$spaces
K_GRID   <- sel$grid$k_grid
MIN_CLUS <- sel$grid$min_clus_cells %||% 100L
MIN_POS  <- sel$min_pos_donor %||% 10L   # fallback if the export lacks min_pos_donor (CBLN2: same donors at 10 or 5)

## 1. Inputs -----------------------------------------------------------------------------
sce <- qs::qread(ct_sce_path(CT))
cd  <- as.data.frame(colData(sce)); cd$cell_id <- colnames(sce)
cd$PHF1 <- as.logical(as.character(cd$PHF1))
rm(sce); invisible(gc())
emb <- load_or_build_cbln2_embedding(function() stop("Embedding cache missing for ", CT), need_full = TRUE, ct = CT)
m <- match(cd$cell_id, emb$coords$cell_id); stopifnot(!anyNA(m))
rd_list <- list(
  umap    = { u <- as.matrix(emb$coords[m, c("UMAP_1", "UMAP_2")]); rownames(u) <- cd$cell_id; u },
  harmony = emb$harmony[cd$cell_id, , drop = FALSE])
pos   <- cd$PHF1
donor <- as.character(cd$sample_id)
donor_idx  <- split(seq_along(pos), donor)
n_pos_d    <- vapply(donor_idx, function(i) sum(pos[i]), numeric(1))
obj_donors <- names(n_pos_d)[n_pos_d >= MIN_POS]

## 2. Clusterings + label-free fits from EVERY eligible root (cached) -------------------
clusterings <- list()
for (sp in SPACES) for (k in K_GRID) clusterings[[paste(sp, k)]] <- kmeans_clusters(rd_list[[sp]], k)
configs <- list()
for (key in names(clusterings)) {
  n_cl <- table(clusterings[[key]])
  for (r in names(n_cl)[n_cl >= MIN_CLUS])
    configs[[length(configs) + 1]] <- list(key = key, space = strsplit(key, " ")[[1]][1],
                                           k = as.integer(strsplit(key, " ")[[1]][2]), root = r)
}
spec <- list(spaces = SPACES, k_grid = K_GRID, min_clus = MIN_CLUS, seed_kmeans = SEED_KMEANS,
             approx_points = APPROX_POINTS, cells = cd$cell_id, use_dims = emb$prep$use_dims,
             umap_checksum = sum(rd_list$umap), harmony_checksum = sum(rd_list$harmony))
fits <- NULL
if (!args$refit && file.exists(FITS_CACHE)) {
  cc <- qs::qread(FITS_CACHE)
  if (identical(cc$spec, spec)) { fits <- cc$fits; message("Reusing cached all-root fits: ", FITS_CACHE) }
}
if (is.null(fits)) {
  message(sprintf("Fitting %d label-free trajectories (every eligible root) on %d cores ...", length(configs), args$n_cores))
  t0 <- Sys.time()
  pts <- parallel::mclapply(configs, function(cf) {
    sds <- suppressMessages(run_slingshot(rd_list[[cf$space]], clusterings[[cf$key]], cf$root, cf$key))
    if (is.null(sds)) NULL else as.numeric(avg_pt(sds))
  }, mc.cores = args$n_cores, mc.preschedule = FALSE)
  fits <- Map(function(cf, p) c(cf, list(pt = p)), configs, pts)
  message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  qs::qsave(list(spec = spec, fits = fits), FITS_CACHE)
}
ok <- vapply(fits, function(f) !is.null(f$pt) && !anyNA(f$pt), logical(1))
fit_id <- vapply(fits, function(f) paste(f$key, f$root), "")
names(fits) <- fit_id

## Pre-compute, per fit: within-donor ranks (objective donors), global ranks (tie-break), and
## equal-width bins (slope statistic) -- a permutation then costs only sums.
prep <- lapply(fits[ok], function(f) {
  list(rw = lapply(donor_idx[obj_donors], function(i) rank(f$pt[i])),
       rg = rank(f$pt),
       bin = cut(f$pt, args$n_bins, labels = FALSE))
})
auc_r <- function(r, lab) { n1 <- sum(lab); n0 <- length(lab) - n1; (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
score <- function(fid, lab) {             # donor-mean AUC + pooled AUC for one fit
  p <- prep[[fid]]
  c(obj = mean(vapply(obj_donors, function(d) auc_r(p$rw[[d]], lab[donor_idx[[d]]]), numeric(1))),
    pooled = auc_r(p$rg, lab))
}
slope_stat <- function(fid, lab) {        # donor-mean slope of % PHF1+ across bins (pp per bin)
  b <- prep[[fid]]$bin
  mean(vapply(obj_donors, function(d) {
    i <- donor_idx[[d]]
    pct <- 100 * tapply(lab[i], factor(b[i], levels = seq_len(args$n_bins)), mean)
    x <- seq_len(args$n_bins)[!is.na(pct)]
    unname(coef(lm(pct[!is.na(pct)] ~ x))[2])
  }, numeric(1)))
}
## The selection, exactly as step 1: per clustering the fewest-PHF1+ root, then the best objective.
select_traj <- function(lab) {
  cand <- vapply(names(clusterings), function(key) {
    root <- pick_cluster(cluster_table(clusterings[[key]], lab, donor), "fraction", MIN_CLUS, "min")
    paste(key, root)
  }, "")
  cand <- unname(cand[cand %in% names(prep)])   # unnamed, so the score rows are named by fit id
  sc <- t(vapply(cand, score, numeric(2), lab = lab))
  best <- cand[order(-sc[, "obj"], -sc[, "pooled"])[1]]
  list(fit = best, auc = sc[best, "obj"], slope = slope_stat(best, lab))
}

## 3. Observed (must reproduce step 1) ----------------------------------------------------
obs <- select_traj(pos)
exp_id <- paste(sel$space, sel$k, sel$root)
if (!identical(unname(obs$fit), exp_id) || abs(obs$auc - sel$objective_value) > 1e-9)
  stop(sprintf("Selection on the real labels gives %s (AUC %.6f), but step 1 chose %s (AUC %.6f). Re-run step 1 or --refit.",
               obs$fit, obs$auc, exp_id, sel$objective_value))
message(sprintf("Observed: %s, donor-mean AUC %.4f, slope %+.4f pp/bin -- matches step 1", obs$fit, obs$auc, obs$slope))

## 4. Permutations ------------------------------------------------------------------------
message(sprintf("Selection-adjusted null: %d within-donor permutations ...", args$n_perm))
set.seed(args$seed)
null <- bind_rows(lapply(seq_len(args$n_perm), function(b) {
  lab <- logical(length(pos))
  for (d in names(donor_idx)) { i <- donor_idx[[d]]; lab[i] <- pos[i][sample.int(length(i))] }
  r <- select_traj(lab)
  tibble(perm = b, fit_chosen = r$fit, auc = r$auc, slope = r$slope)
}))
p_of <- function(o, nl) (1 + sum(nl >= o)) / (1 + length(nl))
res <- tibble(
  statistic = c("auc", "slope"),
  description = c("donor-mean AUC, PHF1+ later in pseudotime (step 1 objective)",
                  "donor-mean slope of % PHF1+ across pseudotime bins (pp per bin)"),
  observed = c(obs$auc, obs$slope),
  null_mean = c(mean(null$auc), mean(null$slope)), null_sd = c(sd(null$auc), sd(null$slope)),
  null_q95 = c(quantile(null$auc, 0.95, names = FALSE), quantile(null$slope, 0.95, names = FALSE)),
  p_selection_adjusted = c(p_of(obs$auc, null$auc), p_of(obs$slope, null$slope)),
  n_perm = args$n_perm, n_donors = length(obj_donors), min_pos_donor = MIN_POS,
  chosen_trajectory = obs$fit)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.table(res, file.path(OUT_DIR, paste0("stats_", PFX, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(null %>% mutate(type = "null"),
                      tibble(perm = NA_integer_, fit_chosen = obs$fit, auc = obs$auc, slope = obs$slope, type = "observed")),
            file.path(OUT_DIR, paste0("source_data_", PFX, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

## 5. Plot + log ------------------------------------------------------------------------------
nd <- null %>% pivot_longer(c(auc, slope), names_to = "statistic", values_to = "value") %>%
  mutate(statistic = factor(statistic, levels = c("auc", "slope"),
                            labels = c("Donor-mean AUC (PHF1+ later)", "Slope of % PHF1+ (pp per bin)")))
od <- res %>% mutate(statistic = factor(statistic, levels = c("auc", "slope"),
                                        labels = levels(nd$statistic)),
                     label = sprintf("p = %s", ifelse(p_selection_adjusted <= 1 / (1 + n_perm),
                                                      sprintf("< %.3g", 1 / n_perm), sprintf("%.3g", p_selection_adjusted))))
p <- ggplot(nd, aes(value)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(data = od, aes(xintercept = observed), colour = PHF1_COLOUR, linewidth = 0.5) +
  geom_text(data = od, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3, size = 2.2) +
  facet_wrap(~ statistic, scales = "free_x") +
  labs(x = "Best trajectory under shuffled labels", y = "Permutations") +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), strip.background = element_blank(),
        plot.margin = margin(4, 8, 4, 4))
ggsave(file.path(OUT_DIR, paste0("plot_", PFX, ".pdf")), p, width = 12, height = 5, units = "cm", device = "pdf")

old_w <- options(width = 160)
sink(file.path(OUT_DIR, paste0("stats_", PFX, ".txt")))
cat(CT, ": selection-adjusted permutation test of tangle % along the trajectory\n", strrep("=", 78), "\n", sep = "")
cat("Script: R/cbln2_trajectory_selection_null.R   Run:", format(Sys.time()), "\n")
cat(sprintf("Search re-run per permutation: spaces %s x k %s, root = fewest-PHF1+ cluster (>= %d cells); %d label-free fits, %d failed.\n",
            paste(SPACES, collapse = "/"), paste(K_GRID, collapse = ","), MIN_CLUS, length(fits), sum(!ok)))
cat(sprintf("Objective donors (>= %d PHF1+): %s\n", MIN_POS, paste(obj_donors, collapse = ", ")))
cat(sprintf("PHF1 labels shuffled within donor (each donor keeps its tangle count); B = %d; seed %d.\n\n", args$n_perm, args$seed))
cat("Observed selection reproduces step 1:", obs$fit, "\n\n")
print(as.data.frame(res %>% select(statistic, observed, null_mean, null_sd, null_q95, p_selection_adjusted)),
      digits = 4, row.names = FALSE)
cat("\np = (1 + #{null >= observed}) / (1 + B), one-sided; smallest possible p =", signif(1 / (1 + args$n_perm), 3), "\n")
cat("Trajectory chosen under the null (top 10):\n"); print(head(sort(table(null$fit_chosen), decreasing = TRUE), 10))
cat("\n--- sessionInfo ---\n"); print(sessionInfo())
sink(); options(old_w)
message(sprintf("AUC %.3f (null %.3f, p %.3g); slope %+.3f pp/bin (null %+.3f, p %.3g) -> %s",
                obs$auc, mean(null$auc), res$p_selection_adjusted[1], obs$slope, mean(null$slope),
                res$p_selection_adjusted[2], OUT_DIR))
