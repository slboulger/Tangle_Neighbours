#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cbln2_trajectory_optimise.R
#
# Figure panels: 2B - upstream step
# Trajectory settings search: exports the selected space, k and root that
# cbln2_slingshot_trajectory.R fits.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# cbln2_trajectory_optimise.R
#
# Step 1 of the tangle-associated trajectory (modelled on Foiani et al. 2026, "Construction
# of tau-associated trajectory"): pick the root, then choose the trajectory settings under
# which pseudotime best follows tangle-bearing (PHF1+) %. The chosen settings are exported
# for R/cbln2_slingshot_trajectory.R (step 2); genes are then modelled against pseudotime by
# R/deg_dream_linear_pseudotime_phf1.r (step 3).
#
#   Foiani                                   here
#   early cell = most pTau-negative sample   root = the fewest-PHF1+ k-means cluster (>= 100 cells)
#   palantir on diffusion components         slingshot on the CBLN2 embedding
#   fixed settings                           space x k chosen by the objective below
#
# SEARCH (pre-declared): space = UMAP (2-D figure embedding) | harmony (the 10 elbow-chosen
# PCs); k = 4, 5, 6, 8, 10, 12 k-means clusters (seed 100, nstart 10, iter.max 50). Root
# fixed by the rule within each clustering. Slingshot defaults, approx_points 300,
# curve-weighted mean pseudotime over all lineages. 12 fits.
# OBJECTIVE: mean over donors with >= 5 PHF1+ cells (--min_pos_donor) of the per-donor rank
# AUC (P(pseudotime of a PHF1+ cell > a PHF1- cell from the same donor)). Ties -> pooled AUC.
# For CBLN2 the donors entering the objective are the same at a threshold of 5 or 10.
#
# Because the settings are selected on this objective, R/cbln2_trajectory_selection_null.R
# repeats the whole selection on within-donor permuted PHF1 labels to give the
# selection-adjusted p for tangle-bearing neurons lying later along pseudotime.
#
# OUTPUT  plots/cbln2_trajectory_optimise/  (prefix cbln2_trajopt_)
#   plot_..._landscape.pdf      objective per config (k, by space), chosen marked
#   plot_..._selected_umap.pdf  chosen pseudotime on the figure UMAP, PHF1+ ringed
#   plot_..._phf1_by_bin.pdf    % PHF1+ per pseudotime bin, chosen config
#   source_data_<plot>.tsv, stats_cbln2_trajopt.txt (+ sessionInfo), stats_cbln2_trajopt_grid.tsv
#   subcluster_phf1/cbln2_traj_selected_params.rds   read by cbln2_slingshot_trajectory.R
#   subcluster_phf1/cbln2_traj_grid_fits.qs          cached fits (reused while grid + embedding unchanged)
# RUN (an R with slingshot + harmony), after the embedding cache exists:
#   Rscript R/cbln2_trajectory_optimise.R --n_cores <cores allocated>
# --n_cores should be the allocation, not detectCores() (which sees the whole node).

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(argparse)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")
source("R/cbln2_embedding_utils.R")
source("R/cbln2_trajectory_utils.R")

parser <- ArgumentParser(description = "Choose the trajectory that best follows tangle %")
parser$add_argument("--celltype", default = CBLN2_CT,
                    help = "celltype (default CBLN2; files named by ct_tag(), CBLN2 paths unchanged)")
parser$add_argument("--spaces", default = "umap,harmony")
parser$add_argument("--k_grid", default = "4,5,6,8,10,12")
parser$add_argument("--min_clus_cells", type = "integer", default = 100L)
parser$add_argument("--min_pos_donor", type = "integer", default = 5L,
                    help = "donors with >= this many PHF1+ cells enter the objective (CBLN2's donor set is the same at 5 or 10)")
parser$add_argument("--n_bins", type = "integer", default = 10L)
parser$add_argument("--n_cores", type = "integer", default = max(1L, min(8L, parallel::detectCores() - 1L)))
parser$add_argument("--refit", action = "store_true", default = FALSE, help = "ignore the cached fits")
args <- parser$parse_args()
SPACES <- strsplit(args$spaces, ",")[[1]]
K_GRID <- as.integer(strsplit(args$k_grid, ",")[[1]])
stopifnot(all(SPACES %in% c("umap", "harmony")))

CT          <- args$celltype
TAG         <- ct_tag(CT)
OUT_DIR     <- file.path("plots", paste0(TAG, "_trajectory_optimise"))
PFX         <- paste0(TAG, "_trajopt")
FITS_CACHE  <- file.path(EMB_CACHE_DIR, paste0(TAG, "_traj_grid_fits.qs"))
PARAMS_FILE <- file.path(EMB_CACHE_DIR, paste0(TAG, "_traj_selected_params.rds"))
PHF1_COLOUR <- "#D55E00"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_tsv_ <- function(x, f) write.table(x, file.path(OUT_DIR, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
save_pdf <- function(p, name, w, h) ggsave(file.path(OUT_DIR, paste0("plot_", PFX, "_", name, ".pdf")), p,
                                           width = w, height = h, units = "cm", device = "pdf")
theme_xy <- function() theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1), plot.margin = margin(4, 4, 4, 4))

## 1. Inputs ------------------------------------------------------------------------------
sce <- qs::qread(ct_sce_path(CT))
cd  <- as.data.frame(colData(sce)); cd$cell_id <- colnames(sce)
cd$PHF1 <- as.logical(as.character(cd$PHF1))
rm(sce); invisible(gc())
emb <- load_or_build_cbln2_embedding(function() stop("Embedding cache missing: run R/plot_cbln2_umap_phf1_p95.R --celltype '", CT, "' first."),
                                     need_full = TRUE, ct = CT)
if (sum(tapply(cd$PHF1, cd$sample_id, sum) >= args$min_pos_donor) < 3)
  stop("Fewer than 3 donors with >= ", args$min_pos_donor, " PHF1+ cells in ", CT,
       " -- the objective needs >= 3; lower --min_pos_donor (and say so).")
if (!requireNamespace("slingshot", quietly = TRUE)) stop("slingshot is not installed in this R.")
m <- match(cd$cell_id, emb$coords$cell_id); stopifnot(!anyNA(m))
rd_list <- list(
  umap    = { u <- as.matrix(emb$coords[m, c("UMAP_1", "UMAP_2")]); rownames(u) <- cd$cell_id; u },
  harmony = emb$harmony[cd$cell_id, , drop = FALSE])
pos   <- cd$PHF1
donor <- as.character(cd$sample_id)
n_pos_d    <- tapply(pos, donor, sum)
obj_donors <- names(n_pos_d)[n_pos_d >= args$min_pos_donor]

## 2. Fits: root = fewest-PHF1+ cluster in each clustering (cached) -----------------------
grid_spec <- list(spaces = SPACES, k_grid = K_GRID, min_clus_cells = args$min_clus_cells,
                  seed_kmeans = SEED_KMEANS, approx_points = APPROX_POINTS, cells = cd$cell_id,
                  use_dims = emb$prep$use_dims, umap_checksum = sum(rd_list$umap),
                  harmony_checksum = sum(rd_list$harmony), root_rule = "fraction",
                  phf1_checksum = sum(which(pos)))
fits <- NULL
if (!args$refit && file.exists(FITS_CACHE)) {
  cc <- qs::qread(FITS_CACHE)
  if (identical(cc$spec, grid_spec)) { fits <- cc$fits; message("Reusing cached fits: ", FITS_CACHE) }
}
if (is.null(fits)) {
  configs <- list()
  for (sp in SPACES) for (k in K_GRID) {
    cl   <- kmeans_clusters(rd_list[[sp]], k)
    root <- pick_cluster(cluster_table(cl, pos, donor), "fraction", args$min_clus_cells, "min")
    configs[[length(configs) + 1]] <- list(space = sp, k = k, root = root, cluster = unname(cl))
  }
  message(sprintf("Fitting %d configurations on %d cores ...", length(configs), args$n_cores))
  t0 <- Sys.time()
  fits <- parallel::mclapply(configs, function(cf) {
    sds <- suppressMessages(run_slingshot(rd_list[[cf$space]], setNames(cf$cluster, cd$cell_id),
                                          cf$root, paste(cf$space, cf$k)))
    if (is.null(sds)) return(c(cf, list(ok = FALSE)))
    c(cf, list(ok = TRUE, pt = as.numeric(avg_pt(sds)),
               lineages = paste(map_chr(slingshot::slingLineages(sds), paste, collapse = ">"), collapse = " | "),
               n_lineages = length(slingshot::slingLineages(sds)), fit_note = attr(sds, "fit_note"),
               curves = if (cf$space == "umap") lapply(slingshot::slingCurves(sds), function(cv) cv$s[cv$ord, , drop = FALSE])))
  }, mc.cores = args$n_cores, mc.preschedule = FALSE)
  message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  qs::qsave(list(spec = grid_spec, fits = fits), FITS_CACHE)
}
ok <- vapply(fits, function(f) isTRUE(f$ok) && !anyNA(f$pt), logical(1))
n_fail <- sum(!ok); fits <- fits[ok]
cfg_id <- vapply(fits, function(f) sprintf("%s_k%d_root%s", f$space, f$k, f$root), "")

## 3. Objective and choice ----------------------------------------------------------------
auc_mat <- vapply(fits, function(f) vapply(obj_donors, function(d) {
  i <- donor == d; auc_fn(f$pt[i], pos[i]) }, numeric(1)), numeric(length(obj_donors)))
auc_mat <- matrix(auc_mat, nrow = length(obj_donors), dimnames = list(obj_donors, cfg_id))
grid_tab <- tibble(config = cfg_id, space = map_chr(fits, "space"), k = map_int(fits, "k"),
                   root = map_chr(fits, "root"),
                   root_frac_phf1 = map_dbl(fits, ~ mean(pos[.x$cluster == .x$root])),
                   n_lineages = map_int(fits, "n_lineages"), lineages = map_chr(fits, "lineages"),
                   fit_note = map_chr(fits, "fit_note"),
                   donor_mean_auc = colMeans(auc_mat),
                   n_donors_auc_gt_0.5 = colSums(auc_mat > 0.5),
                   pooled_auc = map_dbl(fits, ~ auc_fn(.x$pt, pos)),
                   monotonicity_rho = map_dbl(fits, function(f) {
                     b <- cut(f$pt, args$n_bins, labels = FALSE); pr <- tapply(pos, b, mean)
                     suppressWarnings(cor(as.numeric(names(pr)), pr, method = "spearman")) })) %>%
  bind_cols(as_tibble(t(auc_mat), .name_repair = ~ paste0("auc_", obj_donors)))
best <- order(-grid_tab$donor_mean_auc, -grid_tab$pooled_auc)[1]
grid_tab$selected <- seq_len(nrow(grid_tab)) == best
write_tsv_(grid_tab, paste0("stats_", PFX, "_grid.tsv"))

sel <- fits[[best]]
saveRDS(list(celltype = CT, min_pos_donor = args$min_pos_donor,
             space = sel$space, k = sel$k, root = sel$root, use_dims = emb$prep$use_dims,
             seed_kmeans = SEED_KMEANS, approx_points = APPROX_POINTS,
             root_frac_phf1 = grid_tab$root_frac_phf1[best], cluster_sizes = c(table(sel$cluster)),
             objective = sprintf("mean per-donor AUC, donors >= %d PHF1+", args$min_pos_donor), objective_value = grid_tab$donor_mean_auc[best],
             grid = list(spaces = SPACES, k_grid = K_GRID, min_clus_cells = args$min_clus_cells, n_configs = length(fits)),
             embedding_cache = emb$prep$cache, created = format(Sys.time()),
             script = "R/cbln2_trajectory_optimise.R"), PARAMS_FILE)

## 4. Plots -------------------------------------------------------------------------------
p_land <- ggplot(grid_tab, aes(factor(k), donor_mean_auc, group = space, colour = space)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.4) + geom_point(size = 1.3) +
  geom_point(data = grid_tab[grid_tab$selected, ], shape = 21, size = 3, stroke = 0.6, colour = PHF1_COLOUR) +
  scale_colour_manual(values = c(umap = "#0072B2", harmony = "#009E73"), name = "Fitting space") +
  labs(x = "k (k-means clusters)", y = "Donor-mean AUC (PHF1+ later)") + theme_xy()
save_pdf(p_land, "landscape", 9, 6)
write_tsv_(grid_tab %>% select(config, space, k, root, root_frac_phf1, n_lineages, donor_mean_auc,
                               pooled_auc, monotonicity_rho, selected),
           paste0("source_data_", PFX, "_landscape.tsv"))

df_sel <- tibble(cell_id = cd$cell_id, sample_id = donor, PHF1 = pos,
                 UMAP_1 = rd_list$umap[, 1], UMAP_2 = rd_list$umap[, 2], pt = sel$pt)
p_umap <- ggplot(df_sel %>% arrange(pt), aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = pt), size = 0.2, stroke = 0) +
  { if (!is.null(sel$curves)) lapply(sel$curves, function(s)
      geom_path(data = as.data.frame(s) %>% setNames(c("UMAP_1", "UMAP_2")), colour = "black", linewidth = 0.45)) } +
  geom_point(data = df_sel[df_sel$PHF1, ], shape = 21, colour = PHF1_COLOUR, fill = NA, size = 1.1, stroke = 0.3) +
  scale_colour_viridis_c(option = "cividis", name = "Pseudotime") +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_xy()
save_pdf(p_umap, "selected_umap", 8.5, 6.5)
write_tsv_(df_sel, paste0("source_data_", PFX, "_selected_umap.tsv"))

bins_df <- df_sel %>% mutate(bin = cut(pt, args$n_bins, labels = FALSE)) %>%
  group_by(bin) %>% summarise(n = n(), n_phf1 = sum(PHF1), pct_phf1 = 100 * mean(PHF1), .groups = "drop")
p_bin <- ggplot(bins_df, aes(bin, pct_phf1)) + geom_line(linewidth = 0.5) + geom_point(size = 1) +
  scale_x_continuous(breaks = seq_len(args$n_bins)) +
  labs(x = "Pseudotime bin (equal width)", y = "% tangle-bearing") + theme_xy()
save_pdf(p_bin, "phf1_by_bin", 7.5, 5)
write_tsv_(bins_df, paste0("source_data_", PFX, "_phf1_by_bin.tsv"))

## 5. Log ---------------------------------------------------------------------------------
old_w <- options(width = 180)
sink(file.path(OUT_DIR, paste0("stats_", PFX, ".txt")))
cat(CT, " trajectory: choose settings that best follow tangle % (R/cbln2_trajectory_optimise.R)\n", strrep("=", 78), "\n", sep = "")
cat("Run:", format(Sys.time()), " slingshot", as.character(packageVersion("slingshot")), "\n")
cat("Args:", paste(names(args), unlist(lapply(args, paste, collapse = ",")), sep = "=", collapse = "  "), "\n")
cat(sprintf("Cells %d, PHF1+ %d. Objective donors (>= %d PHF1+): %s\n", length(pos), sum(pos),
            args$min_pos_donor, paste(obj_donors, collapse = ", ")))
cat(sprintf("Root = fewest-PHF1+ cluster (>= %d cells) in each clustering. %d configurations fitted, %d failed.\n\n",
            args$min_clus_cells, length(fits), n_fail))
print(as.data.frame(grid_tab %>% arrange(desc(donor_mean_auc)) %>%
                      select(config, root_frac_phf1, n_lineages, lineages, donor_mean_auc, n_donors_auc_gt_0.5,
                             pooled_auc, monotonicity_rho, selected)), digits = 3, row.names = FALSE)
cat("\nChosen:", cfg_id[best], " -> exported", PARAMS_FILE, "\n")
cat("\n--- sessionInfo ---\n"); print(sessionInfo())
sink(); options(old_w)
message(sprintf("Chosen %s: donor-mean AUC %.3f (%d/%d donors > 0.5) -> %s", cfg_id[best],
                grid_tab$donor_mean_auc[best], grid_tab$n_donors_auc_gt_0.5[best], length(obj_donors), PARAMS_FILE))
