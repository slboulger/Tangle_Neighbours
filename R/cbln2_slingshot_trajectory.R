#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cbln2_slingshot_trajectory.R
#
# Figure panels: 2B-2D - upstream step
# Fits the selected trajectory and writes the per-cell pseudotime table read by 2B-2D,
# 3C and the Table S4 model.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# cbln2_slingshot_trajectory.R
#
# Step 2 of the tangle-associated trajectory through Exc-IT-L2-3-CBLN2-HOPX (modelled on
# Foiani et al. 2026): fit the trajectory chosen by R/cbln2_trajectory_optimise.R (step 1),
# draw it, describe how tangle-bearing (PHF1+) cells sit along it, and write per-cell
# pseudotime for the gene model R/deg_dream_linear_pseudotime_phf1.r (step 3).
#
# METHOD: k-means labels on the fitting space as MST vertices (seed 100, nstart 10,
# iter.max 50), slingshot defaults with approx_points 300, the run_slingshot() fallback
# chain, curve-weighted mean pseudotime over all lineages, a sensitivity fit in the other
# space, a reversal fit from the opposite end, and composition across 10 equal-width
# pseudotime bins.
#
# SETTINGS: space, k and root come from subcluster_phf1/cbln2_traj_selected_params.rds
# (step 1); the k-means labels are re-derived and asserted to match the export. Without an
# export (or with --params none) the default rule is used: UMAP, k = 4, fewest-PHF1+ root.
#
# EMBEDDING: exactly the UMAP of R/plot_cbln2_umap_phf1_p95.R via R/cbln2_embedding_utils.R
# (SCT, harmony on sample_id, nPC from the elbow of the PC spectrum by the project's
# min(co1, co2) rule = 10; the elbow plot below shows how).
#
# The AUCs here are descriptive; the selection-adjusted p for tangle-bearing neurons lying
# later along pseudotime comes from R/cbln2_trajectory_selection_null.R. Pseudotime orders
# transcriptional states, not time.
#
# OUTPUT  plots/cbln2_slingshot_trajectory/  (prefix cbln2_traj_)
#   plot_..._elbow.pdf             PCA spectrum, co1 / co2 / chosen nPC
#   plot_..._umap_clusters.pdf     k-means clusters, % PHF1+, root / end, trajectory
#   plot_..._umap_pseudotime.pdf   pseudotime (cividis), trajectory, PHF1+ ringed
#   plot_..._umap_pseudotime_interactive.html  the same in plotly (hover: cell, donor, pt, ...)
#   plot_..._phf1_by_bin.pdf       % PHF1+ per pseudotime bin, pooled + per donor
#   plot_..._braak_by_bin.pdf      Braak composition per bin
#   plot_..._donor_auc.pdf         per-donor AUC (PHF1+ later) with t CI
#   source_data_<plot>.tsv; source_data_cbln2_traj_percell.tsv (READ BY STEP 3)
#   stats_cbln2_traj.txt (+ sessionInfo); stats_cbln2_traj_{clusters,lineages,donor_auc}.tsv
#
# RUN (an R with slingshot), after step 1:
#   Rscript R/cbln2_slingshot_trajectory.R

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(argparse)
})
setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")               # fig_theme, braak_palette, braak_line_palette, braak_levels
source("R/cbln2_embedding_utils.R")  # the shared CBLN2 embedding
source("R/cbln2_trajectory_utils.R") # clustering, slingshot, AUC helpers (shared with step 1)

parser <- ArgumentParser(description = "Fit and describe a celltype's tangle-associated trajectory")
parser$add_argument("--celltype", default = CBLN2_CT,
                    help = "celltype (default CBLN2; files named by ct_tag(), CBLN2 paths unchanged)")
parser$add_argument("--min_pos_donor", type = "integer", default = 5L,
                    help = "donors with >= this many PHF1+ cells summarise the AUC; OVERRIDDEN by step 1's export when it is used")
parser$add_argument("--params", default = "auto",
                    help = "auto = subcluster_phf1/<tag>_traj_selected_params.rds if present; none = UMAP k=4 fewest-PHF1+ root; or a path")
parser$add_argument("--min_clus_cells", type = "integer", default = 100L)
parser$add_argument("--n_bins", type = "integer", default = 10L)
args <- parser$parse_args()

CT            <- args$celltype
TAG           <- ct_tag(CT)
OUTPUT_DIR    <- file.path("plots", paste0(TAG, "_slingshot_trajectory"))
PFX           <- paste0(TAG, "_traj")
stopifnot(file.path(OUTPUT_DIR, paste0("source_data_", PFX, "_percell.tsv")) == ct_percell_file(CT))
PHF1_COLOUR   <- "#D55E00"
CLUSTER_COLOURS <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7", "#000000", "#D55E00",
                     "#882255", "#44AA99", "#999933", "#AA4499", "#F0E442")  # Okabe-Ito + Tol muted, k <= 12
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_tsv_ <- function(x, f) write.table(x, file.path(OUTPUT_DIR, f), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
save_pdf <- function(p, name, w, h) ggsave(file.path(OUTPUT_DIR, paste0("plot_", PFX, "_", name, ".pdf")), p,
                                           width = w, height = h, units = "cm", device = "pdf")
theme_xy <- function() theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1), plot.margin = margin(4, 4, 4, 4))

## 1. Load ----------------------------------------------------------------------------------
sce <- qs::qread(ct_sce_path(CT))
cd  <- as.data.frame(colData(sce)); cd$cell_id <- colnames(sce)
cd$PHF1 <- as.logical(as.character(cd$PHF1))
rm(sce); invisible(gc())
stopifnot(all(as.character(cd$celltype) == CT), !anyNA(cd$PHF1))
emb <- load_or_build_cbln2_embedding(function() stop("Embedding cache missing: run R/plot_cbln2_umap_phf1_p95.R --celltype '", CT, "' first."),
                                     need_full = TRUE, ct = CT)
if (!requireNamespace("slingshot", quietly = TRUE)) stop("slingshot is not installed in this R.")
m <- match(cd$cell_id, emb$coords$cell_id); stopifnot(!anyNA(m))
umap_mat <- as.matrix(emb$coords[m, c("UMAP_1", "UMAP_2")]); rownames(umap_mat) <- cd$cell_id
harm_mat <- emb$harmony[cd$cell_id, , drop = FALSE]
use_dims <- emb$prep$use_dims
pos   <- cd$PHF1
donor <- as.character(cd$sample_id)

## 2. Elbow plot (how nPC was chosen) -------------------------------------------------------
pcs <- pick_pcs(emb$stdev)
stopifnot(identical(as.integer(pcs$min_pc), as.integer(emb$prep$pcs_elbow$min_pc)))
elbow_df <- tibble(PC = seq_along(emb$stdev), stdev = emb$stdev, pct_stdev = pcs$percent_stdv,
                   cumulative_pct = pcs$cumulative, used = seq_along(emb$stdev) %in% use_dims)
cut_df <- tibble(PC = c(pcs$co2, pcs$co1),
                 rule = c(sprintf("knee (co2): last drop > 0.1 -> PC %d", pcs$co2),
                          sprintf("90%% cumulative (co1) -> PC %d", pcs$co1)))
p_elbow <- ggplot(elbow_df, aes(PC, pct_stdev)) +
  geom_line(linewidth = 0.3, colour = "grey40") +
  geom_point(aes(fill = used), shape = 21, size = 1.4, stroke = 0.25) +
  geom_vline(data = cut_df, aes(xintercept = PC, linetype = rule), linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "white"),
                    labels = c(`TRUE` = sprintf("used (1:%d)", max(use_dims)), `FALSE` = "not used"), name = NULL) +
  scale_linetype_manual(values = c("dashed", "dotted"), name = NULL) +
  labs(x = "Principal component", y = "% of total SD") +
  theme_xy() + theme(legend.position = "inside", legend.position.inside = c(0.45, 0.72),
                     legend.background = element_blank(), legend.key.height = grid::unit(7, "pt"))
save_pdf(p_elbow, "elbow", 8, 6)
write_tsv_(elbow_df, paste0("source_data_", PFX, "_elbow.tsv"))

## 3. Settings, clusters, root, end ---------------------------------------------------------
params_file <- switch(args$params, auto = file.path(EMB_CACHE_DIR, paste0(TAG, "_traj_selected_params.rds")),
                      none = NA_character_, args$params)
sel <- if (!is.na(params_file) && file.exists(params_file)) readRDS(params_file) else NULL
if (!is.na(params_file) && args$params != "auto" && is.null(sel)) stop("--params file not found: ", params_file)
if (!is.null(sel) && !identical(sel$celltype %||% CBLN2_CT, CT)) stop("Step 1 export is for ", sel$celltype, ", not ", CT)
MIN_POS_DONOR <- sel$min_pos_donor %||% args$min_pos_donor   # always step 1's threshold when its export is used
fit_space <- sel$space %||% "umap"
k_use     <- sel$k %||% 4L
rd_main   <- if (fit_space == "umap") umap_mat else harm_mat
rd_other  <- if (fit_space == "umap") harm_mat else umap_mat
other_space <- if (fit_space == "umap") sprintf("harmony 1:%d", max(use_dims)) else "umap"

cl <- kmeans_clusters(rd_main, k_use)
clus_tab <- cluster_table(cl, pos, donor)
root <- sel$root %||% pick_cluster(clus_tab, "fraction", args$min_clus_cells, "min")
if (!is.null(sel) && (!identical(c(table(cl)), sel$cluster_sizes) ||
                      abs(clus_tab$frac_phf1[clus_tab$cluster == root] - sel$root_frac_phf1) > 1e-12))
  stop("k-means labels do not reproduce step 1's clustering; re-run R/cbln2_trajectory_optimise.R")
endc <- pick_cluster(clus_tab, "fraction", args$min_clus_cells, "max")
clus_tab <- clus_tab %>% mutate(role = case_when(cluster == root ~ "root", cluster == endc ~ "end (reversal)", TRUE ~ ""))
write_tsv_(clus_tab, paste0("stats_", PFX, "_clusters.tsv"))

## 4. Fits: chosen, reversal, other-space sensitivity ---------------------------------------
sds   <- run_slingshot(rd_main, cl, root, fit_space)
stopifnot(!is.null(sds))
pt    <- as.numeric(avg_pt(sds))
if (!is.null(sel)) stopifnot(abs(mean(vapply(names(which(tapply(pos, donor, sum) >= MIN_POS_DONOR)), function(d)
  auc_fn(pt[donor == d], pos[donor == d]), numeric(1))) - sel$objective_value) < 1e-9)   # = step 1's fit
pt_rev <- as.numeric(avg_pt(run_slingshot(rd_main, cl, endc, paste(fit_space, "reversal"))) %||% rep(NA_real_, length(pt)))
cl_other   <- kmeans_clusters(rd_other, k_use)
root_other <- pick_cluster(cluster_table(cl_other, pos, donor), "fraction", args$min_clus_cells, "min")
sds_other  <- run_slingshot(rd_other, cl_other, root_other, other_space)
pt_other   <- as.numeric(avg_pt(sds_other) %||% rep(NA_real_, length(pt)))

lineages <- slingshot::slingLineages(sds)
pt_lin   <- slingshot::slingPseudotime(sds, na = TRUE)
w_lin    <- slingshot::slingCurveWeights(sds)
lin_tab <- tibble(lineage = names(lineages), clusters = map_chr(lineages, paste, collapse = " -> "),
                  terminal = map_chr(lineages, ~ .x[length(.x)]),
                  n_cells_weight_gt_0.5 = colSums(w_lin > 0.5),
                  auc_phf1_later = map_dbl(colnames(pt_lin), ~ auc_fn(pt_lin[, .x], pos))) %>%
  left_join(clus_tab %>% select(terminal = cluster, terminal_frac_phf1 = frac_phf1), by = "terminal")
write_tsv_(lin_tab, paste0("stats_", PFX, "_lineages.tsv"))

df <- tibble(cell_id = cd$cell_id, sample_id = donor, Braak = factor(as.character(cd$Braak), levels = braak_levels),
             PHF1 = pos, phf1_intensity_p95_z = cd$phf1_intensity_p95_z, nCount_RNA = cd$nCount_RNA,
             UMAP_1 = umap_mat[, 1], UMAP_2 = umap_mat[, 2], cluster = cl,
             pt = pt, pt_other_space = pt_other, pt_reversal = pt_rev) %>%
  bind_cols(as_tibble(pt_lin) %>% setNames(paste0("pt_", colnames(pt_lin)))) %>%
  bind_cols(as_tibble(w_lin)  %>% setNames(paste0("w_",  colnames(w_lin))))
df$bin <- cut(df$pt, breaks = args$n_bins, labels = FALSE)
write_tsv_(df %>% mutate(fit_space = fit_space, k = k_use, root = root),
           paste0("source_data_", PFX, "_percell.tsv"))

## 5. Description: PHF1+ along pseudotime, bins, what else tracks the axis -------------------
d_auc <- donor_auc(df$pt, df$PHF1, df$sample_id) %>% left_join(distinct(df, donor = sample_id, Braak), by = "donor")
auc_rows <- bind_rows(summ_auc(d_auc, 1, sprintf("%s fit, donors >= 1 PHF1+", fit_space)),
                      summ_auc(d_auc, MIN_POS_DONOR, sprintf("%s fit, donors >= %d PHF1+", fit_space, MIN_POS_DONOR)),
                      summ_auc(donor_auc(df$pt_other_space, df$PHF1, df$sample_id), MIN_POS_DONOR,
                               sprintf("%s fit (sensitivity), donors >= %d PHF1+", other_space, MIN_POS_DONOR)))
write_tsv_(d_auc, paste0("stats_", PFX, "_donor_auc.tsv"))
bin_tab <- df %>% group_by(bin) %>%
  summarise(n = n(), n_phf1 = sum(PHF1), pct_phf1 = 100 * mean(PHF1), pt_lo = min(pt), pt_hi = max(pt), .groups = "drop")
head_donors <- d_auc$donor[d_auc$n_phf1 >= MIN_POS_DONOR]
bin_donor <- df %>% filter(sample_id %in% head_donors) %>% group_by(sample_id, Braak, bin) %>%
  summarise(n = n(), n_phf1 = sum(PHF1), pct_phf1 = 100 * mean(PHF1), .groups = "drop")
bin_braak <- df %>% count(bin, Braak) %>% group_by(bin) %>% mutate(frac = n / sum(n)) %>% ungroup()
d_pt <- df %>% group_by(sample_id, Braak) %>% summarise(median_pt = median(pt), .groups = "drop")
kw   <- kruskal.test(median_pt ~ Braak, data = d_pt)
eps2 <- unname((kw$statistic - nlevels(droplevels(d_pt$Braak)) + 1) / (nrow(d_pt) - nlevels(droplevels(d_pt$Braak))))

## 6. Plots ---------------------------------------------------------------------------------
## Trajectory on the UMAP: slingshot's own curves if it was fitted in UMAP; otherwise each
## lineage as the median UMAP position of its cells (weight > 0.5) in 40 quantile bins of
## that lineage's pseudotime -- a display path, not the fit.
curve_df <- if (fit_space == "umap") {
  bind_rows(imap(slingshot::slingCurves(sds), function(cv, i) {
    s <- cv$s[cv$ord, , drop = FALSE]; tibble(lineage = names(lineages)[i], UMAP_1 = s[, 1], UMAP_2 = s[, 2]) }))
} else {
  bind_rows(lapply(seq_along(lineages), function(i) {
    keep <- which(w_lin[, i] > 0.5 & !is.na(pt_lin[, i]))
    tibble(lineage = names(lineages)[i], b = dplyr::ntile(pt_lin[keep, i], 40),
           UMAP_1 = umap_mat[keep, 1], UMAP_2 = umap_mat[keep, 2]) %>%
      group_by(lineage, b) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") %>%
      arrange(b) %>% select(-b) }))
}
write_tsv_(curve_df, paste0("source_data_", PFX, "_curves.tsv"))
pos_df <- df[df$PHF1, ]

cent <- df %>% group_by(cluster) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") %>%
  left_join(clus_tab, by = "cluster") %>%
  mutate(label = sprintf("%s%s\n%.2f%% PHF1+", cluster,
                         ifelse(cluster == root, " root", ifelse(cluster == endc, " end", "")), 100 * frac_phf1))
p_cl <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(colour = cluster), size = 0.2, stroke = 0) +
  geom_point(data = pos_df, colour = "black", size = 0.5, stroke = 0) +
  geom_path(data = curve_df, aes(group = lineage), colour = "black", linewidth = 0.4) +
  geom_label(data = cent, aes(label = label), size = 2, label.padding = unit(1.2, "pt"), lineheight = 0.9) +
  scale_colour_manual(values = setNames(rep(CLUSTER_COLOURS, length.out = length(unique(cl))), sort(unique(cl))),
                      guide = "none") +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_xy()
save_pdf(p_cl, "umap_clusters", 8, 7)
write_tsv_(df %>% select(cell_id, sample_id, PHF1, UMAP_1, UMAP_2, cluster), paste0("source_data_", PFX, "_umap_clusters.tsv"))

p_pt <- ggplot() +
  geom_point(data = df %>% arrange(pt), aes(UMAP_1, UMAP_2, colour = pt), size = 0.2, stroke = 0) +
  geom_path(data = curve_df, aes(UMAP_1, UMAP_2, group = lineage), colour = "black", linewidth = 0.45) +
  geom_point(data = pos_df, aes(UMAP_1, UMAP_2, shape = "Tangle-bearing"), colour = PHF1_COLOUR,
             fill = NA, size = 1.1, stroke = 0.3) +
  scale_colour_viridis_c(option = "cividis", name = "Pseudotime") +
  scale_shape_manual(values = c("Tangle-bearing" = 21), name = NULL) +
  guides(colour = guide_colourbar(order = 1), shape = guide_legend(order = 2, override.aes = list(size = 1.6, stroke = 0.4))) +
  labs(x = "UMAP 1", y = "UMAP 2") + theme_xy() + theme(legend.key.height = grid::unit(14, "pt"))
save_pdf(p_pt, "umap_pseudotime", 8.5, 6.5)
write_tsv_(df %>% select(cell_id, sample_id, Braak, PHF1, UMAP_1, UMAP_2, cluster, pt),
           paste0("source_data_", PFX, "_umap_pseudotime.tsv"))

p_bin <- ggplot() +
  geom_line(data = bin_donor, aes(bin, pct_phf1, group = sample_id, colour = Braak), linewidth = 0.3, alpha = 0.8) +
  geom_line(data = bin_tab, aes(bin, pct_phf1), linewidth = 0.6) +
  geom_point(data = bin_tab, aes(bin, pct_phf1), size = 1) +
  scale_colour_manual(values = braak_line_palette, name = sprintf("Donor (>= %d PHF1+),\nby Braak", MIN_POS_DONOR)) +
  scale_x_continuous(breaks = seq_len(args$n_bins)) +
  labs(x = sprintf("Pseudotime bin (1-%d, equal width)", args$n_bins), y = "% tangle-bearing") +
  theme_xy() + theme(legend.key.height = grid::unit(7, "pt"))
save_pdf(p_bin, "phf1_by_bin", 8.5, 5.5)
write_tsv_(bind_rows(bin_tab %>% mutate(level = "pooled"), bin_donor %>% mutate(level = "donor", Braak = as.character(Braak))),
           paste0("source_data_", PFX, "_phf1_by_bin.tsv"))

p_braak <- ggplot(bin_braak, aes(bin, frac, fill = Braak)) + geom_col(width = 0.9) +
  scale_fill_manual(values = braak_palette, name = "Braak") + scale_x_continuous(breaks = seq_len(args$n_bins)) +
  labs(x = sprintf("Pseudotime bin (1-%d)", args$n_bins), y = "Fraction of cells") + theme_xy()
save_pdf(p_braak, "braak_by_bin", 7.5, 5)
write_tsv_(bin_braak, paste0("source_data_", PFX, "_braak_by_bin.tsv"))

hd <- auc_rows[auc_rows$donor_set == sprintf("%s fit, donors >= %d PHF1+", fit_space, MIN_POS_DONOR), ]
p_auc <- ggplot(d_auc %>% filter(n_phf1 >= 1), aes(x = 1, y = auc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  geom_jitter(aes(size = n_phf1, fill = Braak), shape = 21, width = 0.12, height = 0, stroke = 0.25) +
  geom_errorbar(data = hd, aes(x = 1.45, ymin = ci_lo, ymax = ci_hi), inherit.aes = FALSE, width = 0.08, linewidth = 0.35) +
  geom_point(data = hd, aes(x = 1.45, y = mean_auc), inherit.aes = FALSE, size = 1.4) +
  scale_fill_manual(values = braak_palette, name = "Braak") + scale_size_area(max_size = 3.5, name = "PHF1+ cells") +
  scale_x_continuous(breaks = c(1, 1.45), labels = c("donors", sprintf("mean,\n>= %d", MIN_POS_DONOR)), limits = c(0.75, 1.65)) +
  labs(x = NULL, y = "AUC (PHF1+ later in pseudotime)") + theme_xy() + theme(legend.key.height = grid::unit(7, "pt"))
save_pdf(p_auc, "donor_auc", 7, 6)
write_tsv_(d_auc, paste0("source_data_", PFX, "_donor_auc.tsv"))

html_file <- NA_character_
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
  df$hover <- with(df, sprintf("%s<br>donor %s (Braak %s)<br>PHF1: %s<br>pseudotime %.2f (bin %d)<br>cluster %s<br>p95 z %.2f",
                               cell_id, sample_id, Braak, ifelse(PHF1, "tangle-bearing", "-"), pt, bin, cluster, phf1_intensity_p95_z))
  cividis_scale <- lapply(seq(0, 1, length.out = 11), function(v) list(v, viridisLite::cividis(11)[round(v * 10) + 1]))
  pi <- plotly::plot_ly(type = "scattergl", mode = "markers") %>%
    plotly::add_trace(data = df %>% arrange(pt), x = ~UMAP_1, y = ~UMAP_2, text = ~hover, hoverinfo = "text",
                      name = "cells", showlegend = FALSE,
                      marker = list(color = ~pt, colorscale = cividis_scale, size = 3, colorbar = list(title = "Pseudotime"))) %>%
    plotly::add_trace(data = df[df$PHF1, ], x = ~UMAP_1, y = ~UMAP_2, text = ~hover, hoverinfo = "text",
                      name = "Tangle-bearing", marker = list(symbol = "circle-open", color = PHF1_COLOUR, size = 9,
                                                             line = list(width = 1.5, color = PHF1_COLOUR)))
  for (ln in unique(curve_df$lineage))
    pi <- plotly::add_trace(pi, data = curve_df[curve_df$lineage == ln, ], x = ~UMAP_1, y = ~UMAP_2, type = "scatter",
                            mode = "lines", name = ln, hoverinfo = "name", line = list(color = "black", width = 2))
  pi <- plotly::layout(pi, xaxis = list(title = "UMAP 1", zeroline = FALSE, showgrid = FALSE),
                       yaxis = list(title = "UMAP 2", zeroline = FALSE, showgrid = FALSE, scaleanchor = "x"),
                       legend = list(orientation = "h", x = 0, y = -0.12), plot_bgcolor = "white")
  html_file <- file.path(OUTPUT_DIR, paste0("plot_", PFX, "_umap_pseudotime_interactive.html"))
  htmlwidgets::saveWidget(pi, normalizePath(html_file, mustWork = FALSE),
                          selfcontained = requireNamespace("rmarkdown", quietly = TRUE) && isTRUE(rmarkdown::pandoc_available()),
                          libdir = "lib", title = paste(CT, "slingshot pseudotime"))
}

## 7. Log -----------------------------------------------------------------------------------
old_w <- options(width = 160)
sink(file.path(OUTPUT_DIR, paste0("stats_", PFX, ".txt")))
cat(CT, " tangle-associated trajectory (R/cbln2_slingshot_trajectory.R; step 2 of 3)\n", strrep("=", 78), "\n", sep = "")
cat("Run:", format(Sys.time()), " slingshot", as.character(packageVersion("slingshot")), "\n")
cat("Settings:", if (!is.null(sel)) sprintf("from %s (created %s)", params_file, sel$created) else
      "default rule (no step-1 export)", sprintf("-> space %s, k %d, root %s\n", fit_space, k_use, root))
cat("HTML:", html_file, "\n")
cat(sprintf("Cells %d (%d PHF1+, %d donors)\n\n", nrow(df), sum(df$PHF1), n_distinct(df$sample_id)))
cat(sprintf("nPC: knee co2 = PC %d, 90%%-cumulative co1 = PC %d -> min = %d; used 1:%d\n\n",
            pcs$co2, pcs$co1, pcs$min_pc, max(use_dims)))
cat("--- Clusters (k-means k =", k_use, "on", fit_space, ", seed 100) ---\n")
print(as.data.frame(clus_tab), digits = 3, row.names = FALSE)
cat("\n--- Lineages ---\n"); print(as.data.frame(lin_tab), digits = 3, row.names = FALSE)
cat(sprintf("\nSpearman rho: pseudotime vs %s-space fit %.3f; vs reversal from cluster %s %.3f\n",
            other_space, spearman(df$pt, df$pt_other_space), endc, spearman(df$pt, df$pt_reversal)))
cat("\n--- PHF1+ along pseudotime (descriptive; selection-adjusted p: R/cbln2_trajectory_selection_null.R) ---\n")
cat(sprintf("Pooled AUC (cell level): %.3f\n", auc_fn(df$pt, df$PHF1)))
print(as.data.frame(auc_rows), digits = 3, row.names = FALSE)
print(as.data.frame(d_auc), digits = 3, row.names = FALSE)
cat("\n--- Bins (equal width) ---\n"); print(as.data.frame(bin_tab), digits = 3, row.names = FALSE)
cat("\n--- What else tracks the axis ---\n")
cat(sprintf("rho(pseudotime, log nCount_RNA) = %.3f; rho(pseudotime, p95 z) in PHF1- cells = %.3f\n",
            spearman(df$pt, log(df$nCount_RNA)), spearman(df$pt[!df$PHF1], df$phf1_intensity_p95_z[!df$PHF1])))
cat(sprintf("Donor median pseudotime by Braak: Kruskal-Wallis H = %.2f, p = %.3g, epsilon^2 = %.2f (n = 3/stage)\n",
            kw$statistic, kw$p.value, eps2))
print(as.data.frame(d_pt), digits = 3, row.names = FALSE)
cat("\nNext: Rscript R/deg_dream_linear_pseudotime_phf1.r --celltype", CT, "(dgenv) reads", ct_percell_file(CT), "\n")
cat("\n--- sessionInfo ---\n"); print(sessionInfo())
sink(); options(old_w)
message(sprintf("Done: %s k=%d root %s; donor-mean AUC %.3f -> %s", fit_space, k_use, root, hd$mean_auc, OUTPUT_DIR))
