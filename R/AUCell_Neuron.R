# ---------------------------------------------------------------------------
# AUCell_Neuron.R
#
# Upstream pipeline - step 8 of 9 | Figure 1C (neuronal subtype marker heatmap)
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
#
# WHAT THIS DOES
#   Round 2 of cell-type annotation, restricted to cells called neuronal in round 1 and
#   scored against neuronal subtype marker sets from the same reference.
#
#   The two rounds are NOT identically parameterised - see the marker-selection and
#   threshold columns of the label-transfer supplementary tables.
#
# INPUTS   AUCell_1/seu.rds, neuronal snRNA-seq reference
# OUTPUTS  AUCell_Neuron/ - gene sets, thresholds, per-cell calls, seu.rds
# ---------------------------------------------------------------------------
# AUCell label transfer from snRNA-seq (sn_seu) to CosMx (seu) — SCT-aware + proportion-adjusted thresholds
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(AUCell)
  library(dplyr)
  library(readr)
})

# ---------------------------
# I/O and parameters
# ---------------------------
setwd("<PROJECT_ROOT>/phf1_v2")
out_dir <- "AUCell_Neuron/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sn_seu <- readRDS("subclustering_round2_corrected/seu_sn_neuron_ref.RDS")  # snRNA ref (SCT)
seu    <- readRDS("AUCell_1/seu.rds")                                                   # CosMx (SCT)

seu <- seu[ ,
            grepl(
              "Glutamatergic|GABAergic|Neuron|Exc|Inh",
              seu$AUCell_label_SCT_propAdj 
            )
]

x <- as.character(sn_seu$subcluster)

sn_seu$broad_subcluster <- case_when(
  # Excitatory — IT L2/3 cloud (CUX2/CALB1/PCP4 + L2 RELN-BMPR1B)
  str_detect(x, "^(Exc-L2-3-CUX2-CALB1|Exc-L2-CUX2-PDGFD|Exc-L3-PCP4-CALB1|Exc-L2-RELN-BMPR1B)$") ~ "Exc-IT-L2-3",
  # Excitatory — IT L4/5 (RORB programs)
  str_detect(x, "^(Exc-L5-RORB-TLL1|Exc-L5-RORB-TPBG)$") ~ "Exc-IT-L3-5",
  # Excitatory — ET L5 (BCL11B/CTIP2)
  str_detect(x, "^Exc-L5-BCL11B-ADRA1A$") ~ "Exc-ET-L5",
  # Excitatory — CT L6 (TLE4; includes L5/6 TLE4-NXPH2)
  str_detect(x, "^(Exc-L6-TLE4-SULF1|Exc-L5-6-TLE4-NXPH2)$") ~ "Exc-CT-L6",
  # Excitatory — IT L6 (THEMIS/CDH13)
  str_detect(x, "^Exc-L6-THEMIS-CDH13$") ~ "Exc-IT-L6",
  
  # Inhibitory — PVALB (merge MYO5B + UNC5B)
  str_detect(x, "^(Inh-PVALB-MYO5B|Inh-PVALB-UNC5B)$") ~ "Inh-PVALB",
  # Inhibitory — SST (include SST-NPY)
  str_detect(x, "^(Inh-SST-NPY)$") ~ "Inh-SST",
  # Inhibitory — VIP
  str_detect(x, "^(Inh-VIP-RELN)$") ~ "Inh-VIP",
  # Inhibitory — CGE RELN family (merge RELN,SST into LAMP5/RELN)
  str_detect(x, "^(Inh-LAMP5-RELN|Inh-RELN,SST)$") ~ "Inh-LAMP5",
  
  TRUE ~ "Other"
)

# Choose the label column in the snRNA reference
label_col <- dplyr::coalesce(
  intersect(c("broad_subcluster"), colnames(sn_seu@meta.data))[1],
  "subcluster"
)
message("Using snRNA label column: ", label_col)

# Gene set and DE parameters
min_genes_per_set <- 10
top_n_markers     <- 60      # top-N up markers per label
fdr_cutoff        <- 0.05
min_pct           <- 0.20
logfc_thresh      <- 1
ncores_rank       <- max(1, parallel::detectCores() - 1)

# If you already have a markers CSV from the same reference + SCT, set path here to skip FindAllMarkers
markers_csv <- NA_character_   # e.g., "snrna_findallmarkers_SCT.csv"

# ---------------------------
# Shared features (both already subset to CosMx panel, but we check)
# ---------------------------
DefaultAssay(sn_seu) <- "SCT"
DefaultAssay(seu)    <- "SCT"

shared <- intersect(rownames(sn_seu[["SCT"]]), rownames(seu[["SCT"]]))
stopifnot(length(shared) >= 100)
message("Shared genes (SCT): ", length(shared))

# ---------------------------
# Build marker gene sets per snRNA label (SCT assay)
# ---------------------------
if (!is.na(markers_csv) && file.exists(markers_csv)) {
  message("Reading markers from CSV: ", markers_csv)
  markers <- read_csv(markers_csv, show_col_types = FALSE)
  # Expect columns: cluster, gene, avg_log2FC, p_val_adj
  markers <- markers %>%
    mutate(gene = as.character(gene)) %>%
    filter(avg_log2FC > 0, p_val_adj <= fdr_cutoff, gene %in% shared) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n_markers, with_ties = FALSE) %>%
    ungroup()
} else {
  message("Computing markers from SCT assay in the snRNA reference …")
  Idents(sn_seu) <- sn_seu[[label_col]][,1]
  sn_seu_filt <- subset(sn_seu, features = shared)
  markers <- FindAllMarkers(
    sn_seu_filt,
    only.pos        = TRUE,
    min.pct         = min_pct,
    logfc.threshold = logfc_thresh,
    test.use        = "wilcox"
  )
  markers <- markers %>%
    filter(p_val_adj <= fdr_cutoff, gene %in% shared) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n_markers, with_ties = FALSE) %>%
    ungroup()
}
stopifnot(nrow(markers) > 0)

# Build gene sets (list of character vectors)
label_names <- sort(unique(markers$cluster))
gene_sets <- lapply(label_names, function(cl) unique(markers$gene[markers$cluster == cl]))
names(gene_sets) <- label_names
gene_sets <- gene_sets[sapply(gene_sets, length) >= min_genes_per_set]
stopifnot(length(gene_sets) > 0)
saveRDS(gene_sets, file.path(out_dir, "gene_sets_from_SCT_snRNA.rds"))
message("Gene sets created for labels: ", paste(names(gene_sets), collapse = ", "))

# ---------------------------
# AUCell rankings on CosMx — use non-negative matrix
# ---------------------------
get_sct_counts <- function(obj) {
  if (!is.null(obj[["SCT"]]@counts) && nrow(obj[["SCT"]]@counts) > 0) {
    as.matrix(obj[["SCT"]]@counts)
  } else if (!is.null(obj[["RNA"]]) && nrow(obj[["RNA"]]@counts) > 0) {
    as.matrix(obj[["RNA"]]@counts)
  } else {
    m <- as.matrix(obj[["SCT"]]@data); m[m < 0] <- 0; m
  }
}

expr_cosmx <- get_sct_counts(seu)
expr_cosmx <- expr_cosmx[intersect(rownames(expr_cosmx), shared), , drop = FALSE]

message("Building AUCell rankings (SCT counts) …")
rankings <- AUCell_buildRankings(expr_cosmx, plotStats = FALSE)

message("Calculating AUCell AUCs …")
auc <- AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(0.06 * nrow(rankings)))
auc_mat <- t(getAUC(auc))   # cells x labels
stopifnot(nrow(auc_mat) == ncol(seu))
saveRDS(auc_mat, file.path(out_dir, "cosmx_AUCell_matrix_SCT.rds"))

# ---------------------------
# Thresholds per label (from CosMx AUC) with robust fallback
# ---------------------------
get_label_threshold <- function(auc_vec, fallback_q = 0.90) {
  th <- tryCatch({
    thr <- AUCell_exploreThresholds(auc_vec, plotHist = FALSE, assign = FALSE)
    if (!is.null(thr[[1]]$selected)) thr[[1]]$selected else tail(thr[[1]]$thresholds, 1)
  }, error = function(e) NA_real_)
  if (is.na(th)) stats::quantile(auc_vec, probs = fallback_q, na.rm = TRUE) else th
}

thr_tbl <- data.frame(
  label     = colnames(auc_mat),
  threshold = sapply(colnames(auc_mat), function(lbl) get_label_threshold(auc_mat[, lbl], fallback_q = 0.90)),
  row.names = colnames(auc_mat)
)

# Clamp thresholds to avoid extremes
thr_floor   <- 0.027
thr_ceiling <- max(0.99, quantile(as.numeric(auc_mat), 0.99, na.rm = TRUE))
thr_tbl$threshold <- pmin(pmax(thr_tbl$threshold, thr_floor), thr_ceiling)

# ---------------------------
# Proportion adjustment (based on snRNA reference label proportions)
# ---------------------------
# Compute priors from sn_seu
pri_tab <- prop.table(table(sn_seu[[label_col]][,1]))
p_raw <- as.numeric(pri_tab[colnames(auc_mat)])  # align to AUCell label order
names(p_raw) <- colnames(auc_mat)
# Fill missing with mean prior to avoid NA
if (any(is.na(p_raw))) {
  warning("Missing priors for: ", paste(names(p_raw)[is.na(p_raw)], collapse = ", "), " — using mean prior")
  p_raw[is.na(p_raw)] <- mean(p_raw, na.rm = TRUE)
}
# Normalize and set softness (gamma controls influence on thresholds)
p_norm <- p_raw / max(p_raw, na.rm = TRUE)
gamma  <- 0.45               # 0 = no influence, 0.5 = sqrt-like, 1 = linear
p_soft <- p_norm^gamma

# Adjust thresholds: common labels (high p_soft) get slightly lower thresholds; rare get higher
alpha  <- 0.11              # magnitude of adjustment; try 0.10–0.20
mu     <- mean(p_soft, na.rm = TRUE)
thr_adj <- thr_tbl$threshold - alpha * (p_soft - mu)
# Enforce floors/ceilings
thr_adj <- pmin(pmax(thr_adj, thr_floor), thr_ceiling)
thr_tbl$threshold_adj <- thr_adj

write.csv(thr_tbl, file.path(out_dir, "AUCell_label_thresholds_SCT_with_priors.csv"), row.names = FALSE)

# Optional absolute AUC floor to avoid selecting labels on tiny thresholds
score_floor <- 0.027

# ---------------------------
# Assign labels by ratio = AUC / threshold_adj
# If multiple pass (ratio >= 1 and AUC >= score_floor), choose the highest ratio; else "Unassigned".
# ---------------------------
thr_map <- setNames(thr_tbl$threshold_adj, thr_tbl$label)

R <- sweep(auc_mat, 2, thr_map[colnames(auc_mat)], "/")
pass <- (R >= 1) & (auc_mat >= matrix(score_floor, nrow(auc_mat), ncol(auc_mat)))

R_sel <- R
R_sel[!pass] <- -Inf

best_idx   <- max.col(R_sel, ties.method = "first")
best_lab   <- colnames(R_sel)[best_idx]
best_ratio <- R_sel[cbind(seq_len(nrow(R_sel)), best_idx)]
none_pass  <- !is.finite(best_ratio)

final_label <- ifelse(none_pass, "Unassigned Neuron", best_lab)

# Attach to Seurat and save
seu$AUCell_neuron_label_SCT_propAdj      <- final_label
seu$AUCell_neuron_best_ratio_SCT_propAdj <- ifelse(none_pass, NA_real_, best_ratio)

out_calls <- data.frame(
  cell        = colnames(seu),
  final_label = final_label,
  best_ratio  = seu$AUCell_best_ratio_SCT_propAdj,
  stringsAsFactors = FALSE
)
write.csv(out_calls, file.path(out_dir, "cosmx_AUCell_labels_SCT_propAdj.csv"), row.names = FALSE)
saveRDS(seu, file.path(out_dir, "seu.rds"))

message("Assigned (non-Unassigned) cells: ", sum(final_label != "Unassigned Neuron"),
        " / ", length(final_label))
print(round(prop.table(table(final_label)), 3))
print(round(prop.table(table(sn_seu$broad_subcluster)), 3))


# Heatmap

marker_genes <- c( # from Franjic et al., 2022)
  # Excitatory markers
  "CALB1",    # calbindin+
  "IL1RAPL2",
  "LAMA3",
  "PDGFD",
  "RELN",     # Layer II reelin+
  "BCL11B",
  "BMPR1B",   # Layer II
  "RORB",     # IT
  "PCP4",
  "CDH13",
  "RGS12",
  "ADRA1A",
  "TLE4",     # Corticothalamic
  "CCN2",
  
  
  # Inhibitory markers
  "PVALB",    # Fast-spiking
  "SST",      # Martinotti-like
  "VIP",      # Disinhibitory
  "NPY",      
  "UNC5B",   # Chandelier
  "MYO5B",   # Basket
  "CALB2",
  
  # General markers
  "SLC17A7",  # Excitatory
  "CAMK2A",
  "GAD1" ,    # Inhibitory
  "SLC32A1"  # Vesicular GABA transporter
)

#  heatmap

avg_expr <- AverageExpression(
  seu,
  features = marker_genes,
  group.by = "AUCell_neuron_label_SCT_propAdj",
  assays = "SCT",
  slot = "data"
)$SCT

# Z-score per gene 
mat_z <- t(scale(t(avg_expr)))

# Hierarchical clustering
gene_hc <- hclust(dist(mat_z, method = "euclidean"))

gene_order <- rownames(mat_z)[gene_hc$order]

avg_df <- as.data.frame(avg_expr) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "cluster", values_to = "mean_expression") |>
  group_by(gene) |>
  mutate(z_expr = as.numeric(scale(mean_expression))) |>
  ungroup() |>
  mutate(
    gene = factor(gene, levels = gene_order)
  )

ggplot(avg_df, aes(cluster, gene, fill = z_expr)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Z-score"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Mean gene expression per Cluster",
    x = "Predicted Celltype",
    y = "Gene"
  )

ggsave(filename = paste0(out_dir,"/heatmap_AUCell.png"),
              width = 14, height = 10, dpi = 300)
