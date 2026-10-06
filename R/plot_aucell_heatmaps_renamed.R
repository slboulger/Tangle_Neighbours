# ---------------------------------------------------------------------------
# plot_aucell_heatmaps_renamed.R
#
# Figure panels: 1B, 1C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# Replot AUCell label-transfer heatmaps using renamed celltypes from seu_PHF1.rds.
# Same gene sets / averaging / z-score / gene hclust / ggplot style as
# AUCell_1.R and AUCell_Neuron.R. Only deliberate change: larger axis text so
# the heatmaps can be shrunk on the figure page.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")  # fig_theme + palettes

seu <- readRDS("seu_PHF1.rds")

message("Celltype tally on seu_PHF1.rds:")
print(table(seu$celltype))

# Output dirs for source-data + stats triple (alongside the existing PDF locations)
dir.create("plots/aucell_heatmaps", showWarnings = FALSE, recursive = TRUE)

# Shared helper: build the avg-expression -> z-score -> hclust-ordered long df,
# then make and save the same ggplot tile heatmap used in AUCell_1.R / _Neuron.R.
# Also writes the source-data TSV and stats log (output triple).
make_heatmap <- function(obj, marker_genes, group_col, out_path, slug,
                         cluster_levels = NULL, left_margin_pt = 38) {
  avg_expr <- AverageExpression(
    obj,
    features = marker_genes,
    group.by = group_col,
    assays   = "SCT",
    slot     = "data"
  )$SCT

  mat_z      <- t(scale(t(avg_expr)))
  gene_hc    <- hclust(dist(mat_z, method = "euclidean"))
  gene_order <- rownames(mat_z)[gene_hc$order]

  avg_df <- as.data.frame(avg_expr) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "cluster", values_to = "mean_expression") |>
    group_by(gene) |>
    mutate(z_expr = as.numeric(scale(mean_expression))) |>
    ungroup() |>
    mutate(gene = factor(gene, levels = gene_order))

  # Order the x-axis (celltype/cluster) to match the canonical palettes.R
  # ordering used across the rest of the figures. Any label not present in the
  # supplied order is appended (never dropped) so no column is silently lost.
  if (!is.null(cluster_levels)) {
    present <- unique(as.character(avg_df$cluster))
    extra   <- setdiff(present, cluster_levels)
    if (length(extra) > 0) {
      message("  note: clusters absent from supplied order, appended at end: ",
              paste(extra, collapse = ", "))
    }
    x_order <- c(intersect(cluster_levels, present), extra)
    avg_df  <- mutate(avg_df, cluster = factor(cluster, levels = x_order))
  }

  p <- ggplot(avg_df, aes(cluster, gene, fill = z_expr)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "Z-score"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic(base_size = 8) +
    fig_theme +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    # Ensure axis text is size 7 (keeps the 45deg x-angle from fig_theme).
    # Override the left margin per-plot: fig_theme reserves l = 38pt for long
    # rotated x labels (needed for neuron subtypes, not for short broad labels).
    theme(axis.text = element_text(size = 7, colour = "black"),
          plot.margin = margin(t = 4, r = 6, b = 10, l = left_margin_pt,
                               unit = "pt"))

  ggsave(filename = out_path, plot = p,
         width = 8, height = 9.38, units = "cm",
         device = "pdf")
  message("Saved: ", out_path)

  # ---- source data: long-format (gene, celltype, mean_expression, z_expr) ----
  src_path <- file.path("plots/aucell_heatmaps", paste0("source_data_", slug, ".tsv"))
  src_df <- avg_df |>
    mutate(gene = as.character(gene)) |>
    arrange(gene, cluster)
  write.table(src_df, src_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: ", src_path)

  # ---- stats log: descriptive (no formal test) ----
  log_path <- file.path("plots/aucell_heatmaps", paste0("stats_", slug, ".txt"))
  sink(log_path)
  cat("AUCell label-transfer heatmap — descriptive (no formal statistical test).\n")
  cat("Slug:        ", slug, "\n")
  cat("Grouped by:  ", group_col, "\n")
  cat("Date:        ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Cells per group:\n")
  print(table(obj@meta.data[[group_col]]))
  cat("\nGene clustering: hclust on euclidean distance of per-gene z-scored means.\n")
  cat("Gene order (top -> bottom on plot, i.e. dendrogram leaf order):\n")
  cat(paste(gene_order, collapse = ", "), "\n\n")
  cat("Marker gene set (input order):\n")
  cat(paste(marker_genes, collapse = ", "), "\n\n")
  cat("sessionInfo():\n")
  print(sessionInfo())
  sink()
  message("Saved: ", log_path)
}

# ---- 1) Broad heatmap (AUCell_1.R style) ----
broad_marker_genes <- c(
  # Neurons
  "RBFOX3","SNAP25","SYT1","MAP2",
  "SLC17A7","CAMK2A",
  "GAD1","GAD2","SLC6A1",

  # Astro
  "GFAP","ALDH1L1","AQP4","SLC1A2","SLC1A3","GLUL","SPARCL1",

  # Micro
  "P2RY12","CX3CR1","TYROBP","C1QA",

  # Oligo lineage
  "OLIG2","PDGFRA","CSPG4",
  "MBP","PLP1","MOG",

  # Vascular
  "PECAM1","VWF","KDR","CLDN5",
  "COL1A1","DCN"
)

# Prefer the existing broad AUCell column; otherwise derive from `celltype`.
if ("AUCell_label_SCT_propAdj" %in% colnames(seu@meta.data)) {
  broad_col <- "AUCell_label_SCT_propAdj"
  message("Broad heatmap: grouping by AUCell_label_SCT_propAdj")
} else if ("celltype" %in% colnames(seu@meta.data)) {
  ct <- as.character(seu$celltype)
  seu$broad_celltype_for_plot <- dplyr::case_when(
    grepl("^Exc", ct) ~ "Glutamatergic",
    grepl("^Inh", ct) ~ "GABAergic",
    TRUE              ~ ct
  )
  broad_col <- "broad_celltype_for_plot"
  message("Broad heatmap: AUCell_label_SCT_propAdj absent; derived broad labels from `celltype`.")
} else {
  stop("Neither `AUCell_label_SCT_propAdj` nor `celltype` is present on seu_PHF1.rds.")
}

make_heatmap(
  obj          = seu,
  marker_genes = broad_marker_genes,
  group_col    = broad_col,
  out_path     = "AUCell_1/heatmap_AUCell_renamed.pdf",
  slug         = "broad_heatmap",
  # Neuron families first, then glia_order, then catch-alls — matches the
  # broad-label ordering implied by palettes.R (celltype_palette / glia_order).
  cluster_levels = c("Glutamatergic", "GABAergic",
                     glia_order, "Unassigned", "Unassigned Neuron"),
  # Short broad labels don't need the extra left gap fig_theme reserves.
  left_margin_pt = 2
)

# ---- 2) Neuron heatmap (AUCell_Neuron.R style) ----
neuron_marker_genes <- c( # from Franjic et al., 2022
  # Excitatory markers
  "CALB1",
  "IL1RAPL2",
  "LAMA3",
  "PDGFD",
  "RELN",
  "BCL11B",
  "BMPR1B",
  "RORB",
  "PCP4",
  "CDH13",
  "RGS12",
  "ADRA1A",
  "TLE4",
  "CCN2",

  # Inhibitory markers
  "PVALB",
  "SST",
  "VIP",
  "NPY",
  "UNC5B",
  "MYO5B",
  "CALB2",

  # General markers
  "SLC17A7",
  "CAMK2A",
  "GAD1",
  "SLC32A1"
)

# Subset to neurons via the renamed `celltype` column (where the new names live)
neuron_cells <- grepl("^(Exc|Inh)-", as.character(seu$celltype)) |
                seu$celltype == "Unassigned Neuron"
message("Neuron cells: ", sum(neuron_cells), " / ", length(neuron_cells))

seu_neuron <- seu[, neuron_cells]
seu_neuron$celltype <- droplevels(factor(seu_neuron$celltype))

make_heatmap(
  obj          = seu_neuron,
  marker_genes = neuron_marker_genes,
  group_col    = "celltype",
  out_path     = "AUCell_Neuron/heatmap_AUCell_renamed.pdf",
  slug         = "neuron_heatmap",
  # Canonical neuron_order (Exc layers then Inh subtypes) + Unassigned Neuron.
  cluster_levels = c(neuron_order, "Unassigned Neuron")
)
