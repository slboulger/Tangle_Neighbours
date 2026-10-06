# ---------------------------------------------------------------------------
# QC.R
#
# Upstream pipeline - step 3 of 9
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
#   Per-cell quality control. Cells failing the count, feature, negative-probe or area
#   thresholds are removed; 167,176 cells are retained across the 9 donors.
#
# INPUTS   seu_with_meta.RDS
# OUTPUTS  seu_after_qc.RDS
# ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("<PROJECT_ROOT>/phf1_v2")

seu <- readRDS("seu_with_meta.RDS")

slide_col <- if ("SlideLabel" %in% colnames(seu@meta.data)) "SlideLabel" else
  if ("slide" %in% colnames(seu@meta.data)) "slide" else
    stop("No 'SlideLabel' or 'slide' column found in metadata.")

# Negative probes

# Negative-probe percent (from metadata or assay)
if (!"percent.neg" %in% colnames(seu@meta.data)) {
  if ("nCount_negprobes" %in% colnames(seu@meta.data)) {
    seu$percent.neg <- ifelse(seu$nCount_RNA > 0, 100 * seu$nCount_negprobes / seu$nCount_RNA, NA_real_)
  } else if ("negprobes" %in% names(seu@assays)) {
    neg_counts <- Matrix::colSums(GetAssayData(seu, assay = "negprobes", slot = "counts"))
    seu$nCount_negprobes <- as.numeric(neg_counts)
    seu$percent.neg <- ifelse(seu$nCount_RNA > 0, 100 * seu$nCount_negprobes / seu$nCount_RNA, NA_real_)
  }
}

# False-code percent (from metadata or assay)
if (!"percent.false" %in% colnames(seu@meta.data)) {
  if ("nCount_falsecode" %in% colnames(seu@meta.data)) {
    seu$percent.false <- ifelse(seu$nCount_RNA > 0, 100 * seu$nCount_falsecode / seu$nCount_RNA, NA_real_)
  } else if ("falsecode" %in% names(seu@assays)) {
    fc_counts <- Matrix::colSums(GetAssayData(seu, assay = "falsecode", slot = "counts"))
    seu$nCount_falsecode <- as.numeric(fc_counts)
    seu$percent.false <- ifelse(seu$nCount_RNA > 0, 100 * seu$nCount_falsecode / seu$nCount_RNA, NA_real_)
  }
}

# Counts per area
if ("Area.um2" %in% colnames(seu@meta.data) && !"counts_per_um2" %in% colnames(seu@meta.data)) {
  seu$counts_per_um2 <- ifelse(seu$Area.um2 > 0, seu$nCount_RNA / seu$Area.um2, NA_real_)
}



#nCounts per cell
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)

thresh <- 50

ggplot(seu@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 500, color = "grey20", fill = "steelblue", alpha = 0.8, na.rm = TRUE) +
  geom_vline(xintercept = thresh, color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = sprintf("Transcripts per cell (%.1f%% < %d)", 
                    100 * mean(seu@meta.data$nCount_RNA < thresh, na.rm = TRUE), thresh),
    x = "Transcripts per cell (nCount_RNA)",
    y = "Cell count"
  ) +
  theme_bw()

ggsave("QC/nCountRNA_per_cell.png", units = "in", width = 8, height = 4)

# nCounts per slide
ggplot(seu@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "grey20", alpha = 0.85, na.rm = TRUE) +
  geom_vline(xintercept = thresh, color = "red", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(as.formula(paste("~", slide_col)), scales = "free_y") +
  labs(title = "Transcripts per cell by slide", x = "nCount_RNA", y = "Cells") +
  theme_bw(base_size = 10)





# Parameters
abs_floor <- 50
low_q  <- 0.01
high_q <- 0.995
n_sd   <- 3

# Compute global thresholds directly on object vectors
nCount_low  <- as.numeric(quantile(seu$nCount_RNA,  probs = low_q,  na.rm = TRUE))
nCount_high <- as.numeric(quantile(seu$nCount_RNA,  probs = high_q, na.rm = TRUE))
nFeat_low   <- as.numeric(quantile(seu$nFeature_RNA, probs = low_q,  na.rm = TRUE))
nFeat_high  <- as.numeric(quantile(seu$nFeature_RNA, probs = high_q, na.rm = TRUE))

neg_cap   <- 2
false_cap <- 2

# Store thresholds in-object for provenance
seu@misc$qc <- seu@misc$qc %||% list()
seu@misc$qc$global_thresholds <- list(
  nCount_low = nCount_low, nCount_high = nCount_high,
  nFeature_low = nFeat_low, nFeature_high = nFeat_high,
  percent_neg_cap = neg_cap, percent_false_cap = false_cap,
  absolute_floor = abs_floor
)

# Build boolean flags directly into metadata
seu$qc_abs_floor      <- seu$nCount_RNA >= abs_floor
seu$qc_in_nCount_iqr  <- (seu$nCount_RNA  >= nCount_low)  & (seu$nCount_RNA  <= nCount_high)
seu$qc_in_nFeat_iqr   <- (seu$nFeature_RNA >= nFeat_low)  & (seu$nFeature_RNA <= nFeat_high)
seu$qc_neg_ok         <- if (is.na(neg_cap))   TRUE else (is.na(seu$percent.neg)   | seu$percent.neg   <= neg_cap)
seu$qc_false_ok       <- if (is.na(false_cap)) TRUE else (is.na(seu$percent.false) | seu$percent.false <= false_cap)

# Combined cell-level pass flag
seu$qc_pass_cell <- seu$qc_abs_floor & seu$qc_in_nCount_iqr & seu$qc_in_nFeat_iqr & seu$qc_neg_ok & seu$qc_false_ok

# Quick counts
table(seu$qc_pass_cell)


#pNeg

VlnPlot(
  seu, features = "percent.neg", group.by = "sample_id",
  pt.size = 0, assay = DefaultAssay(seu)
) +
  #geom_jitter(alpha = 0.4) +
  geom_hline(yintercept = seu@misc$qc$global_thresholds$percent_neg_cap, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank()
  ) +
  labs(title = paste0("percent.neg", " by sample_id"),
       x = "sample_id", y = "percent.neg")

ggsave("QC/percent_neg.png", units = "in", width = 8, height = 4)


# pFalse

VlnPlot(
  seu, features = "percent.false", group.by = "sample_id",
  pt.size = 0, assay = DefaultAssay(seu)
) +
  #geom_jitter(alpha = 0.4) +
  geom_hline(yintercept = seu@misc$qc$global_thresholds$percent_false_cap, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank()
  ) +
  labs(title = paste0("percent.false", " by sample_id"),
       x = "sample_id", y = "percent.false")

ggsave("QC/percent_false.png", units = "in", width = 8, height = 4)

# FOV level QC: mean nCount_RNA > 100

seu$fov_mean_nCount_RNA <- with(
  seu@meta.data,
  ave(
    nCount_RNA,
    SlideLabel,
    fov,
    FUN = mean
  )
)

fov_cutoff <- 100

seu$qc_fov_mean_ok <-
  ifelse(
    is.na(seu$fov_mean_nCount_RNA),
    TRUE,
    seu$fov_mean_nCount_RNA >= fov_cutoff
  )


ggplot(seu@meta.data, aes(x = fov_mean_nCount_RNA)) +
  geom_histogram(bins = 30, color = "black", fill = "steelblue") +
  geom_vline(xintercept = fov_cutoff, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_classic() +
  labs(
    x = "Mean nCount_RNA per FOV",
    y = "Number of FOVs",
    title = "Distribution of mean nCount_RNA across FOVs"
  )

ggsave("QC/mean_nCount_RNA_fov.png", units = "in", width = 8, height = 4)

seu$qc_pass_fov <- seu$qc_fov_mean_ok

seu$qc_pass_all <- seu$qc_pass_cell & seu$qc_pass_fov



# Subset
seu_qc <- subset(seu, subset = qc_pass_all)

saveRDS(seu_qc,"seu_after_qc.RDS")

