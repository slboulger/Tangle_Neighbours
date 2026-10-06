# ---------------------------------------------------------------------------
# plot_dist_to_phf1_braak_celltype.R
#
# Figure panels: 6A - upstream step
# Writes the per-cell distances, glia summary and leave-one-out contrasts
# (plots/dist_to_phf1/) that plot_dist_to_phf1_by_glia_vs_null.R reads.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# Mean distance to the nearest PHF1+ neuron, as a readout of how cells sit
# relative to ptau pathology. Two panels, each with the output triple
# (plot PDF + source_data TSV + stats log):
#
# (a) dist_to_phf1_by_braak    — per-donor mean distance (over ALL PHF1-negative
#     cells) vs Braak stage. Expectation: as Braak rises, PHF1+ neurons become
#     denser, so the average cell sits closer -> mean distance falls.
# (b) dist_to_phf1_by_celltype — mean distance compared across all celltypes
#     (which celltypes are positioned closest to ptau pathology).
#
# Data semantics:
#   dist_to_phf1_um = Euclidean distance (um) to the nearest PHF1+ NEURON within
#   the same sample. It is NA for PHF1+ cells, so filtering to non-NA distance
#   restricts to the PHF1-negative population automatically.
#
# Stats:
# - Panel (a): Braak is a donor-level variable (each donor in exactly one stage),
#   so the comparison is donor-independent -> Kruskal-Wallis + Dunn (BH) on the
#   9 per-donor means. n=3/group is underpowered; non-significant is not evidence
#   of absence.
# - Panel (b): celltype is donor-paired (each donor contributes cells to every
#   celltype) -> cell-level linear mixed model dist ~ celltype + (1|sample_id) with
#   emmeans leave-one-out contrasts (each celltype vs the average of all others,
#   BH-adjusted). The unit drawn on both panels is the donor.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(emmeans)
  library(rstatix)
  library(ggpubr)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

out_dir <- "plots/dist_to_phf1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Micron axis label rendered via plotmath so the mu glyph draws from the font's
# symbol table under the default pdf device.
y_lab_dist <- expression("Distance to PHF1+ neuron (" * mu * "m)")

# Distance axis is log10-scaled; ticks sit at log positions but are labelled with
# the actual micron values (not log values).
dist_breaks <- c(1, 3, 10, 30, 100, 300, 1000, 3000)

# -------------------------------------------------------------------
# Load + sanity check
# -------------------------------------------------------------------
seu <- readRDS("seu_PHF1.rds")

stopifnot("dist_to_phf1_um column missing - run label_phf1_neighbours.r first" =
            "dist_to_phf1_um" %in% colnames(seu@meta.data))

message("table(sample_id, Braak):"); print(table(seu$sample_id, seu$Braak))
message("dist_to_phf1_um summary (all cells):"); print(summary(seu$dist_to_phf1_um))
message("NA dist by PHF1 status (NA expected only for PHF1+):")
print(table(PHF1 = seu$PHF1, dist_NA = is.na(seu$dist_to_phf1_um)))

# Coerce types we rely on; keep only PHF1-negative cells with a defined distance.
md <- seu@meta.data |>
  mutate(
    sample_id        = as.character(sample_id),
    celltype         = as.character(celltype),
    Braak            = factor(as.character(Braak), levels = braak_levels),
    PHF1             = as.logical(PHF1),
    dist_to_phf1_um  = as.numeric(dist_to_phf1_um)
  ) |>
  select(sample_id, Braak, celltype, PHF1, dist_to_phf1_um)

md_dist <- md |> filter(!is.na(dist_to_phf1_um))

# The distance axis is log10-scaled; flag any non-positive distances (they would be
# dropped by the log transform). Cell centroids differ, so this is expected to be 0.
message("Non-positive distances (dropped on log axis): ",
        sum(md_dist$dist_to_phf1_um <= 0, na.rm = TRUE))

# Cell-level distances underlying the violins on every panel below (all
# PHF1-negative cells). Written once and shared; the per-plot source_data files
# hold the per-donor means (the dots / crossbar).
write.table(
  md_dist |> select(sample_id, Braak, celltype, dist_to_phf1_um),
  file.path(out_dir, "source_data_dist_to_phf1_cells.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =============================================================================
# Panel (a) — per-donor mean distance to PHF1+ neuron vs Braak (KW + Dunn)
# =============================================================================
slug_a <- "dist_to_phf1_by_braak"

# Per-donor mean distance over ALL PHF1-negative cells (donor = inferential unit)
dist_by_donor <- md_dist |>
  group_by(sample_id, Braak) |>
  summarise(mean_dist = mean(dist_to_phf1_um), n_cells = n(), .groups = "drop")

# Kruskal-Wallis + Dunn (BH) on the per-donor means (the inferential unit)
kw_a <- dist_by_donor |> kruskal_test(mean_dist ~ Braak)
dunn_a <- dist_by_donor |> dunn_test(mean_dist ~ Braak, p.adjust.method = "BH")

# Donor-level mean +/- SEM (drawn as a crossbar over the cell-level violin)
donor_summ_a <- dist_by_donor |>
  group_by(Braak) |>
  summarise(m = mean(mean_dist), sem = sd(mean_dist) / sqrt(n()), .groups = "drop") |>
  mutate(xc = as.integer(Braak))

# Best-of-both display: cell-level violin (spread of ALL PHF1-negative cells,
# full range) behind the per-donor mean dots and a donor mean +/- SEM marker.
pa <- ggplot() +
  geom_violin(data = md_dist, aes(x = Braak, y = dist_to_phf1_um, fill = Braak),
              scale = "width", linewidth = 0.2, alpha = 0.5, colour = NA) +
  geom_jitter(data = dist_by_donor, aes(x = Braak, y = mean_dist),
              width = 0.1, height = 0, size = 0.9, colour = "black", alpha = 0.9) +
  geom_errorbar(data = donor_summ_a, aes(x = Braak, ymin = m - sem, ymax = m + sem),
                inherit.aes = FALSE, width = 0.16, linewidth = 0.3) +
  geom_segment(data = donor_summ_a, aes(x = xc - 0.2, xend = xc + 0.2, y = m, yend = m),
               inherit.aes = FALSE, linewidth = 0.4) +
  scale_fill_manual(values = braak_palette, guide = "none") +
  scale_y_log10(breaks = dist_breaks, labels = as.character(dist_breaks),
                expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "Braak", y = y_lab_dist) +
  coord_flip() +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))

# Dunn brackets stacked above the violins if any pairwise comparison is significant
sig_a <- dunn_a |> filter(p.adj < 0.05)
if (nrow(sig_a) > 0) {
  base_a <- max(md_dist$dist_to_phf1_um, na.rm = TRUE)
  sig_a <- sig_a |> mutate(y.position = base_a * (1.03 + 0.07 * (row_number() - 1)))
  pa <- pa + stat_pvalue_manual(
    sig_a, label = "p.adj.signif", coord.flip = TRUE,
    tip.length = 0.005, bracket.size = 0.25, label.size = 2.2
  )
}

ggsave(file.path(out_dir, paste0("plot_", slug_a, ".pdf")), pa,
       width = 8.9, height = 4.5, units = "cm", device = "pdf")

write.table(dist_by_donor, file.path(out_dir, paste0("source_data_", slug_a, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug_a, ".txt")))
cat("Panel (a) — per-donor mean distance to nearest PHF1+ neuron vs Braak\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("dist_to_phf1_um = distance (um) to the nearest PHF1+ neuron in the same\n")
cat("sample; NA for PHF1+ cells, so this is the PHF1-negative population.\n")
cat("Per-donor mean computed over ALL PHF1-negative cells (any celltype).\n\n")
cat("Figure: cell-level violin (all PHF1-negative cells, full y-range) behind the\n")
cat("per-donor mean dots and a donor mean +/- SEM crossbar. The violin shows the\n")
cat("cell-level spread; the dots and crossbar are the donor-level inferential unit\n")
cat("the test below acts on. Cell-level values: source_data_dist_to_phf1_cells.tsv.\n\n")
cat("Kruskal-Wallis is appropriate: each donor is in exactly one Braak group.\n")
cat("n=3 donors per Braak group -- underpowered; non-significant results are not\n")
cat("evidence of absence.\n\n")
cat("Per-donor table (mean distance, um):\n"); print(as.data.frame(dist_by_donor))
cat("\nKruskal-Wallis:\n"); print(kw_a)
cat("\nDunn pairwise (BH-adjusted):\n"); print(dunn_a)
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

dunn_a_out <- as.data.frame(dunn_a) |>
  (\(x) x[, !vapply(x, is.list, logical(1)), drop = FALSE])()
write.table(dunn_a_out, file.path(out_dir, paste0("stats_", slug_a, "_dunn.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Done panel (a).")

# =============================================================================
# Panel (b) — mean distance to PHF1+ neuron across celltypes (lmer + emmeans LOO)
# =============================================================================
slug_b <- "dist_to_phf1_by_celltype"

# Celltype ordering (neurons first, then glia), restricted to celltypes actually
# present among the PHF1-negative cells. "Unassigned" and "Unassigned Neuron" are
# excluded entirely: they are dropped from the figure AND the mixed model, so the
# leave-one-out "average of all other celltypes" is computed only over annotated types.
all_ct_levels <- intersect(
  c(neuron_order, glia_order, "Other"),
  unique(md_dist$celltype)
)

# Cell-level model frame
md_dist_ct <- md_dist |>
  filter(celltype %in% all_ct_levels) |>
  mutate(celltype = factor(celltype, levels = all_ct_levels))

# Display window for the log distance axis on the celltype plots. Start at the 1st
# percentile rather than the absolute minimum (a rare directly-adjacent cell that
# would leave a big empty stretch of low values), and end just past the max. This is
# applied via coord_flip(ylim=) so it zooms the view only -- every cell still feeds
# the violins and the mixed model. Raise the 0.01 quantile if you want to start higher.
dist_lo_ct <- as.numeric(quantile(md_dist_ct$dist_to_phf1_um, 0.01, na.rm = TRUE))
dist_hi_ct <- max(md_dist_ct$dist_to_phf1_um, na.rm = TRUE) * 1.10

# Per-donor x celltype mean distance (the dots / boxes drawn on the figure)
plot_b_df <- md_dist_ct |>
  group_by(sample_id, Braak, celltype) |>
  summarise(mean_dist = mean(dist_to_phf1_um), n_cells = n(), .groups = "drop") |>
  mutate(celltype = factor(celltype, levels = all_ct_levels))

# Cell-level linear mixed model across all celltypes; LOO contrasts (each celltype
# vs the average of all others), BH-adjusted. Guarded so the figure is produced
# even if the model fails to converge.
emm_loo_b <- NULL
m_full <- NULL
m_lrt  <- NULL
emm_pairs_b <- NULL
emm_b <- NULL
try({
  m_full <- lme4::lmer(dist_to_phf1_um ~ celltype + (1 | sample_id), data = md_dist_ct)
  m_null <- lme4::lmer(dist_to_phf1_um ~ 1 + (1 | sample_id), data = md_dist_ct)
  m_lrt  <- anova(m_null, m_full)   # lme4 refits with ML for the comparison

  emm_b <- emmeans(m_full, ~ celltype, lmer.df = "asymptotic")

  # Pairwise contrasts (Tukey) kept for the supplementary table only
  emm_pairs_b <- pairs(emm_b, adjust = "tukey") |>
    summary(infer = TRUE) |>
    as.data.frame()

  # Leave-one-out contrasts: each celltype vs the average of all OTHER celltypes.
  lv <- levels(md_dist_ct$celltype); kk <- length(lv)
  loo <- setNames(lapply(seq_len(kk), function(i) {
    w <- rep(-1 / (kk - 1), kk); w[i] <- 1; w
  }), lv)
  emm_loo_b <- contrast(emm_b, method = loo, adjust = "BH") |>
    summary(infer = TRUE) |>
    as.data.frame() |>
    dplyr::rename(celltype = contrast) |>
    mutate(
      celltype  = factor(celltype, levels = lv),
      sig_label = case_when(
        p.value < 0.0001 ~ "****",
        p.value < 0.001  ~ "***",
        p.value < 0.01   ~ "**",
        p.value < 0.05   ~ "*",
        TRUE             ~ ""
      ),
      # Direction of the LOO contrast (celltype mean minus the average of all
      # others): negative = this celltype sits CLOSER to PHF1+ neurons.
      direction = ifelse(estimate < 0, "closer", "further")
    ) |>
    arrange(celltype)
}, silent = FALSE)

# Asterisk colours encode the direction of the LOO contrast (no legend):
# red = this celltype sits CLOSER to PHF1+ neurons than
# the average celltype, blue = FURTHER. Both the jitter (celltype colour) and the
# asterisks (direction colour) feed a single identity colour scale, so they can
# use different colour mappings without conflicting.
col_closer  <- "#BD0026"   # red  — nearer to ptau pathology than average
col_further <- "#08519C"   # blue — further from ptau pathology than average

# Donor-level mean +/- SEM per celltype (crossbar over the cell-level violin)
donor_summ_b <- plot_b_df |>
  group_by(celltype) |>
  summarise(m = mean(mean_dist), sem = sd(mean_dist) / sqrt(n()), .groups = "drop") |>
  mutate(xc = as.integer(celltype))

# Best-of-both display: cell-level violin (spread of all PHF1-negative cells of
# each celltype, full range) behind per-donor mean dots + donor mean +/- SEM.
pb <- ggplot() +
  geom_violin(data = md_dist_ct, aes(x = celltype, y = dist_to_phf1_um, fill = celltype),
              scale = "width", linewidth = 0.2, alpha = 0.5, colour = NA) +
  geom_jitter(data = plot_b_df,
              aes(x = celltype, y = mean_dist,
                  colour = unname(celltype_palette[as.character(celltype)])),
              width = 0.12, height = 0, size = 0.9, alpha = 0.95) +
  geom_errorbar(data = donor_summ_b, aes(x = celltype, ymin = m - sem, ymax = m + sem),
                inherit.aes = FALSE, width = 0.16, linewidth = 0.3) +
  geom_segment(data = donor_summ_b, aes(x = xc - 0.2, xend = xc + 0.2, y = m, yend = m),
               inherit.aes = FALSE, linewidth = 0.4) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_y_log10(breaks = dist_breaks, labels = as.character(dist_breaks),
                expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = y_lab_dist) +
  coord_flip(ylim = c(dist_lo_ct, dist_hi_ct)) +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))

# Direction-coloured LOO asterisks, aligned in a row just above the tallest violin
if (!is.null(emm_loo_b)) {
  y_ast_b <- max(md_dist_ct$dist_to_phf1_um, na.rm = TRUE) * 1.03
  ast_b <- emm_loo_b |>
    select(celltype, sig_label, direction) |>
    mutate(y_position = y_ast_b,
           dir_col    = ifelse(direction == "closer", col_closer, col_further))
  pb <- pb + geom_text(data = ast_b,
                       aes(x = celltype, y = y_position, label = sig_label,
                           colour = dir_col),
                       inherit.aes = FALSE, size = 2.6, vjust = 0.5, hjust = 0.5,
                       fontface = "bold")
}

ggsave(file.path(out_dir, paste0("plot_", slug_b, ".pdf")), pb,
       width = 8.9, height = 6.5, units = "cm", device = "pdf")

write.table(plot_b_df, file.path(out_dir, paste0("source_data_", slug_b, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug_b, ".txt")))
cat("Panel (b) — mean distance to nearest PHF1+ neuron across celltypes\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Figure: cell-level violin (all PHF1-negative cells per celltype, full y-range)\n")
cat("behind per-donor mean dots + a donor mean +/- SEM crossbar. The violin shows\n")
cat("cell-level spread; dots/crossbar are the donor-level unit. Cell-level values:\n")
cat("source_data_dist_to_phf1_cells.tsv.\n")
cat("Stats: cell-level linear mixed model\n")
cat("  dist_to_phf1_um ~ celltype + (1 | sample_id)\n")
cat("with emmeans leave-one-out contrasts (each celltype vs the average of all\n")
cat("others, BH-adjusted) driving the asterisks above each box.\n\n")
cat("N cells (PHF1-negative):", nrow(md_dist_ct),
    "| N donors:", length(unique(md_dist_ct$sample_id)),
    "| N celltypes:", length(all_ct_levels), "\n\n")
if (!is.null(m_full)) {
  cat("Full model summary:\n"); print(summary(m_full))
  cat("\nNull vs full LRT:\n"); print(m_lrt)
  cat("\nEMMs (estimated marginal mean distance, um):\n"); print(summary(emm_b, infer = TRUE))
  cat("\nLeave-one-out contrasts (each celltype vs the AVERAGE OF ALL OTHER\n")
  cat("celltypes; BH-adjusted). These drive the asterisks on the figure.\n")
  print(emm_loo_b)
  cat("\nPairwise contrasts (Tukey-adjusted) -- supplementary table only.\n")
  print(emm_pairs_b)
} else {
  cat("(mixed model failed; see message log)\n")
}
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

if (!is.null(emm_loo_b))
  write.table(emm_loo_b, file.path(out_dir, paste0("stats_", slug_b, "_loo.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(emm_pairs_b))
  write.table(emm_pairs_b, file.path(out_dir, paste0("stats_", slug_b, "_emmeans.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

message("Done panel (b).")

# =============================================================================
# Panel (b) variant — identical data/stats, celltypes ordered by ascending distance
# =============================================================================
slug_b2 <- "dist_to_phf1_by_celltype_ordered"

# Order celltypes low -> high by the MEDIAN of the per-donor mean distances (the
# central tendency of the dots / crossbar).
ct_order_dist <- plot_b_df |>
  group_by(celltype) |>
  summarise(median_dist = median(mean_dist), .groups = "drop") |>
  arrange(median_dist) |>
  pull(celltype) |>
  as.character()

plot_b2_df <- plot_b_df |>
  mutate(celltype = factor(as.character(celltype), levels = ct_order_dist))
md_dist_ct2 <- md_dist_ct |>
  mutate(celltype = factor(as.character(celltype), levels = ct_order_dist))
donor_summ_b2 <- plot_b2_df |>
  group_by(celltype) |>
  summarise(m = mean(mean_dist), sem = sd(mean_dist) / sqrt(n()), .groups = "drop") |>
  mutate(xc = as.integer(celltype))

pb2 <- ggplot() +
  geom_violin(data = md_dist_ct2, aes(x = celltype, y = dist_to_phf1_um, fill = celltype),
              scale = "width", linewidth = 0.2, alpha = 0.5, colour = NA) +
  geom_jitter(data = plot_b2_df,
              aes(x = celltype, y = mean_dist,
                  colour = unname(celltype_palette[as.character(celltype)])),
              width = 0.12, height = 0, size = 0.9, alpha = 0.95) +
  geom_errorbar(data = donor_summ_b2, aes(x = celltype, ymin = m - sem, ymax = m + sem),
                inherit.aes = FALSE, width = 0.16, linewidth = 0.3) +
  geom_segment(data = donor_summ_b2, aes(x = xc - 0.2, xend = xc + 0.2, y = m, yend = m),
               inherit.aes = FALSE, linewidth = 0.4) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_y_log10(breaks = dist_breaks, labels = as.character(dist_breaks),
                expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = y_lab_dist) +
  coord_flip(ylim = c(dist_lo_ct, dist_hi_ct)) +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))

# Same direction-coloured LOO asterisks, aligned above the tallest violin
if (!is.null(emm_loo_b)) {
  y_ast_b2 <- max(md_dist_ct2$dist_to_phf1_um, na.rm = TRUE) * 1.03
  ast_b2 <- emm_loo_b |>
    select(celltype, sig_label, direction) |>
    mutate(celltype   = factor(as.character(celltype), levels = ct_order_dist),
           y_position = y_ast_b2,
           dir_col    = ifelse(direction == "closer", col_closer, col_further))
  pb2 <- pb2 + geom_text(data = ast_b2,
                         aes(x = celltype, y = y_position, label = sig_label,
                             colour = dir_col),
                         inherit.aes = FALSE, size = 2.6, vjust = 0.5, hjust = 0.5,
                         fontface = "bold")
}

ggsave(file.path(out_dir, paste0("plot_", slug_b2, ".pdf")), pb2,
       width = 8.9, height = 6.5, units = "cm", device = "pdf")

write.table(plot_b2_df, file.path(out_dir, paste0("source_data_", slug_b2, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug_b2, ".txt")))
cat("Panel (b) variant -- distance across celltypes, ordered low -> high mean distance\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Identical data and statistics to 'dist_to_phf1_by_celltype'; only the x-axis\n")
cat("ordering differs (celltypes sorted by ascending median of the per-donor mean\n")
cat("distances, rather than the canonical neuron->glia order). Figure composition\n")
cat("(violin + donor dots + mean +/- SEM crossbar) is identical to the canonical\n")
cat("panel. See stats_dist_to_phf1_by_celltype.txt for the full mixed model,\n")
cat("leave-one-out contrasts and Tukey pairwise tables.\n\n")
cat("Celltype ordering (low -> high median of per-donor mean distance, um):\n")
ord_tbl <- plot_b_df |>
  group_by(celltype) |>
  summarise(median_of_donor_means = median(mean_dist), .groups = "drop") |>
  arrange(median_of_donor_means)
print(as.data.frame(ord_tbl))
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

message("Done panel (b) ordered variant.")

# =============================================================================
# Panel (c) — distance to PHF1+ neuron, GLIA ONLY (no neurons / no unassigned),
#             celltypes ordered low -> high by median distance
# =============================================================================
slug_g <- "dist_to_phf1_by_glia"

# Glia celltypes present among the PHF1-negative cells (neurons and Unassigned/
# Unassigned Neuron excluded entirely).
glia_present <- intersect(glia_order, unique(md_dist$celltype))

md_dist_glia <- md_dist |>
  filter(celltype %in% glia_present) |>
  mutate(celltype = factor(celltype, levels = glia_present))

# Per-donor x glia mean distance (dots / crossbar)
plot_g_df <- md_dist_glia |>
  group_by(sample_id, Braak, celltype) |>
  summarise(mean_dist = mean(dist_to_phf1_um), n_cells = n(), .groups = "drop")

# Order glia low -> high by the median of the per-donor mean distances (matches the
# central tendency of the crossbar), then apply that ordering everywhere.
glia_levels <- plot_g_df |>
  group_by(celltype) |>
  summarise(median_dist = median(mean_dist), .groups = "drop") |>
  arrange(median_dist) |>
  pull(celltype) |>
  as.character()

md_dist_glia <- md_dist_glia |>
  mutate(celltype = factor(as.character(celltype), levels = glia_levels))
plot_g_df <- plot_g_df |>
  mutate(celltype = factor(as.character(celltype), levels = glia_levels))
donor_summ_g <- plot_g_df |>
  group_by(celltype) |>
  summarise(m = mean(mean_dist), sem = sd(mean_dist) / sqrt(n()), .groups = "drop") |>
  mutate(xc = as.integer(celltype))

# Log distance axis display window (glia-specific 1st pctile -> max)
dist_lo_g <- as.numeric(quantile(md_dist_glia$dist_to_phf1_um, 0.01, na.rm = TRUE))
dist_hi_g <- max(md_dist_glia$dist_to_phf1_um, na.rm = TRUE) * 1.10

# Glia-only mixed model + LOO contrasts (each glia vs the average of the OTHER glia)
emm_loo_g <- NULL; m_full_g <- NULL; m_lrt_g <- NULL; emm_pairs_g <- NULL; emm_g <- NULL
try({
  m_full_g <- lme4::lmer(dist_to_phf1_um ~ celltype + (1 | sample_id), data = md_dist_glia)
  m_null_g <- lme4::lmer(dist_to_phf1_um ~ 1 + (1 | sample_id), data = md_dist_glia)
  m_lrt_g  <- anova(m_null_g, m_full_g)
  emm_g <- emmeans(m_full_g, ~ celltype, lmer.df = "asymptotic")
  emm_pairs_g <- pairs(emm_g, adjust = "tukey") |> summary(infer = TRUE) |> as.data.frame()
  lv <- levels(md_dist_glia$celltype); kk <- length(lv)
  loo <- setNames(lapply(seq_len(kk), function(i) {
    w <- rep(-1 / (kk - 1), kk); w[i] <- 1; w
  }), lv)
  emm_loo_g <- contrast(emm_g, method = loo, adjust = "BH") |>
    summary(infer = TRUE) |>
    as.data.frame() |>
    dplyr::rename(celltype = contrast) |>
    mutate(
      celltype  = factor(celltype, levels = lv),
      sig_label = case_when(
        p.value < 0.0001 ~ "****",
        p.value < 0.001  ~ "***",
        p.value < 0.01   ~ "**",
        p.value < 0.05   ~ "*",
        TRUE             ~ ""
      ),
      direction = ifelse(estimate < 0, "closer", "further")
    ) |>
    arrange(celltype)
}, silent = FALSE)

pg <- ggplot() +
  geom_violin(data = md_dist_glia, aes(x = celltype, y = dist_to_phf1_um, fill = celltype),
              scale = "width", linewidth = 0.2, alpha = 0.5, colour = NA) +
  geom_jitter(data = plot_g_df,
              aes(x = celltype, y = mean_dist,
                  colour = unname(celltype_palette[as.character(celltype)])),
              width = 0.12, height = 0, size = 0.9, alpha = 0.95) +
  geom_errorbar(data = donor_summ_g, aes(x = celltype, ymin = m - sem, ymax = m + sem),
                inherit.aes = FALSE, width = 0.16, linewidth = 0.3) +
  geom_segment(data = donor_summ_g, aes(x = xc - 0.2, xend = xc + 0.2, y = m, yend = m),
               inherit.aes = FALSE, linewidth = 0.4) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  scale_colour_identity(guide = "none") +
  scale_y_log10(breaks = dist_breaks, labels = as.character(dist_breaks),
                expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = NULL, y = y_lab_dist) +
  coord_flip(ylim = c(dist_lo_g, dist_hi_g)) +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))

# Direction-coloured LOO asterisks (each glia vs the average of the other glia)
if (!is.null(emm_loo_g)) {
  y_ast_g <- max(md_dist_glia$dist_to_phf1_um, na.rm = TRUE) * 1.03
  ast_g <- emm_loo_g |>
    select(celltype, sig_label, direction) |>
    mutate(y_position = y_ast_g,
           dir_col    = ifelse(direction == "closer", col_closer, col_further))
  pg <- pg + geom_text(data = ast_g,
                       aes(x = celltype, y = y_position, label = sig_label,
                           colour = dir_col),
                       inherit.aes = FALSE, size = 2.6, vjust = 0.5, hjust = 0.5,
                       fontface = "bold")
}

ggsave(file.path(out_dir, paste0("plot_", slug_g, ".pdf")), pg,
       width = 6.6, height = 4.5, units = "cm", device = "pdf")

write.table(plot_g_df, file.path(out_dir, paste0("source_data_", slug_g, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug_g, ".txt")))
cat("Panel (c) -- distance to nearest PHF1+ neuron, GLIA ONLY\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Glia compared (ordered on figure low -> high median distance):",
    paste(glia_levels, collapse = ", "), "\n")
cat("Neurons and Unassigned/Unassigned Neuron are excluded entirely.\n\n")
cat("Figure: cell-level violin (all PHF1-negative glia per type, log axis) behind\n")
cat("per-donor mean dots + a donor mean +/- SEM crossbar. Cell-level values are in\n")
cat("source_data_dist_to_phf1_cells.tsv (subset to the glia celltypes).\n\n")
cat("Stats: cell-level linear mixed model\n")
cat("  dist_to_phf1_um ~ celltype + (1 | sample_id)   [fit on GLIA only]\n")
cat("with emmeans leave-one-out contrasts (each glia vs the average of the OTHER\n")
cat("glia, BH-adjusted) driving the asterisks; colour = direction (red = closer to\n")
cat("PHF1+ neurons than the other glia, blue = further).\n\n")
cat("N cells (PHF1-negative glia):", nrow(md_dist_glia),
    "| N donors:", length(unique(md_dist_glia$sample_id)),
    "| N glia types:", length(glia_levels), "\n\n")
if (!is.null(m_full_g)) {
  cat("Full model summary:\n"); print(summary(m_full_g))
  cat("\nNull vs full LRT:\n"); print(m_lrt_g)
  cat("\nEMMs (estimated marginal mean distance, um):\n"); print(summary(emm_g, infer = TRUE))
  cat("\nLeave-one-out contrasts (each glia vs the AVERAGE OF THE OTHER GLIA; BH).\n")
  print(emm_loo_g)
  cat("\nPairwise contrasts (Tukey-adjusted) -- supplementary.\n")
  print(emm_pairs_g)
} else {
  cat("(mixed model failed; see message log)\n")
}
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

if (!is.null(emm_loo_g))
  write.table(emm_loo_g, file.path(out_dir, paste0("stats_", slug_g, "_loo.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(emm_pairs_g))
  write.table(emm_pairs_g, file.path(out_dir, paste0("stats_", slug_g, "_emmeans.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

message("Done panel (c) glia-only.")
message("All outputs written to: ", out_dir)
