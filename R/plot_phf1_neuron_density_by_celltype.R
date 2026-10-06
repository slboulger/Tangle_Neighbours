# ---------------------------------------------------------------------------
# plot_phf1_neuron_density_by_celltype.R
#
# Figure panels: 2A
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# PHF1+ neuron DENSITY (counts per mm2) per neuron celltype.
# Distinct from the PROPORTION of each celltype that is PHF1+: here we show the
# absolute abundance of PHF1+
# (ptau/tangle-bearing) neurons per unit cortex, so the readout reflects tangle
# burden rather than how much tissue / how many cells each donor contributed.
#
# Denominator: per-donor tissue area = convex hull of that donor's cell centroids
# (x_slide_mm / y_slide_mm), in mm^2. No FOV-size assumption is needed; base R
# only. The hull uses ALL cells so the area exposure
# is identical for every celltype within a donor.
#
# Scope: the 9 named neuron subtypes (neuron_order; excludes "Unassigned Neuron").
# Layout: single pooled panel (all 9 donors); x = celltype, y = density, boxplot +
#   per-donor dots coloured by celltype.
#
# Stats: each donor contributes to every celltype -> donor-paired/nested. We fit a
#   negative-binomial GLMM on the counts with an area offset (density = rate):
#     n_phf1_pos ~ celltype + (1 | sample_id) + offset(log(area_mm2))
#   with emmeans leave-one-out contrasts (each celltype vs the average of all other
#   neuron subtypes, BH-adjusted) driving the asterisks. tryCatch-guarded so the
#   figure still renders if the GLMM fails to converge (counts are sparse).
#
# Output triple: plot PDF + source_data TSV + stats log (+ LOO TSV).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(MASS)      # glmer.nb lives in lme4 but pulls MASS::theta.ml
  library(emmeans)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")

out_dir <- "plots/phf1_neuron_density"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
slug <- "phf1_neuron_density_by_celltype"

# -------------------------------------------------------------------
# Load + slim frame (robust column resolution)
# -------------------------------------------------------------------
pick <- function(md, candidates, what, required = TRUE) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) {
    msg <- sprintf("Could not find a column for '%s'. Tried: %s.", what,
                   paste(candidates, collapse = ", "))
    if (required) stop(msg) else { warning(msg); return(NA_character_) }
  }
  hit[1]
}

seu <- readRDS("seu_PHF1.rds")
meta <- seu@meta.data

col_samp <- pick(meta, c("sample_id", "sampleID"), "sample_id")
col_ct   <- pick(meta, c("celltype"), "celltype")
col_brk  <- pick(meta, c("Braak"), "Braak")
col_phf1 <- pick(meta, c("PHF1"), "PHF1")
col_x    <- pick(meta, c("x_slide_mm", "sdimx"), "x_slide_mm")
col_y    <- pick(meta, c("y_slide_mm", "sdimy"), "y_slide_mm")

md <- tibble(
  sample_id = as.character(meta[[col_samp]]),
  Braak     = factor(as.character(meta[[col_brk]]), levels = braak_levels),
  celltype  = as.character(meta[[col_ct]]),
  PHF1      = meta[[col_phf1]] %in% c(TRUE, "TRUE", "True", 1, "1"),
  x_mm      = as.numeric(meta[[col_x]]),
  y_mm      = as.numeric(meta[[col_y]])
) %>%
  dplyr::filter(!is.na(sample_id), !is.na(celltype),
                is.finite(x_mm), is.finite(y_mm))

message("N cells (finite coords): ", nrow(md))
message("table(sample_id, Braak):"); print(table(md$sample_id, md$Braak))
message("table(celltype, PHF1) [neurons]:")
print(table(md$celltype[md$celltype %in% neuron_order],
            md$PHF1[md$celltype %in% neuron_order]))

# -------------------------------------------------------------------
# Per-donor tissue area (mm^2) via convex hull of ALL cell centroids
# -------------------------------------------------------------------
# Shoelace polygon area; chull returns hull vertices in order.
poly_area <- function(x, y) {
  n <- length(x)
  abs(sum(x * y[c(2:n, 1)] - x[c(2:n, 1)] * y)) / 2
}

hull_area <- md %>%
  group_by(sample_id, Braak) %>%
  summarise(
    area_mm2      = { i <- chull(x_mm, y_mm); poly_area(x_mm[i], y_mm[i]) },
    n_cells_total = n(),
    .groups       = "drop"
  )

# -------------------------------------------------------------------
# PHF1+ neuron counts -> density per donor x celltype (zeros kept)
# -------------------------------------------------------------------
counts <- md %>%
  dplyr::filter(celltype %in% neuron_order, PHF1) %>%
  dplyr::count(sample_id, celltype, name = "n_phf1_pos") %>%
  tidyr::complete(sample_id = unique(md$sample_id), celltype = neuron_order,
                  fill = list(n_phf1_pos = 0)) %>%
  left_join(hull_area, by = "sample_id") %>%
  mutate(
    celltype    = factor(celltype, levels = neuron_order),
    density_mm2 = n_phf1_pos / area_mm2
  ) %>%
  arrange(celltype, sample_id)

# -------------------------------------------------------------------
# Stats: negative-binomial GLMM on counts with area offset (guarded)
# -------------------------------------------------------------------
emm_loo <- tryCatch({
  m_full <- lme4::glmer.nb(
    n_phf1_pos ~ celltype + (1 | sample_id) + offset(log(area_mm2)),
    data = counts,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  m_null <- lme4::glmer.nb(
    n_phf1_pos ~ 1 + (1 | sample_id) + offset(log(area_mm2)),
    data = counts,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  m_lrt <- anova(m_null, m_full)
  emm   <- emmeans(m_full, ~ celltype)

  # Leave-one-out: each celltype vs the average of all OTHER neuron subtypes.
  lv <- levels(counts$celltype); k <- length(lv)
  loo <- setNames(lapply(seq_len(k), function(i) {
    w <- rep(-1 / (k - 1), k); w[i] <- 1; w
  }), lv)

  loo_tab <- contrast(emm, method = loo, adjust = "BH") %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    dplyr::rename(celltype = contrast) %>%
    mutate(
      celltype  = factor(celltype, levels = lv),
      sig_label = case_when(
        p.value < 0.0001 ~ "****",
        p.value < 0.001  ~ "***",
        p.value < 0.01   ~ "**",
        p.value < 0.05   ~ "*",
        TRUE             ~ ""
      )
    ) %>%
    arrange(celltype)

  list(m_full = m_full, m_lrt = m_lrt, emm = emm, loo_tab = loo_tab)
}, error = function(e) {
  message("NB-GLMM failed (", conditionMessage(e), "); figure is descriptive only.")
  NULL
})

# -------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------
# Horizontal layout: density on x, neuron celltypes on y (canonical order top->bottom).
p <- ggplot(counts, aes(x = density_mm2, y = celltype)) +
  geom_boxplot(aes(fill = celltype), alpha = 0.55, outlier.shape = NA,
               linewidth = 0.3, width = 0.6) +
  geom_jitter(aes(colour = celltype), width = 0, height = 0.12,
              size = 0.7, alpha = 0.95) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  scale_colour_manual(values = celltype_palette, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_discrete(limits = rev(neuron_order)) +
  labs(x = expression("PHF1+ neurons per " * mm^2), y = NULL) +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x  = element_text(angle = 0, hjust = 0.5, vjust = 1),
        plot.margin  = margin(t = 4, r = 10, b = 4, l = 4, unit = "pt"))

# Asterisks to the right of each row (only if the GLMM fit)
if (!is.null(emm_loo)) {
  x_top <- max(counts$density_mm2, na.rm = TRUE)
  x_pad <- max(x_top * 0.05, 1e-4)
  per_ct_max <- counts %>%
    group_by(celltype) %>%
    summarise(box_top = max(density_mm2, na.rm = TRUE), .groups = "drop")
  asterisk_df <- emm_loo$loo_tab %>%
    dplyr::select(celltype, sig_label) %>%
    left_join(per_ct_max, by = "celltype") %>%
    mutate(x_position = box_top + x_pad)
  p <- p + geom_text(data = asterisk_df,
                     aes(x = x_position, y = celltype, label = sig_label),
                     inherit.aes = FALSE, size = 2.6, hjust = 0)
}

ggsave(file.path(out_dir, paste0("plot_", slug, ".pdf")), p,
       width = 3.6, height = 2.0, units = "in", device = "pdf")

# -------------------------------------------------------------------
# Source data: one row per donor x celltype, exactly as drawn
# -------------------------------------------------------------------
src <- counts %>%
  dplyr::select(sample_id, Braak, celltype, n_phf1_pos, area_mm2,
                n_cells_total, density_mm2)
write.table(src, file.path(out_dir, paste0("source_data_", slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------------------------------------------
# Stats log
# -------------------------------------------------------------------
sink(file.path(out_dir, paste0("stats_", slug, ".txt")))
cat("PHF1+ neuron density (counts per mm^2) per neuron celltype\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Denominator: per-donor tissue area = convex hull of ALL cell centroids\n")
cat("(x_slide_mm / y_slide_mm), shoelace area in mm^2. Density = n_phf1_pos / area.\n")
cat("Only the 9 named neuron subtypes (neuron_order) are shown; Unassigned Neuron\n")
cat("is excluded. Missing donor x celltype combinations are true zeros and kept.\n\n")

cat("N cells (finite coords):", nrow(md), "\n")
cat("N donors:", length(unique(counts$sample_id)), "\n\n")

cat("Per-donor tissue area and total PHF1+ neurons:\n")
donor_tab <- counts %>%
  group_by(sample_id, Braak, area_mm2, n_cells_total) %>%
  summarise(n_phf1_pos_total = sum(n_phf1_pos), .groups = "drop") %>%
  arrange(Braak, sample_id)
print(as.data.frame(donor_tab))
cat("\n")

cat("Stats: cell-group comparison is donor-paired (each donor contributes to every\n")
cat("celltype), so we fit a negative-binomial GLMM on the counts with an area offset:\n")
cat("  n_phf1_pos ~ celltype + (1 | sample_id) + offset(log(area_mm2))\n")
cat("emmeans leave-one-out contrasts (each celltype vs the average of all other\n")
cat("neuron subtypes, BH-adjusted) drive the figure asterisks.\n\n")

cat("NOTES:\n")
cat("- Per-donor PHF1+ neuron totals span single digits to ~180, so many celltype x donor\n")
cat("  densities are 0 (kept as true zeros).\n")
cat("- The convex hull is a consistent, assumption-light denominator rather than an exact\n")
cat("  tissue mask (it includes any concavities of the section).\n\n")

if (!is.null(emm_loo)) {
  cat("Full NB-GLMM summary:\n");        print(summary(emm_loo$m_full))
  cat("\nNull vs full LRT:\n");          print(emm_loo$m_lrt)
  cat("\nEMMs (log rate scale):\n");     print(summary(emm_loo$emm, infer = TRUE))
  cat("\nLeave-one-out contrasts (BH-adjusted; drive the asterisks):\n")
  print(emm_loo$loo_tab)
} else {
  cat("NB-GLMM did not fit; figure is descriptive only (no asterisks).\n")
}

cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

if (!is.null(emm_loo)) {
  write.table(emm_loo$loo_tab,
              file.path(out_dir, paste0("stats_", slug, "_loo.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("All outputs written to: ", out_dir)
