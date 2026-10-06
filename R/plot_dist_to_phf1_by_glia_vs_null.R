#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_dist_to_phf1_by_glia_vs_null.R
#
# Figure panels: 6A
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_dist_to_phf1_by_glia_vs_null.R
#
# OBSERVED-MINUS-NULL variant of plot_dist_to_phf1_by_glia.pdf.
#
# Plot 1 geometry is byte-for-byte the canonical glia panel (Panel (c) of
# R/plot_dist_to_phf1_braak_celltype.R, lines 438-546): observed cell-level
# violins, observed per-donor mean dots, donor mean +/- SEM crossbar, log10
# distance axis, coord_flip. ONLY the annotated stats change.
#
# THE STATISTIC IS WITHIN-CELLTYPE. For each glial type, the mean distance to
# the nearest PHF1+ neuron is compared with the SAME quantity computed under
# each of the 1000 PHF1-label shuffles. Nothing is compared to the other glia:
#
#   T_obs(g)  = mean over donors of that donor's mean distance for glia g
#   T_null(g,b) = same quantity in null draw b
#   delta(g)  = T_obs(g) - median_b T_null(g,b)
#   two-sided empirical p = (1 + #{ |T_null - M| >= |delta| }) / (B + 1)
#
# BH-adjusted across the 6 glia. delta < 0 means that glial type sits closer to
# REAL PHF1+ neurons than to randomly relabelled ones.
#
# The canonical panel instead tests a leave-one-out contrast (each glia vs the
# average of the OTHER glia) against zero, a relative rather than an absolute
# question. The canonical LOO result is carried in the stats files for comparison.
#
# Inputs -- the precompute (run R/null_glia_distance_within_celltype.R first)
# plus already-written canonical outputs. No Seurat load.
#
# Outputs, three triples in plots/dist_to_phf1_glia_vs_null/:
#   plot_dist_to_phf1_by_glia_vs_null   canonical violin panel, null-corrected asterisks
#   plot_dist_to_phf1_glia_null_vs_obs  observed mean inside its null distribution
#   plot_dist_to_phf1_glia_delta_null   the deviation itself, forest of delta

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")

out_dir <- "plots/dist_to_phf1_glia_vs_null"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
slug <- "dist_to_phf1_by_glia_vs_null"

y_lab_dist  <- expression("Distance to PHF1+ neuron (" * mu * "m)")
dist_breaks <- c(1, 3, 10, 30, 100, 300, 1000, 3000)
# Direction colours kept identical to the canonical glia panel for cross-panel
# consistency (red = closer, blue = further).
col_closer  <- "#BD0026"
col_further <- "#08519C"

f_cells <- "plots/dist_to_phf1/source_data_dist_to_phf1_cells.tsv"
f_dots  <- "plots/dist_to_phf1/source_data_dist_to_phf1_by_glia.tsv"
f_loo   <- "plots/dist_to_phf1/stats_dist_to_phf1_by_glia_loo.tsv"
f_draw  <- file.path(out_dir, "null_glia_within_celltype_per_draw.tsv")
f_donor <- file.path(out_dir, "null_glia_within_celltype_donor_draw.tsv")
for (f in c(f_cells, f_dots, f_loo, f_draw, f_donor))
  if (!file.exists(f)) stop("Missing input: ", f,
                            if (grepl("within_celltype", f))
                              "  -- run R/null_glia_distance_within_celltype.R first" else "")

## ---------------------------------------------------------------------------
## Observed data for the plot geometry.
## ---------------------------------------------------------------------------
plot_g_df <- read.delim(f_dots, stringsAsFactors = FALSE) |>
  mutate(sample_id = as.character(sample_id), celltype = as.character(celltype))

# Ordering: ascending median of per-donor means, exactly as the canonical panel.
glia_levels <- plot_g_df |>
  group_by(celltype) |>
  summarise(median_dist = median(mean_dist), .groups = "drop") |>
  arrange(median_dist) |>
  pull(celltype)

md_dist_glia <- read.delim(f_cells, stringsAsFactors = FALSE) |>
  filter(celltype %in% glia_levels) |>
  mutate(celltype = factor(celltype, levels = glia_levels))

plot_g_df <- plot_g_df |> mutate(celltype = factor(celltype, levels = glia_levels))

donor_summ_g <- plot_g_df |>
  group_by(celltype) |>
  summarise(m = mean(mean_dist), sem = sd(mean_dist) / sqrt(n()), .groups = "drop") |>
  mutate(xc = as.integer(celltype))

dist_lo_g <- as.numeric(quantile(md_dist_glia$dist_to_phf1_um, 0.01, na.rm = TRUE))
dist_hi_g <- max(md_dist_glia$dist_to_phf1_um, na.rm = TRUE) * 1.10

## ---------------------------------------------------------------------------
## Within-celltype observed vs null.
## Primary statistic: mean of the per-donor mean distances (what the canonical
## panel's crossbar draws). pooled_cell_mean is carried as a secondary check --
## it lets donors with more cells of a type dominate.
## ---------------------------------------------------------------------------
per_draw <- read.delim(f_draw, stringsAsFactors = FALSE) |>
  mutate(celltype = as.character(celltype)) |>
  filter(celltype %in% glia_levels)

STAT <- "mean_of_donor_means"
per_draw$stat <- per_draw[[STAT]]

obs_stat  <- per_draw |> filter(series == "observed") |> select(celltype, T_obs = stat)
null_stat <- per_draw |> filter(series == "null")     |> select(celltype, draw, T_null = stat)

B_per_ct <- null_stat |> count(celltype, name = "B")
if (length(unique(B_per_ct$B)) != 1L)
  warning("Unequal null draw counts per celltype: ", paste(B_per_ct$B, collapse = ", "))

# Observed SEM across the 9 donors, for the observed marker in plot 2.
obs_sem <- plot_g_df |>
  group_by(celltype) |>
  summarise(obs_sem = sd(mean_dist) / sqrt(n()), n_donors = n(), .groups = "drop") |>
  mutate(celltype = as.character(celltype))

cmp <- null_stat |>
  group_by(celltype) |>
  summarise(B         = n(),
            null_med  = median(T_null),
            null_mean = mean(T_null),
            null_sd   = sd(T_null),
            null_q025 = quantile(T_null, 0.025),
            null_q975 = quantile(T_null, 0.975),
            .groups = "drop") |>
  left_join(obs_stat, by = "celltype") |>
  left_join(obs_sem,  by = "celltype")

# Two-sided empirical p: how extreme is the observed statistic relative to its
# null distribution, both measured as a deviation from the null median.
cmp$p_emp <- vapply(seq_len(nrow(cmp)), function(i) {
  nv <- null_stat$T_null[null_stat$celltype == cmp$celltype[i]]
  d  <- abs(cmp$T_obs[i] - cmp$null_med[i])
  (1 + sum(abs(nv - cmp$null_med[i]) >= d)) / (length(nv) + 1)
}, numeric(1))

cmp <- cmp |>
  mutate(
    delta = T_obs - null_med,
    se_nullmed = 1.2533 * null_sd / sqrt(B),   # Monte-Carlo SE of the null median
    # CI on delta is calibrated on the PERMUTATION null: it must be consistent
    # with p_emp, so it uses the spread of the null statistic, not a model SE.
    se_delta = sqrt(null_sd^2 + se_nullmed^2),
    delta_lo = delta - 1.96 * se_delta,
    delta_hi = delta + 1.96 * se_delta,
    z_null    = delta / null_sd,               # in null SDs
    delta_pct = 100 * delta / null_med,        # % of the null expected distance
    p_emp_bh  = p.adjust(p_emp, method = "BH"),
    sig_label = case_when(
      p_emp_bh < 0.0001 ~ "****",
      p_emp_bh < 0.001  ~ "***",
      p_emp_bh < 0.01   ~ "**",
      p_emp_bh < 0.05   ~ "*",
      TRUE              ~ ""),
    direction = ifelse(delta < 0, "closer", "further"),
    dir_col   = ifelse(direction == "closer", col_closer, col_further),
    celltype  = factor(celltype, levels = glia_levels)) |>
  arrange(celltype)

## Donor-level paired effect size. Each donor contributes an observed mean and a
## null mean (averaged over draws) for the same celltype, so the pair is donor-
## matched: paired Cohen's dz plus the count of donors moving the same way.
donor_draw <- read.delim(f_donor, stringsAsFactors = FALSE) |>
  mutate(celltype = as.character(celltype), sample_id = as.character(sample_id)) |>
  filter(celltype %in% glia_levels)

donor_pair <- donor_draw |>
  group_by(celltype, sample_id) |>
  summarise(obs  = mean_dist[series == "observed"],
            nul  = mean(mean_dist[series == "null"]),
            .groups = "drop") |>
  mutate(diff = obs - nul)

donor_eff <- donor_pair |>
  group_by(celltype) |>
  summarise(n_donors      = n(),
            mean_diff     = mean(diff),
            sd_diff       = sd(diff),
            dz            = mean(diff) / sd(diff),
            n_same_dir    = sum(sign(diff) == sign(mean(diff))),
            p_wilcox      = suppressWarnings(
                              wilcox.test(diff, mu = 0, exact = FALSE)$p.value),
            .groups = "drop") |>
  mutate(p_wilcox_bh = p.adjust(p_wilcox, method = "BH"),
         celltype    = factor(celltype, levels = glia_levels)) |>
  arrange(celltype)

cmp <- cmp |> left_join(donor_eff |> select(-n_donors), by = "celltype")

# Canonical leave-one-out result, for context only.
loo_ctx <- read.delim(f_loo, stringsAsFactors = FALSE) |>
  transmute(celltype = as.character(celltype),
            loo_est = estimate, loo_p_vs_zero = p.value)
cmp <- cmp |> mutate(celltype_chr = as.character(celltype)) |>
  left_join(loo_ctx, by = c("celltype_chr" = "celltype")) |>
  select(-celltype_chr)

# Secondary statistic (pooled cell-level mean), for the stats log.
per_draw$stat2 <- per_draw$pooled_cell_mean
cmp2 <- per_draw |>
  group_by(celltype) |>
  summarise(T_obs_pooled  = stat2[series == "observed"],
            null_med_pool = median(stat2[series == "null"]),
            null_sd_pool  = sd(stat2[series == "null"]),
            .groups = "drop") |>
  mutate(delta_pooled = T_obs_pooled - null_med_pool,
         z_pooled     = delta_pooled / null_sd_pool,
         celltype     = factor(celltype, levels = glia_levels)) |>
  arrange(celltype)

## ---------------------------------------------------------------------------
## Shared summary tables (printed into all three stats files).
## ---------------------------------------------------------------------------
tab1 <- cmp |>
  transmute(celltype,
            obs_um   = round(T_obs, 1),
            null_med = round(null_med, 1),
            null_sd  = round(null_sd, 2),
            null_95  = paste0("[", round(null_q025, 1), ", ", round(null_q975, 1), "]"),
            delta    = round(delta, 2),
            delta_CI = paste0("[", round(delta_lo, 1), ", ", round(delta_hi, 1), "]"),
            p_emp    = signif(p_emp, 3),
            p_emp_BH = signif(p_emp_bh, 3),
            sig      = sig_label,
            direction)

tab2 <- cmp |>
  transmute(celltype,
            delta_um  = round(delta, 2),
            delta_CI  = paste0("[", round(delta_lo, 1), ", ", round(delta_hi, 1), "]"),
            delta_pct = round(delta_pct, 2),
            z_null_SD = round(z_null, 2),
            donor_dz  = round(dz, 2),
            donor_diff_um = round(mean_diff, 2),
            n_same_dir = paste0(n_same_dir, "/", n_donors),
            p_wilcox_BH = signif(p_wilcox_bh, 3))

tab3 <- cmp |>
  transmute(celltype,
            canonical_LOO_est = round(loo_est, 2),
            canonical_p       = signif(loo_p_vs_zero, 3),
            canonical_sig     = case_when(loo_p_vs_zero < 0.0001 ~ "****",
                                          loo_p_vs_zero < 0.001  ~ "***",
                                          loo_p_vs_zero < 0.01   ~ "**",
                                          loo_p_vs_zero < 0.05   ~ "*", TRUE ~ ""),
            within_ct_delta   = round(delta, 2),
            within_ct_p_BH    = signif(p_emp_bh, 3),
            within_ct_sig     = sig_label,
            within_ct_dir     = direction)

hdr_method <- function() {
  cat("Statistic -- WITHIN CELLTYPE, observed vs null. Nothing is compared to the\n")
  cat("other glia:\n")
  cat("  T_obs(g)    = mean over the 9 donors of that donor's mean distance to the\n")
  cat("                nearest PHF1+ neuron, for glial type g (PHF1-negative cells).\n")
  cat("                This is exactly what the crossbar in the canonical panel draws.\n")
  cat("  T_null(g,b) = the same quantity in null draw b.\n")
  cat("  M_g         = median of the B null values.\n")
  cat("  delta(g)    = T_obs(g) - M_g   (um). delta < 0 = that glial type sits CLOSER\n")
  cat("                to REAL PHF1+ neurons than to randomly relabelled ones.\n")
  cat("  p_emp       = (1 + #{ |T_null - M_g| >= |delta| }) / (B + 1)   [two-sided]\n")
  cat("  BH-adjusted across the 6 glia.\n\n")
  cat("Null scheme:\n")
  cat("  PHF1+ neuron labels shuffled within (celltype x sample), so the number of\n")
  cat("  PHF1+ neurons per donor is preserved; distance to the nearest PHF1+ neuron\n")
  cat("  recomputed per draw. Precompute: R/null_glia_distance_within_celltype.R.\n\n")
}

## ---------------------------------------------------------------------------
## PLOT 1 -- canonical glia-panel layout, within-celltype null-corrected asterisks.
## ---------------------------------------------------------------------------
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

y_ast <- max(md_dist_glia$dist_to_phf1_um, na.rm = TRUE) * 1.03
ast <- cmp |> filter(sig_label != "") |> mutate(y_position = y_ast)
if (nrow(ast) > 0) {
  pg <- pg + geom_text(data = ast,
                       aes(x = celltype, y = y_position, label = sig_label, colour = dir_col),
                       inherit.aes = FALSE, size = 2.6, vjust = 0.5, hjust = 0.5,
                       fontface = "bold")
}

ggsave(file.path(out_dir, paste0("plot_", slug, ".pdf")), pg,
       width = 6.6, height = 4.5, units = "cm", device = "pdf")

write.table(plot_g_df, file.path(out_dir, paste0("source_data_", slug, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

eff_cols <- c("celltype", "T_obs", "obs_sem", "null_med", "null_mean", "null_sd",
              "null_q025", "null_q975", "B", "delta", "se_delta", "delta_lo", "delta_hi",
              "delta_pct", "z_null", "mean_diff", "sd_diff", "dz", "n_same_dir",
              "n_donors", "p_wilcox", "p_wilcox_bh", "p_emp", "p_emp_bh", "sig_label",
              "direction", "loo_est", "loo_p_vs_zero")
write.table(cmp[, eff_cols],
            file.path(out_dir, paste0("stats_", slug, "_effectsize.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug, ".txt")))
cat("Distance to nearest PHF1+ neuron, GLIA ONLY -- WITHIN-CELLTYPE observed vs null\n")
cat("==============================================================================\n\n")
cat("Figure geometry is identical to plots/dist_to_phf1/plot_dist_to_phf1_by_glia.pdf\n")
cat("(observed cell-level violins, observed per-donor mean dots, donor mean +/- SEM\n")
cat("crossbar, log10 axis, glia ordered low -> high median of per-donor means).\n")
cat("ONLY the annotated significance differs.\n\n")
hdr_method()
cat("Colour of the asterisks: red = closer to real PHF1+ neurons than to shuffled\n")
cat("ones (delta < 0), blue = further (delta > 0).\n\n")

cat("Effect sizes reported per glia:\n")
cat("  delta with 95% CI -- unstandardised, um. The CI is calibrated on the PERMUTATION\n")
cat("    null: delta +/- 1.96*sqrt(sd_null^2 + se_nullmed^2), sd_null = SD of the B null\n")
cat("    values, se_nullmed = 1.2533*sd_null/sqrt(B) = Monte-Carlo SE of the null median.\n")
cat("    A model SE is deliberately NOT used: it would not be consistent with p_emp.\n")
cat("  delta_pct         -- delta as % of the null expected distance M_g.\n")
cat("  z_null            -- delta / sd_null: how many null SDs from the null centre.\n")
cat("  donor_dz          -- paired Cohen's dz over the 9 donors (mean paired difference\n")
cat("    observed-minus-null / SD of those differences), with n_same_dir = donors moving\n")
cat("    the same way and a Wilcoxon signed-rank p.\n\n")

cat("N cells (PHF1-negative glia):", nrow(md_dist_glia),
    "| N donors:", length(unique(md_dist_glia$sample_id)),
    "| N glia types:", nlevels(cmp$celltype),
    "| null draws B:", unique(cmp$B), "\n")
cat("Resolution limit: p_emp cannot fall below 1/(B+1) =",
    signif(1 / (unique(cmp$B)[1] + 1), 3),
    "-- '****' is unreachable with B =", unique(cmp$B)[1], "draws.\n\n")
cat("Glia order on the figure (low -> high median distance):",
    paste(glia_levels, collapse = ", "), "\n\n")

cat("=== Within-celltype observed vs null (primary: mean of per-donor means) ===\n")
print(as.data.frame(tab1), row.names = FALSE)

cat("\n=== Effect sizes ===\n")
print(as.data.frame(tab2), row.names = FALSE)

cat("\n=== Per-donor paired differences (observed - null, um) ===\n")
print(as.data.frame(donor_pair |>
  mutate(across(c(obs, nul, diff), \(x) round(x, 1))) |>
  pivot_wider(names_from = celltype, values_from = c(obs, nul, diff)) |>
  select(sample_id, starts_with("diff_"))), row.names = FALSE)

cat("\n=== Secondary statistic: pooled cell-level mean (donors weighted by n cells) ===\n")
print(as.data.frame(cmp2 |> mutate(across(where(is.numeric), \(x) round(x, 2)))),
      row.names = FALSE)
cat("Reported as a donor-weighting check: agreement in sign with the primary\n")
cat("statistic is the check.\n")

cat("\n=== Context: what the canonical panel tested ===\n")
print(as.data.frame(tab3), row.names = FALSE)
cat("\nThe canonical panel's asterisks come from a leave-one-out contrast (each glia\n")
cat("vs the average of the OTHER glia) tested against zero. That is a different\n")
cat("question (relative to the other glia rather than to the null); it is shown\n")
cat("here for comparison and is not what this figure annotates.\n")

cat("\nTabular output: stats_", slug, "_effectsize.tsv\n", sep = "")
cat("Precompute inputs:\n  ", f_draw, "\n  ", f_donor, "\n", sep = "")
cat("\n------------------------------------------------------------\n")
print(sessionInfo())
sink()

## ===========================================================================
## PLOT 2 -- where the observed mean sits inside its null distribution.
## ===========================================================================
slug_n <- "dist_to_phf1_glia_null_vs_obs"

null_plot <- null_stat |> mutate(celltype = factor(celltype, levels = glia_levels))

obs_pt <- cmp |>
  mutate(obs_lo = T_obs - obs_sem,
         obs_hi = T_obs + obs_sem,
         xc     = as.integer(celltype))

pn <- ggplot() +
  geom_violin(data = null_plot, aes(x = celltype, y = T_null),
              scale = "width", linewidth = 0.2, colour = NA, fill = "grey70", alpha = 0.75) +
  geom_segment(data = obs_pt, aes(x = xc - 0.22, xend = xc + 0.22,
                                  y = null_med, yend = null_med),
               inherit.aes = FALSE, linewidth = 0.35, colour = "grey25") +
  geom_linerange(data = obs_pt, aes(x = celltype, ymin = obs_lo, ymax = obs_hi,
                                    colour = dir_col),
                 inherit.aes = FALSE, linewidth = 0.4) +
  geom_point(data = obs_pt, aes(x = celltype, y = T_obs, colour = dir_col),
             inherit.aes = FALSE, size = 1.5, shape = 18) +
  scale_colour_identity(guide = "none") +
  labs(x = NULL, y = expression("Mean distance to PHF1+ neuron (" * mu * "m)")) +
  coord_flip() +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))

ggsave(file.path(out_dir, paste0("plot_", slug_n, ".pdf")), pn,
       width = 6.6, height = 4.5, units = "cm", device = "pdf")

write.table(
  bind_rows(
    null_plot |> transmute(celltype = as.character(celltype), series = "null_draw",
                           draw, value = T_null, lo = NA_real_, hi = NA_real_),
    obs_pt |> transmute(celltype = as.character(celltype), series = "observed",
                        draw = NA_character_, value = T_obs, lo = obs_lo, hi = obs_hi),
    obs_pt |> transmute(celltype = as.character(celltype), series = "null_median",
                        draw = NA_character_, value = null_med,
                        lo = null_q025, hi = null_q975)),
  file.path(out_dir, paste0("source_data_", slug_n, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug_n, ".txt")))
cat("Observed glia mean distance inside its null-label distribution\n")
cat("=============================================================\n\n")
cat("Grey violin = the", unique(cmp$B), "null values of T_null(g,b) for that glia (density,\n")
cat("scale = width). Dark crossbar = the null median M_g. Diamond = the OBSERVED\n")
cat("T_obs(g) with +/- 1 SEM across the 9 donors, coloured red (closer to real PHF1+\n")
cat("neurons than to shuffled ones) or blue (further).\n\n")
hdr_method()
cat("=== Within-celltype observed vs null (primary: mean of per-donor means) ===\n")
print(as.data.frame(tab1), row.names = FALSE)
cat("\n=== Effect sizes ===\n")
print(as.data.frame(tab2), row.names = FALSE)
cat("\nFull detail, per-donor differences and the secondary pooled statistic:\n")
cat("  stats_", slug, ".txt\n", sep = "")
cat("\n------------------------------------------------------------\n")
print(sessionInfo())
sink()

## ===========================================================================
## PLOT 3 -- the deviation itself: delta = observed - null median, with its CI.
## ===========================================================================
slug_d <- "dist_to_phf1_glia_delta_null"

# Ordered by the deviation itself (most negative at the bottom), NOT by the
# canonical distance ordering -- this panel is about the ranking of delta.
# ggplot puts factor level 1 nearest the axis origin, which under coord_flip is
# the bottom, so ascending delta gives most-negative-at-bottom.
delta_levels <- cmp |> arrange(delta) |> pull(celltype) |> as.character()

# Significance is carried by COLOUR ALONE (no asterisks): black = not significant
# vs the null, red = significantly closer, blue = significantly further
# (BH-adjusted empirical p < 0.05).
dev_df <- cmp |>
  mutate(celltype = factor(as.character(celltype), levels = delta_levels),
         is_sig   = p_emp_bh < 0.05,
         pt_col   = ifelse(is_sig, dir_col, "black")) |>
  arrange(celltype)

x_pad <- 0.06 * diff(range(c(dev_df$delta_lo, dev_df$delta_hi)))

# Plain forest: the only horizontal line is the 95% CI. No lollipop stem to zero --
# a stem would visually merge with the CI and make every interval look as if it
# touched the no-deviation line.
pd <- ggplot(dev_df, aes(x = celltype, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  geom_linerange(aes(ymin = delta_lo, ymax = delta_hi, colour = pt_col), linewidth = 0.4) +
  geom_point(aes(colour = pt_col), size = 1.5) +
  scale_colour_identity(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  labs(x = NULL,
       y = expression("Distance to PHF1+ neuron, observed - null (" * mu * "m)")) +
  coord_flip(ylim = c(min(dev_df$delta_lo) - x_pad, max(dev_df$delta_hi) + x_pad)) +
  theme_classic(base_size = 8) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        # 45-character title at the standard 8pt overruns the 6.6 cm panel. 7pt
        # keeps it on one line (atop() leaves an ugly interline gap) and matches
        # the tick-label size, so it still reads as part of the same type scale.
        axis.title.x = element_text(size = 7),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 4, unit = "pt"))

ggsave(file.path(out_dir, paste0("plot_", slug_d, ".pdf")), pd,
       width = 6.6, height = 4.5, units = "cm", device = "pdf")

write.table(dev_df[, c(eff_cols, "is_sig", "pt_col")],
            file.path(out_dir, paste0("source_data_", slug_d, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(out_dir, paste0("stats_", slug_d, ".txt")))
cat("Deviation of the within-celltype glia mean distance from its null-label null\n")
cat("===========================================================================\n\n")
cat("Point = delta(g) = T_obs(g) - M_g. Line = 95% CI on delta, calibrated on the\n")
cat("permutation null (1.96*sqrt(sd_null^2 + se_nullmed^2)). Dashed line = no\n")
cat("deviation from the null, i.e. that glial type is no closer to real PHF1+ neurons\n")
cat("than to randomly relabelled ones.\n\n")
cat("Encoding: significance is carried by COLOUR ALONE -- there are no asterisks.\n")
cat("  black = NOT significant vs the null (BH-adjusted empirical p >= 0.05)\n")
cat("  red   = significantly CLOSER to real PHF1+ neurons than to shuffled ones\n")
cat("  blue  = significantly FURTHER\n\n")
cat("Ordering: glia are ordered by delta (most negative at the bottom), NOT by the\n")
cat("canonical low -> high median distance used in the other two panels.\n")
cat("Bottom -> top:", paste(delta_levels, collapse = ", "), "\n\n")
hdr_method()
cat("=== Within-celltype observed vs null (primary: mean of per-donor means) ===\n")
print(as.data.frame(tab1), row.names = FALSE)
cat("\n=== Effect sizes ===\n")
print(as.data.frame(tab2), row.names = FALSE)
cat("\n=== Context: what the canonical panel tested ===\n")
print(as.data.frame(tab3), row.names = FALSE)
cat("\nFull detail, per-donor differences and the secondary pooled statistic:\n")
cat("  stats_", slug, ".txt\n", sep = "")
cat("\n------------------------------------------------------------\n")
print(sessionInfo())
sink()

message("Done: ", out_dir)
print(as.data.frame(tab1), row.names = FALSE)
cat("\n")
print(as.data.frame(tab2), row.names = FALSE)
