#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_phf1_state_markers_forest_panel.R
#
# Figure panels: 2E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_phf1_state_markers_forest_panel.R
#
# Fig 2E (IMC state markers, tangle-bearing vs tangle-free excitatory neurons, subcluster-adjusted)
# redrawn as a forest in the style of the Fig 3E forest (R/imc_phf1_distance_markers_subsets.R):
# one point per marker = the per-marker Set 3 linear mixed model estimate,
#
#   value ~ grp + subcluster + Sex + Age_s + PMI_s + (1 | patient_id)
#
# The estimates are NOT refitted here. They are read verbatim from the output
# (plots/imc_phf1_state_markers_subsets/source_data_exc_subcluster_adj.tsv, fit SUB of
# R/imc_phf1_state_markers_subsets.R), so the points are exactly the 2E numbers (log2FC).
#
# Colour follows the 2E call (padj < 0.05 AND |log2FC| >= 0.25), not 3E's padj-only rule: every
# marker here clears padj, so a padj-only rule would colour all ten red. There is no legend on the
# panel. AT8 (positive control) shares the axis.
#
# BAR = donor spread, NOT a confidence interval: the central 95% (2.5th-97.5th percentile) of the
# 44 per-donor contrasts, on the same log2FC mapping, drawn 3E-style. It is a data display -- the
# random-intercept model assumes one effect for every donor, so it has no donor-effect variance to
# draw. The model's own 95% CI (~ +/-0.02-0.06 log2FC, narrower than the point) is not drawn; it
# is in the source data and stats log. Point and colour are unchanged Set 3.
#
# Donor-level check: the log reports cell-level and donor-level effect sizes side by side.
# Each donor's contrast is a within-donor OLS value ~ grp + subcluster (the Set 3 donor-constant
# terms drop out inside one donor), which also feeds the bar; the log reports dz and the number
# of the 44 donors concordant with the model estimate.
# The cell frame is rebuilt with the 2E filter and asserted against the 2E table's counts.
#
# Outputs the standard TRIPLE under plots/imc_phf1_state_markers_forest_panel/. Sized 8.7 x 6.5 cm.

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2)
})

hpc <- "<PROJECT_ROOT>"
loc <- "<PROJECT_ROOT>"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")

out_dir <- "plots/imc_phf1_state_markers_forest_panel"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
SLUG <- "exc_subcluster_adj_forest"

src_file <- "plots/imc_phf1_state_markers_subsets/source_data_exc_subcluster_adj.tsv"
FDR     <- 0.05     # 2E call, copied from R/imc_phf1_state_markers_subsets.R
LFC_MIN <- 0.25

localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot("PHF1_Otsu" %in% colnames(colData(spe)))

# -------------------------------------------------------------------
# 2E model estimates (the plotted numbers)
# -------------------------------------------------------------------
pool <- read.delim(src_file, stringsAsFactors = FALSE)
stopifnot(all(pool$fit == "SUB: all Exc + subcluster"), !anyNA(pool$logFC.L), !anyNA(pool$logFC.R))
markers <- pool$gene
pool <- pool %>% mutate(called = padj < FDR & abs(logFC) >= LFC_MIN)

# -------------------------------------------------------------------
# Cells -- identical filter to R/imc_phf1_state_markers_subsets.R (donor-level check only)
# -------------------------------------------------------------------
cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)

min_rois_per_patient <- 3
cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < min_rois_per_patient) %>% pull(patient_id)

md <- cd %>%
  filter(!celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
         !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1",
         grepl("^Excitatory", celltype_clusters),
         !is.na(PHF1_Otsu)) %>%
  mutate(patient_id = as.character(patient_id),
         celltype_clusters = as.character(celltype_clusters),
         grp = factor(ifelse(PHF1_Otsu == "PHF1_pos", "PHF1+", "PHF1-"),
                      levels = c("PHF1-", "PHF1+")))
# the rebuilt frame must be the one the 2E estimates came from
stopifnot(nrow(md) == pool$n_cells[1], sum(md$grp == "PHF1+") == pool$n_pos[1],
          n_distinct(md$patient_id) == pool$n_donors[1],
          setequal(markers, setdiff(rownames(spe)[rowData(spe)$marker_class == "state"], "PHF1")))

EM <- assay(spe, "exprs")[markers, md$cell_id, drop = FALSE]

donor_info <- md %>% group_by(patient_id) %>%
  summarise(n_pos = sum(grp == "PHF1+"), n_neg = sum(grp == "PHF1-"), .groups = "drop")
stopifnot(all(donor_info$n_pos >= 1), all(donor_info$n_neg >= 1))

donor_est <- bind_rows(lapply(donor_info$patient_id, function(p) {
  i  <- which(md$patient_id == p)
  d0 <- md[i, c("grp", "celltype_clusters")]
  d0$subcluster <- droplevels(factor(d0$celltype_clusters))
  f  <- if (nlevels(d0$subcluster) > 1) value ~ grp + subcluster else value ~ grp
  bind_rows(lapply(markers, function(mk) {
    d0$value <- as.numeric(EM[mk, i])
    co <- coef(summary(lm(f, data = d0)))
    tibble(patient_id = p, gene = mk, est = co["grpPHF1+", "Estimate"])
  }))
})) %>% left_join(donor_info, by = "patient_id") %>%
  left_join(pool %>% select(gene, mu_neg), by = "gene") %>%
  # same monotone mapping as the model estimate, so the bar and the point share one axis
  mutate(logFC = log2(sinh(mu_neg + est) / sinh(mu_neg)))
if (any(!is.finite(donor_est$logFC)))
  stop("a donor contrast falls below -mu_neg; the sinh back-transform is undefined for it")

donor_sum <- donor_est %>%
  left_join(pool %>% select(gene, est_model = est), by = "gene") %>%
  group_by(gene) %>%
  summarise(n_donors = n(), donor_mean = mean(est), donor_sd = sd(est),
            dz = donor_mean / donor_sd,
            n_concordant = sum(sign(est) == sign(est_model[1])),
            donor_q025 = quantile(logFC, 0.025), donor_q975 = quantile(logFC, 0.975),
            .groups = "drop") %>%
  left_join(pool %>% select(gene, est_model = est, cohens_d), by = "gene")

# -------------------------------------------------------------------
# Figure -- construction copied from make_forest() in R/imc_phf1_distance_markers_subsets.R
# -------------------------------------------------------------------
fpd <- pool %>%
  left_join(donor_sum %>% select(gene, donor_q025, donor_q975), by = "gene") %>%
  transmute(gene, x = logFC, lo = logFC.L, hi = logFC.R, called, donor_q025, donor_q975) %>%
  mutate(gene_lab = factor(gene, levels = rev(gene[order(x)])),
         sig = factor(ifelse(called, "called", "not called"), levels = c("not called", "called")))

fig <- ggplot(fpd, aes(x = x, y = gene_lab, colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
  # bar = central 95% of the per-donor contrasts (not the model CI, which is under the point)
  geom_errorbar(aes(xmin = donor_q025, xmax = donor_q975), orientation = "y", width = 0,
                linewidth = 0.4) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("not called" = "grey65", "called" = "#BD0026"),
                      drop = FALSE, guide = "none") +
  labs(x = expression(Log[2]*" (fold-change), tangle-bearing / tangle-free"), y = NULL) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none", plot.margin = margin(4, 8, 4, 4, unit = "pt"))

ggsave(file.path(out_dir, sprintf("plot_%s.pdf", SLUG)), fig,
       width = 8.7, height = 6.5, units = "cm", device = "pdf")

# -------------------------------------------------------------------
# Source data + stats
# -------------------------------------------------------------------
wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)
wr(pool %>% left_join(donor_sum %>% select(gene, donor_q025, donor_q975), by = "gene") %>%
     select(gene, est, SE, df, t, pval, padj, logFC, logFC.L, logFC.R, cohens_d, called,
            donor_q025, donor_q975, mu_neg, n_cells, n_pos, n_neg, n_donors, re_var_donor,
            re_pct_donor, singular),
   sprintf("source_data_%s.tsv", SLUG))
wr(pool %>% select(gene, est, logFC, logFC.L, logFC.R, cohens_d, pval, padj) %>%
     left_join(donor_sum %>% select(gene, n_donors, dz, n_concordant), by = "gene"),
   sprintf("stats_%s_effectsize.tsv", SLUG))
wr(donor_est %>% select(gene, patient_id, est, logFC, n_pos, n_neg),
   sprintf("stats_%s_donor_contrasts.tsv", SLUG))

sink(file.path(out_dir, sprintf("stats_%s.txt", SLUG)))
title <- "Tangle-bearing vs tangle-free state markers, subcluster-adjusted (IMC, Fig 2E forest)"
cat(title, "\n"); cat(strrep("=", nchar(title)), "\n", sep = "")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("MODEL (Set 3, per marker; estimates read verbatim from ", src_file, "):\n", sep = "")
cat("  value ~ grp + subcluster + Sex + Age_s + PMI_s + (1 | patient_id)\n")
cat("  Outcome asinh intensity, REML (lmerTest). grp = PHF1+ (tangle-bearing) vs PHF1-\n")
cat("  (tangle-free); contrast read directly as grpPHF1+ (PHF1+ minus PHF1-).\n")
cat("  log2FC = log2(sinh(mu_neg + est)/sinh(mu_neg)), mu_neg = observed mean asinh in PHF1-\n")
cat("  cells; sinh is monotone, so the 95% CI maps straight onto the log2FC bounds.\n")
cat(sprintf("  BH across %d markers. Red = padj < %.2f AND |log2FC| >= %.2f (2E call); grey = below\n",
            length(markers), FDR, LFC_MIN))
cat("  that call. Every marker here has padj < 0.05, so grey means below the effect-size floor.\n")
cat("POINT = model estimate. BAR = central 95% (2.5th-97.5th percentile, quantile type 7) of the\n")
cat("  44 within-donor contrasts on the same log2FC mapping -- donor spread, NOT a CI. The model\n")
cat("  95% CI (logFC.L / logFC.R below) is narrower than the point and is not drawn.\n\n")
cat("Cells:", pool$n_cells[1], " PHF1+:", pool$n_pos[1], " PHF1-:", pool$n_neg[1],
    " Donors:", pool$n_donors[1], "(rebuilt frame verified equal)\n\n")

cat("=== Results (points in the figure; logFC.L/R = model 95% CI, not drawn) ===\n")
print(as.data.frame(pool %>% arrange(desc(logFC)) %>%
  select(gene, est, SE, df, t, padj, logFC, logFC.L, logFC.R, cohens_d, called)),
  row.names = FALSE, digits = 3)
cat("\nCalled:", sum(pool$called), "of", nrow(pool), "\n")

cat("\n=== Random effects (donor), model per marker ===\n")
cat("A variance at ~0 with singular TRUE means the donor random-intercept variance is\n")
cat("estimated at zero.\n")
print(as.data.frame(pool %>% select(gene, re_var_donor, re_pct_donor, singular)),
      row.names = FALSE, digits = 4)

cat("\n=== Effect sizes: cell level vs donor level ===\n")
cat("cohens_d = model estimate / sqrt(donor RE variance + residual) -- cell-level.\n")
cat("Donor level: within-donor OLS value ~ grp + subcluster in each donor (subcluster dropped\n")
cat("  where only one is present); dz = mean / SD of the 44 donor contrasts; n_concordant =\n")
cat("  donors whose contrast has the sign of the model estimate. donor_q025/q975 = the bar.\n\n")
print(as.data.frame(donor_sum %>% arrange(desc(est_model)) %>%
  select(gene, est_model, cohens_d, donor_mean, dz, n_concordant, n_donors, donor_q025,
         donor_q975)),
  row.names = FALSE, digits = 3)
cat("\nPHF1+ cells per donor: median", median(donor_info$n_pos), " range",
    paste(range(donor_info$n_pos), collapse = "-"), "\n")

cat("\nNOTES\n-----\n")
cat("- AT8 is a POSITIVE CONTROL: it stains the same tau species as PHF1.\n")
cat("- subcluster is listed as a sensitivity term in docs/MODELS.md Set 3; it is in this model because\n")
cat("  the 2E panel uses the subcluster-adjusted fit.\n")
cat("- PHF1_Otsu is a per-donor Otsu threshold with a hard floor, so tangle-bearing is a\n")
cat("  thresholded call on a continuous stain.\n")
cat("- sample_id (ROI) is not modelled (Set 3).\n")
cat("\n"); print(sessionInfo())
sink()

cat("\n== Forest (2E model) ==\n")
print(as.data.frame(pool %>% arrange(desc(logFC)) %>%
  select(gene, logFC, logFC.L, logFC.R, padj, called)), row.names = FALSE, digits = 3)
cat("\n== Donor concordance ==\n")
print(as.data.frame(donor_sum %>% select(gene, dz, n_concordant, n_donors)),
      row.names = FALSE, digits = 3)
