#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_phf1_distance_markers_subsets.R
#
# Figure panels: 3E
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_phf1_distance_markers_subsets.R
#
# Two SUBSET/COVARIATE variants of the all-excitatory IMC distance-marker model, both fitted with
# the Set 3 model, which adjusts for NEITHER amyloid proximity NOR cortical depth:
#
#   A. RELN+CALB1 only   value ~ dist_z + Sex + Age_s + PMI_s + (1 | patient_id)
#                        restricted to "Excitatory neuron cluster 4 (RELN)" +
#                        "Excitatory neuron cluster 2 (CALB1)" pooled. These are the classic
#                        tangle-vulnerable entorhinal populations, so the question is whether the
#                        gradient is carried by them specifically.
#
#   B. Subcluster-adjusted  value ~ dist_z + celltype_clusters + Sex + Age_s + PMI_s + (1|patient_id)
#                        all four excitatory subclusters, with subcluster as a fixed effect. This
#                        asks whether the gradient survives WITHIN subcluster, i.e. whether it is a
#                        real per-cell shift or just a change in which excitatory subtype sits near
#                        tangles. Same cells as the all-excitatory model, so its dist_z
#                        coefficient is directly comparable to it.
#
# Everything else -- cohort filter, anchors, distance construction, cap, marker set, sign
# convention, arcsinh->log2FC back-transform, figure style -- is shared with the other IMC
# distance panels.
#
# NOTE on dist_z: the canonical transform is log(dist)/sd(log(dist)) with the SD computed
# WITHIN each celltype. Analysis A is a different cell set from the all-excitatory model, so its
# dist_sd is recomputed and differs; its coefficients are therefore in ITS OWN s.d. units and are
# not directly comparable in magnitude. Analysis B uses the all-excitatory cell set, so its dist_sd
# is identical and its coefficients ARE directly comparable. Both SDs are printed.
#
# Outputs the standard TRIPLE under plots/imc_phf1_distance_markers_subsets/, per analysis:
#   plot_<slug>_forest.pdf / plot_<slug>_volcano.pdf / source_data_<slug>.tsv / stats_<slug>.txt

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(ggrepel); library(lmerTest); library(RANN)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/imc_utils.R")   # re_variance_summary(), re_total_sd(): donor RE variance in every log

out_dir <- "plots/imc_phf1_distance_markers_subsets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR          <- 0.05
MAX_DIST     <- 300
EDGE_BUFFER  <- 0
MIN_ANCHORS  <- 1
VULN_CLUSTERS <- c("Excitatory neuron cluster 4 (RELN)",
                   "Excitatory neuron cluster 2 (CALB1)")

localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot("PHF1_Otsu" %in% colnames(colData(spe)),
          "Matched_4G8_40" %in% colnames(colData(spe)))

# -------------------------------------------------------------------
# Cells -- identical backbone to the other IMC distance panels
# -------------------------------------------------------------------
cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)
sc <- spatialCoords(spe)
cd$x <- sc[, "Pos_X"]; cd$y <- sc[, "Pos_Y"]

min_rois_per_patient <- 3
cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < min_rois_per_patient) %>% pull(patient_id)

base <- cd %>%
  filter(!celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
         !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1",
         !is.na(PHF1_Otsu)) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
         celltype_clusters = as.character(celltype_clusters),
         Sex = factor(Sex))
stopifnot(!anyNA(base$Sex), !anyNA(base$Age), !anyNA(base$PMI))
base <- base %>% mutate(Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)))

base <- base %>% group_by(sample_id) %>%
  mutate(edge_dist = pmin(x - min(x), max(x) - x, y - min(y), max(y) - y)) %>%
  ungroup()

stopifnot(all(VULN_CLUSTERS %in% unique(base$celltype_clusters)))

anchors <- base %>% filter(PHF1_Otsu == "PHF1_pos",
                           grepl("neuron", celltype_clusters, ignore.case = TRUE))

nearest_dist <- function(tg, an) {
  out <- rep(NA_real_, nrow(tg))
  for (s in unique(tg$sample_id)) {
    i <- which(tg$sample_id == s)
    a <- an[an$sample_id == s, c("x", "y"), drop = FALSE]
    if (!nrow(a)) next
    out[i] <- RANN::nn2(as.matrix(a), as.matrix(tg[i, c("x", "y")]), k = 1)$nn.dists[, 1]
  }
  out
}

anchor_n <- anchors %>% count(sample_id, name = "n_anchors")
roi_ok   <- anchor_n$sample_id[anchor_n$n_anchors >= MIN_ANCHORS]

# Build a modelled frame for an arbitrary target-cell filter. dist_z is recomputed within the
# returned set ("dist_sd is computed within each celltype").
build_md <- function(target_filter) {
  tg <- base %>% filter(PHF1_Otsu == "PHF1_neg") %>% target_filter()
  tg$dist_um <- nearest_dist(tg, anchors)
  md <- tg %>%
    filter(!is.na(dist_um), dist_um <= MAX_DIST, edge_dist >= EDGE_BUFFER,
           sample_id %in% roi_ok) %>%
    mutate(plaque = factor(ifelse(Matched_4G8_40, "plaque-prox", "plaque-distal"),
                           levels = c("plaque-distal", "plaque-prox")),
           dist_z = log(dist_um) / sd(log(dist_um)))
  if (any(md$dist_um <= 0)) stop("non-positive distance; log() would fail")
  md
}

md_vuln <- build_md(function(d) d %>% filter(celltype_clusters %in% VULN_CLUSTERS))
md_exc  <- build_md(function(d) d %>% filter(grepl("^Excitatory", celltype_clusters)))
md_exc$subcluster <- factor(md_exc$celltype_clusters)
# reference level = the largest subcluster, so the intercept is the modal excitatory neuron
md_exc$subcluster <- relevel(md_exc$subcluster,
                             ref = names(sort(table(md_exc$subcluster), decreasing = TRUE))[1])

markers <- setdiff(rownames(spe)[rowData(spe)$marker_class == "state"], "PHF1")
EM <- assay(spe, "exprs")

# -------------------------------------------------------------------
# Per-marker model. Sign: per s.d. CLOSER to a tangle (positive = higher near tangles).
# -------------------------------------------------------------------
grab <- function(m, term) {
  if (is.null(m)) return(rep(NA_real_, 5))
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(rep(NA_real_, 5))
  c(co[term, "Estimate"], co[term, "Std. Error"], co[term, "df"],
    co[term, "t value"], co[term, "Pr(>|t|)"])
}

fit_set <- function(md, extra_terms = character(0)) {
  em <- EM[markers, md$cell_id, drop = FALSE]
  rhs <- paste(c("dist_z", extra_terms, "Sex", "Age_s", "PMI_s", "(1 | patient_id)"),
               collapse = " + ")
  f_main <- as.formula(paste("value ~", rhs))
  # the with-amyloid counterpart, reported alongside per the with/without convention
  f_plq  <- as.formula(paste("value ~", paste(c("dist_z", extra_terms, "plaque", "Sex", "Age_s",
                                                "PMI_s", "(1 | patient_id)"), collapse = " + ")))

  bind_rows(lapply(markers, function(mk) {
    d <- md %>% transmute(value = as.numeric(em[mk, ]), dist_z, plaque, patient_id,
                          Sex, Age_s, PMI_s,
                          subcluster = if ("subcluster" %in% names(md)) md$subcluster else NA)
    d <- d[is.finite(d$value), ]
    fit <- function(f) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE), error = function(e) NULL)
    m  <- fit(f_main); mp <- fit(f_plq)
    a  <- grab(m, "dist_z"); ap <- grab(mp, "dist_z")
    crit <- stats::qt(0.975, ifelse(is.na(a[3]), Inf, a[3]))
    sd_tot <- re_total_sd(m)                                   # R/imc_utils.R
    rv <- re_variance_summary(m); rv <- rv[rv$group == "patient_id", , drop = FALSE]
    tibble(gene = mk,
           beta_closer = -a[1], SE = a[2], df = a[3], t = -a[4], pval = a[5],
           CI.L = -a[1] - crit * a[2], CI.R = -a[1] + crit * a[2],
           cohens_d = -a[1] / sd_tot,
           beta_closer_plaque_adj = -ap[1], pval_plaque_adj = ap[5],
           AveExpr = mean(d$value), n_cells = nrow(d),
           n_donors = dplyr::n_distinct(d$patient_id),
           re_var_donor = if (nrow(rv)) rv$vcov[1] else NA_real_,
           re_pct_donor = if (nrow(rv)) rv$pct_of_total[1] else NA_real_,
           singular = if (is.null(m)) NA else lme4::isSingular(m))
  })) %>%
    mutate(padj = p.adjust(pval, method = "BH"),
           padj_plaque_adj = p.adjust(pval_plaque_adj, method = "BH"),
           # exact inverse of asinh, so a ratio on the intensity scale is a genuine fold change
           logFC   = log2(sinh(AveExpr + beta_closer) / sinh(AveExpr)),
           logFC.L = log2(sinh(AveExpr + CI.L) / sinh(AveExpr)),
           logFC.R = log2(sinh(AveExpr + CI.R) / sinh(AveExpr))) %>%
    arrange(pval)
}

res_vuln <- fit_set(md_vuln)
res_exc  <- fit_set(md_exc, extra_terms = "subcluster")

# -------------------------------------------------------------------
# Figures -- style shared with the other IMC distance panels
# -------------------------------------------------------------------
XLAB_LOG2FC <- expression(Log[2]*" FC per s.d. nearer tangle")

make_forest <- function(d) {
  fpd <- d %>%
    transmute(gene, x = logFC, lo = logFC.L, hi = logFC.R, padj_use = padj) %>%
    mutate(gene_lab = factor(gene, levels = rev(gene[order(x)])),
           sig = ifelse(padj_use < FDR, "padj < 0.05", "n.s."))
  ggplot(fpd, aes(x = x, y = gene_lab, colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.4) +
    geom_point(size = 1.7) +
    scale_colour_manual(values = c("padj < 0.05" = "#BD0026", "n.s." = "grey65"), name = NULL) +
    labs(x = XLAB_LOG2FC, y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", plot.margin = margin(4, 8, 4, 4, unit = "pt"))
}

make_volcano <- function(dt) {
  dt <- as.data.frame(dt)
  dt <- dt[!is.na(dt$padj) & !is.nan(dt$logFC), ]
  dt$de <- "Not sig"
  dt$de[dt$padj <= FDR & dt$logFC >  0] <- "Up"
  dt$de[dt$padj <= FDR & dt$logFC <= 0] <- "Down"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  dt$lab <- ifelse(dt$de != "Not sig", dt$gene, "")
  if (any(dt$padj == 0)) dt$padj[dt$padj == 0] <- min(dt$padj[dt$padj > 0])
  max_x <- max(abs(dt$logFC), na.rm = TRUE) * 1.15
  max_y <- max(-log10(dt$padj), na.rm = TRUE) * 1.12
  ggplot(dt, aes(x = logFC, y = -log10(padj))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, alpha = 0.5) +
    geom_hline(yintercept = -log10(FDR), linetype = 2, linewidth = 0.3, alpha = 0.5) +
    geom_point(aes(colour = de), size = 1.4, alpha = 0.9, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      aes(label = lab, colour = de), size = 2, max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey45",
      box.padding = 0.45, point.padding = 0.15, force = 4, force_pull = 0.5,
      max.iter = 20000, max.time = 2, seed = 42, na.rm = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = c("Up" = "#DC0000FF", "Down" = "#3C5488FF",
                                   "Not sig" = "grey70"), guide = "none") +
    scale_x_continuous(limits = c(-max_x, max_x), oob = scales::squish) +
    scale_y_continuous(limits = c(0, max_y), oob = scales::squish) +
    labs(x = XLAB_LOG2FC, y = expression("-" * log[10] * " (adjusted p-value)")) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8, unit = "pt"))
}

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

emit <- function(res, md, slug, title, model_txt, extra_notes) {
  ggsave(file.path(out_dir, sprintf("plot_%s_forest.pdf", slug)), make_forest(res),
         width = 8.5, height = 6.4, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, sprintf("plot_%s_volcano.pdf", slug)), make_volcano(res),
         width = 6.3, height = 5.4, units = "cm", device = "pdf")
  wr(res, sprintf("source_data_%s.tsv", slug))
  wr(res %>% select(gene, beta_closer, CI.L, CI.R, cohens_d, logFC, logFC.L, logFC.R, pval, padj),
     sprintf("stats_%s_effectsize.tsv", slug))

  sink(file.path(out_dir, sprintf("stats_%s.txt", slug)))
  cat(title, "\n"); cat(strrep("=", nchar(title)), "\n", sep = "")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("MODEL (Set 3; adjusts for NEITHER amyloid proximity NOR cortical depth):\n")
  cat("  ", model_txt, "\n\n", sep = "")
  cat("  A with-amyloid counterpart is fitted for every marker and reported in the source data as\n")
  cat("  beta_closer_plaque_adj / padj_plaque_adj, matching the with/without convention used by\n")
  cat("  the other IMC distance panels. The FIGURES show the no-amyloid, no-depth model.\n\n")
  cat("Cells:", nrow(md), " Donors:", dplyr::n_distinct(md$patient_id),
      " ROIs:", dplyr::n_distinct(md$sample_id), "\n")
  cat("Subclusters in the modelled set:\n")
  print(sort(table(md$celltype_clusters), decreasing = TRUE))
  cat("\ndist_z = log(dist_um)/sd(log(dist_um)), plain log, /SD, NOT centred (canonical).\n")
  cat(sprintf("  dist_sd for THIS cell set = %.5f, computed within the set.\n",
              sd(log(md$dist_um))))
  cat(sprintf("  Distance cap %d um; distance to nearest PHF1+ NEURON, computed within ROI.\n",
              MAX_DIST))
  cat("Sign: coefficients are per s.d. CLOSER to a tangle, so POSITIVE = higher near tangles.\n")
  cat("logFC = log2( sinh(AveExpr + beta) / sinh(AveExpr) ); the model is fitted on asinh\n")
  cat("  intensity so sinh() is the exact inverse and the ratio is a genuine fold change.\n\n")
  cat(extra_notes, sep = "\n"); cat("\n\n")
  cat("=== Results (BH across the", length(markers), "state markers) ===\n")
  print(as.data.frame(res %>%
    select(gene, beta_closer, CI.L, CI.R, logFC, cohens_d, padj,
           beta_closer_plaque_adj, padj_plaque_adj, n_cells)), row.names = FALSE, digits = 3)
  cat("\nSignificant at padj <", FDR, ":", sum(res$padj < FDR), "of", nrow(res), "\n")
  cat("\n=== Random effects (donor), headline model per marker ===\n")
  cat("A variance at ~0 with singular TRUE means the model does NOT adjust for donor.\n")
  print(as.data.frame(res %>% select(gene, re_var_donor, re_pct_donor, singular)),
        row.names = FALSE, digits = 4)
  cat("\nNOTES\n-----\n")
  cat("- AT8 stains the same tau species as PHF1, so its gradient serves as a positive control\n")
  cat("  for the distance construct.\n")
  cat("- Effect sizes: cohens_d and the 95% CI are reported for every marker.\n")
  cat("- No cortical-depth term (Set 3; depth enters as a sensitivity term in\n")
  cat("  R/imc_covariate_sensitivity.R).\n")
  cat("\n"); print(sessionInfo())
  sink()
}

emit(res_vuln, md_vuln, "reln_calb1",
     "Distance gradient in RELN+ and CALB1+ excitatory neurons (IMC, AD cases)",
     "value ~ dist_z + Sex + Age_s + PMI_s + (1 | patient_id)",
     c("CELL SET: the two tangle-vulnerable entorhinal excitatory populations pooled --",
       paste0("  ", paste(VULN_CLUSTERS, collapse = " + ")), "",
       "Both clusters are named for markers that are TYPE markers, i.e. they were inputs to the",
       "PCA -> Harmony -> Phenograph clustering that defined them. That circularity affects the",
       "cluster DEFINITION, not this analysis: the outcomes modelled here are STATE markers, which",
       "were never clustering inputs, so the gradient is not being read off the variables that",
       "created the subset.", "",
       "dist_sd is recomputed within this subset, so these coefficients are in this subset's own",
       "s.d. units and are NOT directly comparable in magnitude to the all-excitatory panel."))

emit(res_exc, md_exc, "exc_subcluster_adj",
     "Distance gradient in excitatory neurons, subcluster-adjusted (IMC, AD cases)",
     "value ~ dist_z + subcluster + Sex + Age_s + PMI_s + (1 | patient_id)",
     c("CELL SET: all four excitatory subclusters (the all-excitatory cell set), so dist_sd",
       "matches the all-excitatory model and the dist_z coefficient IS directly comparable to it.", "",
       "WHAT THE SUBCLUSTER TERM BUYS: without it, a marker can appear to rise near tangles purely",
       "because the MIX of excitatory subtypes shifts with distance. With subcluster as a fixed",
       "effect the coefficient is the WITHIN-subcluster gradient.",
       paste0("  Reference subcluster: ", levels(md_exc$subcluster)[1])))

cat("\n== A. RELN + CALB1 ==\n")
print(as.data.frame(res_vuln %>% select(gene, beta_closer, logFC, padj, n_cells)),
      row.names = FALSE, digits = 3)
cat("cells:", nrow(md_vuln), " donors:", dplyr::n_distinct(md_vuln$patient_id),
    " sig:", sum(res_vuln$padj < FDR), "\n")
cat("\n== B. All excitatory, subcluster-adjusted ==\n")
print(as.data.frame(res_exc %>% select(gene, beta_closer, logFC, padj, n_cells)),
      row.names = FALSE, digits = 3)
cat("cells:", nrow(md_exc), " donors:", dplyr::n_distinct(md_exc$patient_id),
    " sig:", sum(res_exc$padj < FDR), "\n")
