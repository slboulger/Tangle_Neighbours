#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# imc_phf1_state_markers_subsets.R
#
# Figure panels: 2E - upstream step; also S8B
# Writes the per-marker estimates (plots/imc_phf1_state_markers_subsets/) that
# imc_phf1_state_markers_forest_panel.R draws and imc_covariate_sensitivity.R reads.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# imc_phf1_state_markers_subsets.R
#
# PHF1+ vs PHF1- state markers in excitatory neurons, with subset / covariate variants, mirroring
# what R/imc_phf1_distance_markers_subsets.R does for the distance analysis. Three fits, all on
# the same 11 state markers minus PHF1 itself:
#
#   REF   all 4 Exc subclusters   value ~ grp + Sex + Age_s + PMI_s + (1 | patient_id)
#   SUB   all 4 Exc subclusters   value ~ grp + subcluster + Sex + Age_s + PMI_s + (1|patient_id)
#   VULN  RELN + CALB1 pooled     value ~ grp + Sex + Age_s + PMI_s + (1 | patient_id)
#
# REF is the like-for-like baseline for SUB: the two differ only by the subcluster term, so
# REF vs SUB isolates it.
#
# Sex/Age/PMI match the CosMx counterpart deg_pb_category_phf1.r, whose default confounder set is
# Age,Sex,PostMortemInterval, and the Set 3 specification (docs/MODELS.md).
#
# WHY THE SUBCLUSTER TERM MATTERS: PHF1+ and PHF1- cells are not drawn evenly from the four
# excitatory subclusters, so a marker can separate the two groups simply because tangle-bearing
# cells are disproportionately one subtype. With subcluster as a fixed effect the contrast is the
# WITHIN-subcluster difference.
#
# Conventions: arcsinh outcome, grpPHF1+ read directly off the coefficient table (never
# emmeans::pairs, which returns PHF1- minus PHF1+), sinh back-transform to a genuine fold change,
# BH then a PADJ_FLOOR cap of 1e-300, LFC_MIN 0.25 shading, 6.3 x 5.4 cm volcano.
#
# Outputs the standard TRIPLE under plots/imc_phf1_state_markers_subsets/, per fit.

suppressPackageStartupMessages({
  library(SpatialExperiment); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(ggrepel); library(lmerTest)
})

hpc <- "<PROJECT_ROOT>/phf1_v2"
loc <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(hpc)) hpc else loc)
source("R/palettes.R")
source("R/imc_utils.R")   # re_variance_summary(), re_total_sd(): donor RE variance in every log

out_dir <- "plots/imc_phf1_state_markers_subsets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR        <- 0.05
LFC_MIN    <- 0.25
PADJ_FLOOR <- 1e-300
VULN_CLUSTERS <- c("Excitatory neuron cluster 4 (RELN)",
                   "Excitatory neuron cluster 2 (CALB1)")

localwd <- "<IMC_ROOT>/"
spe <- readRDS(paste0(localwd, "spe.rds"))
stopifnot("PHF1_Otsu" %in% colnames(colData(spe)),
          "Matched_4G8_40" %in% colnames(colData(spe)))

# -------------------------------------------------------------------
# Cells -- same cohort frame as the other IMC scripts
# -------------------------------------------------------------------
cd <- as.data.frame(colData(spe))
cd$cell_id <- colnames(spe)

min_rois_per_patient <- 3
cv_removal <- cd %>% distinct(patient_id, sample_id) %>% count(patient_id) %>%
  filter(n < min_rois_per_patient) %>% pull(patient_id)

md_all <- cd %>%
  filter(!celltype_clusters %in% c("Artefact cluster", "Unassigned cluster"),
         !patient_id %in% cv_removal,
         BraakGroup != "Braak_0_1",
         grepl("^Excitatory", celltype_clusters),
         !is.na(PHF1_Otsu)) %>%
  mutate(patient_id = as.character(patient_id), sample_id = as.character(sample_id),
         celltype_clusters = as.character(celltype_clusters),
         Sex = factor(Sex),
         grp = factor(ifelse(PHF1_Otsu == "PHF1_pos", "PHF1+", "PHF1-"),
                      levels = c("PHF1-", "PHF1+")),
         plaque = factor(ifelse(Matched_4G8_40, "plaque-prox", "plaque-distal"),
                         levels = c("plaque-distal", "plaque-prox")))
stopifnot(!anyNA(md_all$Sex), !anyNA(md_all$Age), !anyNA(md_all$PMI))
md_all <- md_all %>% mutate(Age_s = as.numeric(scale(Age)), PMI_s = as.numeric(scale(PMI)))

md_all$subcluster <- factor(md_all$celltype_clusters)
md_all$subcluster <- relevel(md_all$subcluster,
                             ref = names(sort(table(md_all$subcluster), decreasing = TRUE))[1])

stopifnot(all(VULN_CLUSTERS %in% levels(md_all$subcluster)))
md_vuln <- md_all %>% filter(celltype_clusters %in% VULN_CLUSTERS)

markers <- setdiff(rownames(spe)[rowData(spe)$marker_class == "state"], "PHF1")
EM <- assay(spe, "exprs")

# -------------------------------------------------------------------
# Per-marker fit
# -------------------------------------------------------------------
# Reference intensity for the back-transform: the OBSERVED mean arcsinh in the PHF1- cells, not
# co["(Intercept)"]. With grp as the only fixed effect the intercept would BE the PHF1- fitted
# mean; here the intercept is additionally conditioned on Sex = reference and subcluster =
# reference, which would shift the denominator of the fold change for reasons unrelated to PHF1. The functional form log2(sinh(mu + est)/sinh(mu)) is unchanged.
fit_set <- function(md, extra_terms = character(0), label = "") {
  em <- EM[markers, md$cell_id, drop = FALSE]
  rhs   <- paste(c("grp", extra_terms, "Sex", "Age_s", "PMI_s", "(1 | patient_id)"),
                 collapse = " + ")
  f_main <- as.formula(paste("value ~", rhs))
  f_plq  <- as.formula(paste("value ~", paste(c("grp", extra_terms, "plaque", "Sex", "Age_s",
                                                "PMI_s", "(1 | patient_id)"), collapse = " + ")))

  bind_rows(lapply(markers, function(mk) {
    d <- md %>% transmute(value = as.numeric(em[mk, ]), grp, plaque, patient_id,
                          Sex, Age_s, PMI_s, subcluster)
    d <- d[is.finite(d$value), ]
    fit <- function(f) tryCatch(lmerTest::lmer(f, data = d, REML = TRUE), error = function(e) NULL)
    m <- fit(f_main); mp <- fit(f_plq)
    if (is.null(m)) return(tibble(gene = mk, est = NA_real_))

    co <- summary(m)$coefficients
    # grpPHF1+ is already oriented PHF1+ minus PHF1-; read it directly, never via emmeans::pairs()
    stopifnot("grpPHF1+" %in% rownames(co))
    est <- co["grpPHF1+", "Estimate"]; se <- co["grpPHF1+", "Std. Error"]
    df  <- co["grpPHF1+", "df"];       tv <- co["grpPHF1+", "t value"]
    pv  <- co["grpPHF1+", "Pr(>|t|)"]
    crit <- stats::qt(0.975, df)

    mu_neg <- mean(d$value[d$grp == "PHF1-"])
    stopifnot(sinh(mu_neg) > 0)
    fc_of <- function(delta) sinh(mu_neg + delta) / sinh(mu_neg)

    sd_tot <- re_total_sd(m)                                   # R/imc_utils.R
    rv <- re_variance_summary(m); rv <- rv[rv$group == "patient_id", , drop = FALSE]
    ep <- if (is.null(mp) || !"grpPHF1+" %in% rownames(summary(mp)$coefficients)) {
      c(NA_real_, NA_real_)
    } else summary(mp)$coefficients["grpPHF1+", c("Estimate", "Pr(>|t|)")]

    tibble(gene = mk, est = est, SE = se, df = df, t = tv, pval = pv,
           logFC = log2(fc_of(est)),
           logFC.L = log2(fc_of(est - crit * se)),
           logFC.R = log2(fc_of(est + crit * se)),
           cohens_d = est / sd_tot,
           est_plaque_adj = ep[1], pval_plaque_adj = ep[2],
           mu_neg = mu_neg, n_cells = nrow(d),
           n_pos = sum(d$grp == "PHF1+"), n_neg = sum(d$grp == "PHF1-"),
           n_donors = dplyr::n_distinct(d$patient_id),
           re_var_donor = if (nrow(rv)) rv$vcov[1] else NA_real_,
           re_pct_donor = if (nrow(rv)) rv$pct_of_total[1] else NA_real_,
           singular = lme4::isSingular(m))
  })) %>%
    mutate(padj = pmax(p.adjust(pval, method = "BH"), PADJ_FLOOR),
           padj_plaque_adj = pmax(p.adjust(pval_plaque_adj, method = "BH"), PADJ_FLOOR),
           p_underflow = pval == 0,
           fit = label) %>%
    arrange(pval)
}

res_ref  <- fit_set(md_all,  character(0),   "REF: all Exc")
res_sub  <- fit_set(md_all,  "subcluster",   "SUB: all Exc + subcluster")
res_vuln <- fit_set(md_vuln, character(0),   "VULN: RELN + CALB1")

# -------------------------------------------------------------------
# Volcano
# -------------------------------------------------------------------
volcano_plot <- function(dt, lfc_threshold = LFC_MIN, fdr_threshold = FDR) {
  dt <- as.data.frame(dt); dt <- dt[!is.nan(dt$logFC), ]
  dt <- dt[order(dt$padj), ]
  dt$de <- "Not sig"
  dt$de[dt$padj <= fdr_threshold & dt$logFC >=  lfc_threshold] <- "Up"
  dt$de[dt$padj <= fdr_threshold & dt$logFC <= -lfc_threshold] <- "Down"
  dt$de <- factor(dt$de, levels = c("Up", "Down", "Not sig"))
  dt$lab <- ifelse(dt$de != "Not sig", dt$gene, "")
  if (any(dt$padj == 0)) dt$padj[dt$padj == 0] <- min(dt$padj[dt$padj > 0])
  max_x <- max(abs(dt$logFC), na.rm = TRUE) * 1.15
  max_y <- max(-log10(dt$padj), na.rm = TRUE) * 1.12
  ggplot(dt, aes(x = logFC, y = -log10(padj))) +
    annotate("rect", xmin = -lfc_threshold, xmax = lfc_threshold,
             ymin = -Inf, ymax = Inf, fill = "grey55", alpha = 0.12) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = 2, linewidth = 0.3, alpha = 0.5) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = 2, linewidth = 0.3, alpha = 0.5) +
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
    # ASCII hyphen, not a Unicode minus (default pdf device)
    labs(x = expression(Log[2]*" (fold-change), PHF1+ / PHF1-"),
         y = expression("-" * log[10] * " (adjusted p-value)")) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8, unit = "pt"))
}

make_forest <- function(d) {
  fpd <- d %>%
    transmute(gene, x = logFC, lo = logFC.L, hi = logFC.R, padj_use = padj) %>%
    mutate(gene_lab = ifelse(gene == "AT8", "AT8 (pos. control)", gene)) %>%
    mutate(gene_lab = factor(gene_lab, levels = rev(gene_lab[order(x)])),
           sig = ifelse(padj_use < FDR & abs(x) >= LFC_MIN, "padj < 0.05", "n.s."))
  ggplot(fpd, aes(x = x, y = gene_lab, colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.4) +
    geom_point(size = 1.7) +
    scale_colour_manual(values = c("padj < 0.05" = "#BD0026", "n.s." = "grey65"), name = NULL) +
    labs(x = expression(Log[2]*" (fold-change), PHF1+ / PHF1-"), y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top", plot.margin = margin(4, 8, 4, 4, unit = "pt"))
}

wr <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE)

emit <- function(res, md, slug, title, model_txt, notes) {
  ggsave(file.path(out_dir, sprintf("plot_%s_volcano.pdf", slug)), volcano_plot(res),
         width = 6.3, height = 5.4, units = "cm", device = "pdf")
  ggsave(file.path(out_dir, sprintf("plot_%s_forest.pdf", slug)), make_forest(res),
         width = 8.5, height = 6.4, units = "cm", device = "pdf")
  wr(res, sprintf("source_data_%s.tsv", slug))
  wr(res %>% select(gene, est, logFC, logFC.L, logFC.R, cohens_d, pval, padj),
     sprintf("stats_%s_effectsize.tsv", slug))

  sink(file.path(out_dir, sprintf("stats_%s.txt", slug)))
  cat(title, "\n"); cat(strrep("=", nchar(title)), "\n", sep = "")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("MODEL: ", model_txt, "\n", sep = "")
  cat("  Outcome = assay(spe,'exprs') = asinh(counts/1). REML.\n")
  cat("  Contrast read directly off the coefficient table as grpPHF1+, already oriented\n")
  cat("  PHF1+ minus PHF1-. Deliberately NOT emmeans::pairs(), which returns PHF1- minus PHF1+.\n")
  cat("  A with-amyloid counterpart is fitted per marker and reported as est_plaque_adj /\n")
  cat("  padj_plaque_adj, matching the with/without convention of the other IMC panels.\n\n")
  cat("Cells:", nrow(md), " PHF1+:", sum(md$grp == "PHF1+"), " PHF1-:", sum(md$grp == "PHF1-"),
      " Donors:", dplyr::n_distinct(md$patient_id), "\n")
  cat("Subclusters in the modelled set (PHF1+ / PHF1- split):\n")
  print(table(md$celltype_clusters, md$grp))
  cat("\nlogFC = log2( sinh(mu_neg + est) / sinh(mu_neg) ), mu_neg = observed mean arcsinh in the\n")
  cat("  PHF1- cells. sinh() is the exact inverse of the asinh outcome, so the ratio is a genuine\n")
  cat("  fold change; sinh is monotone so the contrast CI maps straight onto the logFC bounds.\n")
  cat("  Note: sinh(mean of arcsinh) is a pseudo-median, not an arithmetic mean (Jensen).\n")
  cat(sprintf("BH across %d markers, then padj capped at %s (pval left uncapped).\n",
              length(markers), format(PADJ_FLOOR, scientific = TRUE)))
  cat(sprintf("Effect-size floor for calling a marker changed: |log2FC| >= %.2f.\n\n", LFC_MIN))
  cat(notes, sep = "\n"); cat("\n\n")
  cat("=== Results ===\n")
  print(as.data.frame(res %>%
    select(gene, est, logFC, logFC.L, logFC.R, cohens_d, padj,
           est_plaque_adj, padj_plaque_adj, n_cells)), row.names = FALSE, digits = 3)
  cat("\nSignificant at padj <", FDR, "AND |log2FC| >=", LFC_MIN, ":",
      sum(res$padj < FDR & abs(res$logFC) >= LFC_MIN), "of", nrow(res), "\n")
  cat("\n=== Random effects (donor), model per marker ===\n")
  cat("A variance at ~0 with singular TRUE means the donor random-intercept variance is\n")
  cat("estimated at zero.\n")
  print(as.data.frame(res %>% select(gene, re_var_donor, re_pct_donor, singular)),
        row.names = FALSE, digits = 4)
  cat("\nNOTES\n-----\n")
  cat("- AT8 is a POSITIVE CONTROL: it stains the same tau species as PHF1.\n")
  cat("- The donor random-intercept LMM was verified as correctly calibrated by within-donor\n")
  cat("  permutation.\n")
  cat("- PHF1_Otsu is a per-donor Otsu threshold with a hard floor, applied to the PHF1 stain.\n")
  cat("\n"); print(sessionInfo())
  sink()
}

emit(res_ref, md_all, "ref_all_exc",
     "PHF1+ vs PHF1- state markers, all excitatory neurons, covariate-adjusted (IMC)",
     "value ~ grp + Sex + Age_s + PMI_s + (1 | patient_id)",
     c("PURPOSE: baseline for the subcluster-adjusted fit; REF vs SUB isolates the subcluster",
       "term."))

emit(res_sub, md_all, "exc_subcluster_adj",
     "PHF1+ vs PHF1- state markers, subcluster-adjusted (IMC)",
     "value ~ grp + subcluster + Sex + Age_s + PMI_s + (1 | patient_id)",
     c(paste0("Reference subcluster: ", levels(md_all$subcluster)[1]),
       "",
       "The contrast is now the WITHIN-subcluster PHF1+ vs PHF1- difference. A marker that shrinks",
       "substantially between REF and SUB was partly reflecting WHICH excitatory subtype bears",
       "tangles rather than a per-cell change in tangle-bearing neurons."))

emit(res_vuln, md_vuln, "reln_calb1",
     "PHF1+ vs PHF1- state markers, RELN+ and CALB1+ excitatory neurons (IMC)",
     "value ~ grp + Sex + Age_s + PMI_s + (1 | patient_id)",
     c("CELL SET: the two tangle-vulnerable entorhinal excitatory populations pooled --",
       paste0("  ", paste(VULN_CLUSTERS, collapse = " + ")),
       "",
       "RELN and CALB1 are TYPE markers, i.e. inputs to the PCA -> Harmony -> Phenograph pipeline",
       "that defined these clusters. That circularity affects the cluster DEFINITION, not this",
       "contrast: the outcomes modelled here are STATE markers, never clustering inputs."))

# -------------------------------------------------------------------
# Cross-fit comparison, so the three are readable side by side
# -------------------------------------------------------------------
cmp <- res_ref  %>% select(gene, logFC_ref = logFC, padj_ref = padj) %>%
  left_join(res_sub  %>% select(gene, logFC_sub = logFC, padj_sub = padj), by = "gene") %>%
  left_join(res_vuln %>% select(gene, logFC_vuln = logFC, padj_vuln = padj), by = "gene") %>%
  mutate(pct_change_sub_vs_ref = 100 * (logFC_sub - logFC_ref) / abs(logFC_ref)) %>%
  arrange(desc(abs(logFC_ref)))
wr(cmp, "source_data_fit_comparison.tsv")

writeLines(c(
  "PHF1+ vs PHF1- state markers: three fits side by side (IMC, AD cases, excitatory neurons).",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "REF  value ~ grp + Sex + Age_s + PMI_s + (1 | patient_id)                 all 4 Exc subclusters",
  "SUB  value ~ grp + subcluster + Sex + Age_s + PMI_s + (1 | patient_id)    all 4 Exc subclusters",
  "VULN value ~ grp + Sex + Age_s + PMI_s + (1 | patient_id)                 RELN + CALB1 only",
  "",
  "pct_change_sub_vs_ref is the change in log2FC when subcluster is held constant. A large",
  "negative value means the marker's PHF1+/- difference was substantially a composition effect --",
  "it reflected which excitatory subtype bears tangles rather than a per-cell change.",
  "",
  "VULN is a different cell set, so its log2FC is on its own PHF1- baseline (mu_neg differs) and",
  "is comparable in DIRECTION and rough magnitude, not to the third decimal.",
  "",
  "See the per-fit stats logs for model details and notes."),
  file.path(out_dir, "stats_fit_comparison.txt"))

cat("\n== Fit comparison (log2FC, PHF1+ / PHF1-) ==\n")
print(as.data.frame(cmp), row.names = FALSE, digits = 3)
for (nm in c("REF", "SUB", "VULN")) {
  r <- switch(nm, REF = res_ref, SUB = res_sub, VULN = res_vuln)
  cat(sprintf("%-5s n=%d cells, %d PHF1+, sig(|lfc|>=%.2f & padj<%.2f): %d/%d\n",
              nm, r$n_cells[1], r$n_pos[1],
              LFC_MIN, FDR, sum(r$padj < FDR & abs(r$logFC) >= LFC_MIN), nrow(r)))
}
