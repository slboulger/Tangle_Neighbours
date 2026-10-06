# ---------------------------------------------------------------------------
# decompose_tangle_response_phf1.r
#
# Supporting analysis for Fig. 4C-E - not drawn in a panel
# Threshold-based decomposition of the same two signatures; rrho2_decomposition_phf1.r
# cross-tabulates its gene classes against the RRHO2 quadrants.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# =============================================================================
# decompose_tangle_response_phf1.r
#
# Decompose the CELL-AUTONOMOUS (Fig 2) vs NON-CELL-AUTONOMOUS (Fig 3) tangle
# response into (i) a shared / propagating pre-tangle program and (ii) a
# proximity-unique, genuinely non-cell-autonomous program (primary deliverable).
# This is NOT a convergence pass/fail test.
#
# Reuses the two existing per-gene result tables (NO recompute):
#   Fig 2 cell-autonomous : deg/de_cat/de_PHF1/{ct}/{ct}_PHF1TRUEVsPHF1FALSE.tsv
#     effect  = logFC (log2), contrast "PHF1TRUE - PHF1FALSE"
#     SIGN    : POSITIVE = up in the tangle-bearing (PHF1+) cell.  (already correct)
#   Fig 3 non-cell-autonomous : deg/de_linear_distance/{ct}/{ct}_dist_to_phf1_um_scaled.tsv
#     effect  = logFC (dream coefficient on sd-scaled *log* distance to nearest PHF1+ neuron)
#     SIGN    : NEGATIVE raw slope = expression rises as distance shrinks = up NEAR tangles.
#            -> we NEGATE the raw slope so that, after the flip, POSITIVE = tangle-associated
#               (up near tangles) in BOTH analyses. See the explicit comment in load_effects().
#
# Subtype panel: neuronal subtypes that have BOTH a cell-autonomous and a proximity
# table AND a cell-autonomous program (>0 sig pathways). Exc-IT-L2-3-CBLN2-HOPX is the
# FOCUS; the other member of PANEL (below) is the cross-subtype comparison.
#
# Outputs: plots/decompose_tangle_response_phf1/  (reads deg/; does not write to it).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(ggrepel); library(stringr); library(purrr)
})

# ---- project root: HPC canonical (per house convention) with local-mount fallback
#      so the same script runs on the HPC and on the mounted volume. ----
proj_hpc   <- "<PROJECT_ROOT>/phf1_v2"
proj_local <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(proj_hpc)) proj_hpc else proj_local)

source("R/palettes.R")                                  # celltype_palette, fig_theme, neuron_order ...
source("R/pathway_enrichment_enrichR_with_background.r") # pathway_analysis_enrichr(), format_res_table_enrichr()

set.seed(1)                                             # anything stochastic (label permutation)
options(timeout = 120)                                  # enrichR web calls

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a  # base R>=4.4 has this; define for older HPC R

# ------------------------------- config --------------------------------------
FOCUS  <- "Exc-IT-L2-3-CBLN2-HOPX"
# Two-subtype comparison: the vulnerable calbindin population and the L3-5 IT
# neuron. CBLN2 is the FOCUS (gets the apoptosis/CCR/spillover deep-dive); both get the
# core autonomous-vs-non decomposition (correlation, quadrant classes, GO, signed-OR quadrant).
PANEL  <- c("Exc-IT-L2-3-CBLN2-HOPX",       # FOCUS (vulnerable calbindin population)
            "Exc-IT-L3-5-CHGA-IL1RAPL2")

CA_PADJ <- 0.1     # Fig 2 significance: BH p < 0.1 AND ...
CA_LFC  <- 0.25    #                     ... |log2FC| > 0.25
PX_PADJ <- 0.1     # Fig 3 significance: BH p < 0.1 (no fold-change threshold; matches linear pipeline)
N_PERM  <- 10000   # gene-label permutations for the correlation null

GO_DBS  <- c("GO_Biological_Process_2025",
             "GO_Cellular_Component_2025",
             "GO_Molecular_Function_2025")
GO_MIN_SIZE <- 10; GO_MAX_SIZE <- 250; GO_MIN_OVERLAP <- 3   # match running_pathway_enrichment_w_background.R

# Saved under plots/ (house convention: named analysis dir, alongside the other figures).
OUT <- "plots/decompose_tangle_response_phf1"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fig2_path <- function(ct) file.path("deg/de_cat/de_PHF1", ct, paste0(ct, "_PHF1TRUEVsPHF1FALSE.tsv"))
fig3_path <- function(ct) file.path("deg/de_linear_distance", ct, paste0(ct, "_dist_to_phf1_um_scaled.tsv"))

class_palette <- c(
  "shared-concordant" = "#D7301F",
  "proximity-unique"  = "#2171B5",
  "tangle-unique"     = "#FDAE6B",
  "discordant"        = "#6A51A3",
  "ns"                = "grey82"
)
class_levels <- names(class_palette)

# Darker shades of the class colours, used for the repelled gene LABELS so text stays legible
# against the lighter dots (esp. tangle-unique orange).
label_palette <- c(
  "shared-concordant" = "#A50F15",
  "proximity-unique"  = "#08519C",
  "tangle-unique"     = "#A63603",
  "discordant"        = "#4A1486",
  "ns"                = "grey40")

# -------------------------- curated gene modules -----------------------------
# Apoptosis annotation. SOURCE: GO:0043065 (positive regulation of apoptotic process)
# vs GO:0043066 (negative regulation of apoptotic process), refined with the canonical
# BCL2-family / caspase / IAP literature (Singh, Letai & Olson, Nat Rev Mol Cell Biol 2019).
apoptosis_pro <- c(
  "BAX","BAK1","BAD","BID","BBC3","PMAIP1","BCL2L11","BIK","BMF","HRK",          # BH3-only / effectors
  "CASP3","CASP6","CASP7","CASP8","CASP9","CASP2","APAF1","CYCS","DIABLO",       # apoptosome / caspases
  "TP53","FAS","FASLG","FADD","TNFRSF10A","TNFRSF10B","TNFRSF1A","TRADD",        # death-receptor axis / p53
  "BNIP3","BNIP3L","AIFM1","ENDOG","PARP1","GSDME","BAX","STK17A","PDCD5")
apoptosis_anti <- c(
  "BCL2","BCL2L1","MCL1","BCL2L2","BCL2A1","BCL2L10",                            # pro-survival BCL2 family
  "BIRC2","BIRC3","BIRC5","XIAP","NAIP","CFLAR",                                 # IAPs / FLIP
  "AKT1","BAG1","BAG3","NFKB1","NFKBIA","CDKN1A","GADD45B","HSPA1A","HSPA1B",
  "SOD1","CLU","FAIM2","TNFAIP3")
# Cell-cycle-reentry / developmental-neurogenesis module (report MKI67 explicitly).
ccr_cellcycle <- c(
  "MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCND1","CCND2","CCNE1","CCNA2",
  "CDK1","CDK2","CDK4","CDK6","CDC20","CDC6","MCM2","MCM3","MCM5","MCM6",
  "AURKA","AURKB","BUB1","BUB1B","FOXM1","E2F1","PLK1","TYMS","RRM2","CDKN2A")
ccr_neurogenesis <- c(
  "DCX","SOX2","SOX4","SOX11","NES","NEUROD1","NEUROD2","TUBB3","ASCL1",
  "PAX6","EOMES","NCAM1","STMN1","STMN2","GAP43","TBR1","PROX1","CALB2")

# Glial / vascular markers (canonical class markers) -- spillover / ambient-RNA sanity check.
glial_markers <- unique(c(
  # Astro
  "GFAP","ALDH1L1","AQP4","SLC1A2","SLC1A3","GLUL","SPARCL1","S100B","VIM",
  # Micro
  "P2RY12","CX3CR1","TYROBP","C1QA","C1QB","C1QC","CSF1R","AIF1","ITGAM","PTPRC",
  # Oligo lineage / OPC
  "OLIG1","OLIG2","PDGFRA","CSPG4","MBP","PLP1","MOG","MAG","SOX10","CLDN11",
  # Vascular (Endo / VLMC)
  "PECAM1","VWF","KDR","CLDN5","FLT1","COL1A1","COL1A2","DCN","PDGFRB"))

# =============================================================================
# Task 1 -- load, sign-align, restrict to common tested-gene background
# =============================================================================
load_effects <- function(ct) {
  ca <- read.delim(fig2_path(ct), check.names = FALSE)
  px <- read.delim(fig3_path(ct), check.names = FALSE)

  ca2 <- ca %>% transmute(gene = as.character(gene),
                          ca_log2FC = as.numeric(logFC),
                          ca_pval   = as.numeric(pval),
                          ca_padj   = as.numeric(padj))

  # ---- SIGN FLIP (Fig 3) -------------------------------------------------
  # `logFC` is the dream coefficient on (sd-scaled log) DISTANCE. A NEGATIVE slope
  # means expression goes UP as distance DECREASES, i.e. higher NEAR tangles.
  # We negate it so that, like the Fig 2 log2FC, POSITIVE = tangle-associated
  # (up near tangles). prox_slope is in logCPM per SD of log-distance (the tested
  # coefficient; its padj applies to it). prox_slope_per_um is the same flip on the
  # per-unit back-transform, kept for reference only.
  px2 <- px %>% transmute(gene = as.character(gene),
                          prox_slope        = -as.numeric(logFC),
                          prox_slope_per_um = -as.numeric(logFC_per_unit),
                          prox_t            = -as.numeric(t),   # proximity-signed moderated t (+ = up near)
                          prox_pval         = as.numeric(pval),
                          prox_padj         = as.numeric(padj))

  n_ca <- nrow(ca2); n_px <- nrow(px2)
  joined <- inner_join(ca2, px2, by = "gene")

  # Per-analysis enrichment backgrounds (the tested/expressed universe of EACH analysis,
  # matching the published enrichR convention). Used so a gene set is enriched against the
  # universe of the analysis that defines its membership -- NOT the common intersection,
  # which is self-enriched and collapses odds ratios for the large cell-autonomous sets.
  attr(joined, "bg_ca") <- ca2$gene   # all cell-autonomously (Fig 2) tested genes
  attr(joined, "bg_px") <- px2$gene   # all proximity (Fig 3) tested genes

  attr(joined, "prov") <- tibble(
    celltype       = ct,
    n_tested_ca    = n_ca,
    n_tested_px    = n_px,
    n_common       = nrow(joined),
    n_dropped_ca_only = n_ca - nrow(joined),   # tested cell-autonomously but not in proximity bg
    n_dropped_px_only = n_px - nrow(joined))   # tested in proximity but not cell-autonomously
  joined
}

# significance flags + quadrant class (Task 3)
classify <- function(df) {
  df %>% mutate(
    ca_sig    = !is.na(ca_padj)   & ca_padj   < CA_PADJ & abs(ca_log2FC) > CA_LFC,
    prox_sig  = !is.na(prox_padj) & prox_padj < PX_PADJ,
    same_sign = sign(ca_log2FC) == sign(prox_slope),
    class = factor(case_when(
      ca_sig & prox_sig &  same_sign ~ "shared-concordant",
      ca_sig & prox_sig & !same_sign ~ "discordant",
      prox_sig & !ca_sig             ~ "proximity-unique",
      ca_sig & !prox_sig             ~ "tangle-unique",
      TRUE                           ~ "ns"), levels = class_levels))
}

# =============================================================================
# Task 2 -- Spearman/Pearson + gene-label permutation p (Spearman)
# =============================================================================
# Spearman via Pearson-on-ranks: rank once, permute the ranks -> fast & identical to
# cor(method="spearman") (average-rank tie handling).
perm_spearman <- function(x, y, n_perm = N_PERM) {
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  rx <- rank(x); ry <- rank(y)
  obs <- suppressWarnings(cor(rx, ry))
  if (!is.finite(obs) || length(rx) < 3)
    return(list(spearman = obs, perm_p = NA_real_, n = length(rx)))
  perm <- replicate(n_perm, suppressWarnings(cor(rx, sample(ry))))
  list(spearman = obs,
       perm_p   = (1 + sum(abs(perm) >= abs(obs), na.rm = TRUE)) / (1 + n_perm),
       n        = length(rx))
}

corr_summary <- function(df, ct) {
  mk <- function(d, set) {
    sp <- perm_spearman(d$ca_log2FC, d$prox_slope)
    tibble(celltype = ct, set = set, n = sp$n,
           spearman = sp$spearman,
           pearson  = suppressWarnings(cor(d$ca_log2FC, d$prox_slope,
                                           method = "pearson", use = "complete.obs")),
           perm_p   = sp$perm_p)
  }
  bind_rows(
    mk(df, "all_shared_genes"),
    mk(dplyr::filter(df, ca_sig | prox_sig), "union_sig_either"))
}

# =============================================================================
# Task 4 -- GO ORA via the existing enrichR engine (per-celltype expressed-gene
# background = the common tested-gene set), returning a tidy combined table.
# =============================================================================
run_go <- function(genes, background) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) < GO_MIN_OVERLAP) return(tibble())
  res <- tryCatch(
    pathway_analysis_enrichr(interest_gene = genes,
                             enrichment_database = GO_DBS,
                             min_size = GO_MIN_SIZE, max_size = GO_MAX_SIZE,
                             min_overlap = GO_MIN_OVERLAP, background = background),
    error = function(e) { message("  enrichR failed: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(tibble())
  res$plot <- NULL
  dfs <- purrr::keep(res, is.data.frame)
  if (length(dfs) == 0) return(tibble())
  bind_rows(lapply(dfs, as_tibble)) %>% dplyr::select(-any_of("clusters"))
}

# Build the signed-OR quadrant. Each enriched GO term is placed at
#   x = cell-autonomous signed odds ratio (+OR if enriched among PHF1-UP genes,
#       -OR if enriched among PHF1-DOWN genes)
#   y = proximity signed odds ratio       (+OR if enriched among UP-NEAR genes,
#       -OR if among DOWN-NEAR genes)
# Directional GO is run on the two PUBLISHED signatures (all significant genes per
# analysis, split by direction), so the axes use the same significance as Fig 2 / Fig 3
# and the OR units match the lab's standard signed-OR dotplots. A term absent from a
# dimension sits on that axis (signed OR = 0). If a term is enriched in BOTH directions
# of one analysis (rare), the stronger-FDR direction wins. The dot's fill (FDR) and size
# (overlap) come from the term's strongest (min-FDR) directional enrichment. The four
# quadrants then read directly as shared-concordant / discordant / proximity-only /
# tangle-only.
signed_or_positions <- function(ca_up, ca_down, prox_up, prox_down) {
  tag <- function(g, dim, sgn) if (nrow(g)) dplyr::transmute(
    g, geneset, description, database, dim = dim,
    signed_or = sgn * odds_ratio, FDR, overlap) else tibble()
  all <- bind_rows(tag(ca_up, "ca", 1), tag(ca_down, "ca", -1),
                   tag(prox_up, "prox", 1), tag(prox_down, "prox", -1))
  if (nrow(all) == 0) return(tibble())
  best <- all %>% group_by(geneset, description, dim) %>%   # stronger-FDR direction per dim
    slice_min(FDR, n = 1, with_ties = FALSE) %>% ungroup()
  ca_part <- best %>% filter(dim == "ca") %>%
    transmute(geneset, description, database, x_or = signed_or, ca_FDR = FDR, ca_overlap = overlap)
  px_part <- best %>% filter(dim == "prox") %>%
    transmute(geneset, description, database, y_or = signed_or, px_FDR = FDR, px_overlap = overlap)
  full_join(ca_part, px_part, by = c("geneset", "description")) %>%
    mutate(database = dplyr::coalesce(database.x, database.y),
           x_or = dplyr::coalesce(x_or, 0), y_or = dplyr::coalesce(y_or, 0),
           FDR  = pmin(ca_FDR, px_FDR, na.rm = TRUE),
           overlap = ifelse(!is.na(ca_FDR) & (is.na(px_FDR) | ca_FDR <= px_FDR),
                            ca_overlap, px_overlap),
           region = dplyr::case_when(
             x_or != 0 & y_or != 0 & sign(x_or) == sign(y_or) ~ "shared-concordant",
             x_or != 0 & y_or != 0 & sign(x_or) != sign(y_or) ~ "discordant",
             x_or == 0 & y_or != 0                            ~ "proximity-only",
             x_or != 0 & y_or == 0                            ~ "tangle-only",
             TRUE                                             ~ "other")) %>%
    dplyr::select(geneset, description, database, x_or, y_or, overlap, FDR, region,
                  ca_FDR, px_FDR, ca_overlap, px_overlap)
}

# ------------------------- small triple-writer -------------------------------
# Plot PDF, source-data TSV, and stats TXT are written as SEPARATE files: no stats
# text is drawn on the figure itself (it lives in stats_<slug>.txt). Sizes in cm.
write_triple <- function(plot, source_df, stats_lines, dir, slug,
                          width_cm, height_cm) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0("plot_", slug, ".pdf")), plot,
         width = width_cm, height = height_cm, units = "cm", device = "pdf")
  write.table(source_df, file.path(dir, paste0("source_data_", slug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(stats_lines, file.path(dir, paste0("stats_", slug, ".txt")))
}

# =============================================================================
# Per-subtype decomposition (Tasks 1-4 + the two plots), run for every panel member
# =============================================================================
decompose_subtype <- function(ct) {
  message("== ", ct, " ==")
  ct_dir <- file.path(OUT, ct)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  je     <- load_effects(ct)
  prov   <- attr(je, "prov")
  joined <- classify(je)

  # ---- gene-level joined CSV (Task 1/3 deliverable) ----
  joined_out <- joined %>%
    transmute(gene, ca_log2FC, ca_pval, ca_padj,
              prox_slope, prox_slope_per_um, prox_pval, prox_padj,
              ca_sig, prox_sig, class) %>%
    arrange(class, desc(abs(prox_slope)))
  write.csv(joined_out,
            file.path(ct_dir, paste0("joined_gene_level_", ct, ".csv")),
            row.names = FALSE)

  # ---- quadrant class counts (Task 3) ----
  class_counts <- joined %>% count(class, .drop = FALSE, name = "n_genes")

  # ---- correlations + permutation (Task 2) ----
  corr <- corr_summary(joined, ct)

  # ---- GO on ALL FOUR relationship classes (Task 4), each split by direction ----
  # (proximity split by proximity direction; tangle-unique split by cell-autonomous
  #  direction since proximity is not significant there).
  # Enrichment background is per-analysis (see load_effects): tangle-unique against the
  # cell-autonomous universe; all proximity-defined classes against the proximity universe.
  bg_ca <- attr(je, "bg_ca"); bg_px <- attr(je, "bg_px")
  pick_bg <- function(nm) if (grepl("^tangle_unique", nm)) bg_ca else bg_px
  sets <- list(
    proximity_unique_up    = joined %>% filter(class == "proximity-unique",  prox_slope > 0) %>% pull(gene),
    proximity_unique_down  = joined %>% filter(class == "proximity-unique",  prox_slope < 0) %>% pull(gene),
    tangle_unique_up       = joined %>% filter(class == "tangle-unique",     ca_log2FC  > 0) %>% pull(gene),
    tangle_unique_down     = joined %>% filter(class == "tangle-unique",     ca_log2FC  < 0) %>% pull(gene),
    shared_concordant_up   = joined %>% filter(class == "shared-concordant", prox_slope > 0) %>% pull(gene),
    shared_concordant_down = joined %>% filter(class == "shared-concordant", prox_slope < 0) %>% pull(gene),
    discordant_upnear      = joined %>% filter(class == "discordant",        prox_slope > 0) %>% pull(gene),
    discordant_downnear    = joined %>% filter(class == "discordant",        prox_slope < 0) %>% pull(gene))

  go_tbls <- lapply(names(sets), function(nm) {
    message("  GO [class]: ", nm, " (", length(sets[[nm]]), " genes)")
    g <- run_go(sets[[nm]], pick_bg(nm))
    if (nrow(g)) { g$set <- nm; write.table(
      g, file.path(ct_dir, paste0("go_", nm, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE) }
    g
  })
  names(go_tbls) <- names(sets)

  # ---- directional GO on the PUBLISHED signatures -> signed-OR quadrant figure ----
  # ca_* against the cell-autonomous universe; prox_* against the proximity universe.
  dir_sets <- list(
    ca_up    = list(g = joined %>% filter(ca_sig,   ca_log2FC  > 0) %>% pull(gene), bg = bg_ca),
    ca_down  = list(g = joined %>% filter(ca_sig,   ca_log2FC  < 0) %>% pull(gene), bg = bg_ca),
    prox_up  = list(g = joined %>% filter(prox_sig, prox_slope > 0) %>% pull(gene), bg = bg_px),
    prox_down= list(g = joined %>% filter(prox_sig, prox_slope < 0) %>% pull(gene), bg = bg_px))
  dir_go <- lapply(names(dir_sets), function(nm) {
    message("  GO [dir]: ", nm, " (", length(dir_sets[[nm]]$g), " genes)")
    run_go(dir_sets[[nm]]$g, dir_sets[[nm]]$bg)
  })
  names(dir_go) <- names(dir_sets)
  pos <- signed_or_positions(dir_go$ca_up, dir_go$ca_down, dir_go$prox_up, dir_go$prox_down)

  list(ct = ct, joined = joined, joined_out = joined_out, prov = prov,
       class_counts = class_counts, corr = corr, go = go_tbls, dir_go = dir_go,
       pos = pos, ct_dir = ct_dir)
}

# ---------------------------- scatter (Task 3) -------------------------------
make_scatter <- function(d, ct, label_genes) {
  # A trend line summarises the association; the Spearman rho / p live in stats_<slug>.txt and
  # panel_correlations.tsv (not drawn on the plot).
  lab_df <- d %>% filter(gene %in% label_genes)
  ggplot(d, aes(ca_log2FC, prox_slope)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = c(-CA_LFC, CA_LFC), linetype = "dashed",
               linewidth = 0.25, colour = "grey75") +
    # points first, then the trend line (in front of dots), then labels (in front of line)
    geom_point(aes(colour = class), size = 0.7, alpha = 0.75) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey25", linewidth = 0.5) +
    ggrepel::geom_text_repel(data = lab_df, aes(label = gene),
                             colour = label_palette[as.character(lab_df$class)],
                             size = 2, max.overlaps = Inf, min.segment.length = 0,
                             segment.size = 0.2, show.legend = FALSE) +
    scale_colour_manual(values = class_palette, name = NULL, drop = FALSE) +
    labs(x = expression("Cell-autonomous log"[2] * "FC (+ = up in PHF1+ cell)"),
         y = "Non-cell-autonomous slope\n(logCPM per SD; + = up near tangles)") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8),
          aspect.ratio = 1,                     # square plotting panel
          legend.position = "right") +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)))
}

# ------------------- rank-rank map (gene-level concordance grid) --------------
# Genes are ranked on x by cell-autonomous log2FC and on y by proximity slope; each gene is one
# coloured square at (its x-rank, its y-rank) -- i.e. the "same-gene" cells of a gene x gene grid.
# Concordant genes lie on the diagonal (up/up = top-right, down/down = bottom-left); discordant on
# the anti-diagonal (top-left / bottom-right); proximity-unique form vertical arms (middle x,
# extreme y); tangle-unique form horizontal arms (extreme x, middle y). Colour = class.
make_rank_heatmap <- function(d) {
  dd <- d %>% mutate(rx = rank(ca_log2FC, ties.method = "first"),
                     ry = rank(prox_slope, ties.method = "first"),
                     class = factor(class, levels = class_levels))
  ggplot() +
    geom_abline(slope = nrow(dd) / nrow(dd), intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey70") +
    geom_point(data = dplyr::filter(dd, class == "ns"),  aes(rx, ry),
               colour = "grey88", size = 0.35, shape = 15) +
    geom_point(data = dplyr::filter(dd, class != "ns"),  aes(rx, ry, colour = class),
               size = 0.5, shape = 15) +
    scale_colour_manual(values = class_palette, name = NULL,
                        breaks = setdiff(class_levels, "ns")) +
    coord_equal() +
    labs(x = expression("Cell-autonomous log"[2] * "FC rank (low to high)"),
         y = "Proximity-slope rank (low to high)") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8), legend.position = "right") +
    guides(colour = guide_legend(override.aes = list(size = 1.8)))
}

# --------------------- signed-OR pathway quadrant dotplot (Task 4 headline) ----
# Axes are odds ratios (matching the lab's signed-OR dotplots): x = cell-autonomous
# signed OR, y = proximity signed OR. All enriched terms are drawn as dots; to keep the
# figure readable we LABEL every jointly-enriched term (shared-concordant + discordant --
# the off-axis points, which are rare) plus the top few most-significant proximity-only and
# tangle-only terms. Label colour encodes the region.
# --------- area-proportional 2-set Euler diagram (dependency-free) ------------
# Circle AREAS are proportional to set sizes; the centre distance is solved numerically
# so the lens (overlap) AREA is proportional to the shared count. n_a / n_b are the FULL
# set sizes (each already includes the shared items); n_shared is the intersection.
# Degrades gracefully to a single circle if one set is empty.
# shared_breakdown (optional c(concordant, discordant)) subdivides the overlap label so it
# does not overstate "shared": the intersection is genes/terms significant in BOTH analyses,
# but only the concordant fraction moves the SAME way in both (a shared program); discordant
# items move oppositely. When supplied, the lens is labelled "<total> (<conc> conc / <disc> disc)".
euler2 <- function(n_a, n_b, n_shared, label_a, label_b,
                   fill_a = "#D94801", fill_b = "#2171B5", shared_breakdown = NULL) {
  circ <- function(xc, r, grp) { t <- seq(0, 2 * pi, length.out = 240)
    data.frame(x = xc + r * cos(t), y = r * sin(t), grp = grp) }

  if (n_a <= 0 || n_b <= 0) {                        # one set empty -> single circle
    big  <- if (n_a >= n_b) list(n = n_a, lab = label_a, fill = fill_a, other = label_b)
            else            list(n = n_b, lab = label_b, fill = fill_b, other = label_a)
    r <- sqrt(max(big$n, 1) / pi)
    return(ggplot() +
      geom_polygon(data = circ(0, r, "x"), aes(x, y), fill = big$fill,
                   alpha = 0.45, colour = "black", linewidth = 0.3) +
      annotate("text", x = 0, y = 0, label = format(big$n, big.mark = ","), size = 3) +
      annotate("text", x = 0, y = r * 1.18,
               label = paste0(big$lab, "  (", big$other, " = 0)"),
               size = 2.5, fontface = "bold") +
      coord_equal(clip = "off") + theme_void(base_size = 8) +
      theme(plot.margin = margin(12, 12, 12, 12)))
  }

  r_a <- sqrt(n_a / pi); r_b <- sqrt(n_b / pi)
  lens <- function(d) {
    if (d >= r_a + r_b) return(0)
    if (d <= abs(r_a - r_b)) return(pi * min(r_a, r_b)^2)
    r_a^2 * acos((d^2 + r_a^2 - r_b^2) / (2 * d * r_a)) +
      r_b^2 * acos((d^2 + r_b^2 - r_a^2) / (2 * d * r_b)) -
      0.5 * sqrt((-d + r_a + r_b) * (d + r_a - r_b) * (d - r_a + r_b) * (d + r_a + r_b))
  }
  d <- if (n_shared <= 0) r_a + r_b + 0.08 * (r_a + r_b)         # disjoint (small gap)
       else if (n_shared >= min(n_a, n_b)) abs(r_a - r_b)         # one set subset of other
       else uniroot(function(x) lens(x) - n_shared,
                    lower = abs(r_a - r_b) + 1e-6, upper = r_a + r_b - 1e-6)$root
  xa <- 0; xb <- d
  polys <- rbind(circ(xa, r_a, "a"), circ(xb, r_b, "b"))
  x_share <- (r_a - r_b + d) / 2                                 # radical-line x (lens centre)
  shared_lab <- if (!is.null(shared_breakdown) && length(shared_breakdown) == 2)
    sprintf("%s\n(%s conc / %s disc)", format(n_shared, big.mark = ","),
            format(shared_breakdown[1], big.mark = ","),
            format(shared_breakdown[2], big.mark = ","))
  else format(n_shared, big.mark = ",")
  count_df <- data.frame(
    x = c(xa - 0.45 * r_a, x_share, xb + 0.45 * r_b), y = 0,
    lab = c(format(n_a - n_shared, big.mark = ","), shared_lab,
            format(n_b - n_shared, big.mark = ",")))
  title_df <- data.frame(x = c(xa - 0.4 * r_a, xb + 0.4 * r_b),
                         y = max(r_a, r_b) * 1.18, lab = c(label_a, label_b))
  ggplot() +
    geom_polygon(data = polys, aes(x, y, fill = grp, group = grp),
                 alpha = 0.45, colour = "black", linewidth = 0.3) +
    geom_text(data = count_df, aes(x, y, label = lab), size = 3) +
    geom_text(data = title_df, aes(x, y, label = lab), size = 2.5, fontface = "bold") +
    scale_fill_manual(values = c(a = fill_a, b = fill_b), guide = "none") +
    coord_equal(clip = "off") + theme_void(base_size = 8) +
    theme(plot.margin = margin(12, 12, 12, 12))
}

# ---- paired signed-OR plot of the SHARED pathways (enriched in BOTH measures) ----
# Each shared GO term is a row with two points: its signed odds ratio in the PHF1+/- measure
# and in the distance measure, joined by a line. Concordant terms keep both points on the same
# side of 0; discordant terms straddle 0 (the connecting line crosses the dashed zero axis).
# Custom legend key: draw a circle AND a triangle side by side at the mapped size, so each
# Overlap key row reads "circle, triangle, <number>" left-to-right (the sizing applies to both
# the PHF1 (circle) and distance (triangle) points).
draw_key_ct <- function(data, params, size) {
  ppm <- 72.27 / 25.4                                  # ggplot size (mm) -> grid fontsize (pt)
  fs  <- (data$size %||% 1.5) * ppm
  grid::gTree(children = grid::gList(
    grid::pointsGrob(0.20, 0.5, pch = 16, gp = grid::gpar(col = "grey20", fontsize = fs)),
    grid::pointsGrob(0.62, 0.5, pch = 17, gp = grid::gpar(col = "grey20", fontsize = fs))))
}

make_shared_pathway_plot <- function(shared_path) {
  if (nrow(shared_path) == 0) return(NULL)
  ord  <- shared_path$description[order(shared_path$distance_signed_OR)]
  # long form: one row per (pathway x measure), each with its OWN signed OR and gene overlap
  long <- bind_rows(
    shared_path %>% transmute(description, region, measure = "PHF1",
                              signed_or = phf1_signed_OR, overlap = phf1_overlap),
    shared_path %>% transmute(description, region, measure = "distance",
                              signed_or = distance_signed_OR, overlap = distance_overlap)) %>%
    mutate(description = factor(description, levels = ord))
  # Overlap size mapping fixed to match plot_dotplot_pathway_cat_phf1.r so circle sizes are
  # comparable across the project's pathway figures; show only the breaks spanning this plot.
  overlap_limits <- c(3, 40)
  overlap_breaks <- c(3, 5, 10, 20, 30)
  ov_min <- min(long$overlap, na.rm = TRUE); ov_max <- max(long$overlap, na.rm = TRUE)
  lo_i <- suppressWarnings(max(which(overlap_breaks <= ov_min)))
  hi_i <- suppressWarnings(min(which(overlap_breaks >= ov_max)))
  if (!is.finite(lo_i)) lo_i <- 1
  if (!is.finite(hi_i)) hi_i <- length(overlap_breaks)
  size_breaks <- overlap_breaks[lo_i:hi_i]
  ggplot(long, aes(signed_or, description)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    geom_line(aes(group = description, colour = region), linewidth = 0.5) +
    # real points: sized by overlap; contribute the Contrast (shape) + region (colour) legends
    geom_point(aes(shape = measure, colour = region, size = overlap),
               show.legend = c(size = FALSE, shape = TRUE, colour = TRUE)) +
    # invisible dummy: supplies ONLY the Overlap size legend, drawn as circle+triangle per row
    # (explicitly excluded from the shape/colour guides so its glyph does not leak into them)
    geom_point(aes(size = overlap), shape = 16, colour = "grey20", alpha = 0,
               show.legend = c(size = TRUE, shape = FALSE, colour = FALSE),
               key_glyph = draw_key_ct) +
    scale_colour_manual(values = c("shared-concordant" = "#D7301F", "discordant" = "#6A51A3"),
                        name = NULL) +
    scale_shape_manual(values = c(PHF1 = 16, distance = 17), name = "Contrast") +
    scale_size(name = "Overlap", range = c(2, 6), limits = overlap_limits, breaks = size_breaks,
               guide = guide_legend(order = 1, override.aes = list(alpha = 1))) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 28)) +
    labs(x = "Odds ratio (UP +, DOWN -)", y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8), legend.position = "right",
          legend.spacing.y   = grid::unit(4, "pt"),
          legend.box.spacing = grid::unit(3, "pt"),
          legend.key.width   = grid::unit(40, "pt"),
          legend.margin      = margin(1, 1, 1, 1)) +
    guides(shape  = guide_legend(order = 2, override.aes = list(size = 2.5)),
           colour = guide_legend(order = 3))
}

# ---- Set-level test: is the cell-autonomous PHF1+ signature recapitulated over proximity? ----
# This is the rigorous analogue of the "PHF1 module score vs distance" model. Genes are ranked
# by the proximity moderated t-statistic (+ = up near tangles); we then test whether the PHF1+
# UP and DOWN signatures (cell-autonomous DEGs) are coordinately enriched in that ranking.
#   - cameraPR : COMPETITIVE (signature vs the rest of the transcriptome; allows for inter-gene
#                correlation) -- the direct analogue of a background-subtracted module score.
#   - fgsea    : competitive enrichment score + leading-edge genes driving the convergence.
#   - t-test vs 0 : SELF-CONTAINED (absolute shift of the signature's proximity stats).
# The competitive-vs-self-contained contrast is informative: competitive-significant but not
# self-contained = the signature rises RELATIVE to a global downshift (as a module score would).
# n_core caps the PHF1+ signature at the strongest markers per direction (by |log2FC|), so it
# behaves like a marker-based module score and GSEA stays meaningful (the full cell-autonomous
# DEG set is ~half the transcriptome, which is invalid for fgsea though fine for cameraPR).
signature_convergence <- function(jf, ct, out_dir, n_core = 150) {
  ranks <- jf$prox_t; names(ranks) <- jf$gene
  ranks <- ranks[is.finite(ranks)]
  core <- function(up) {
    g <- jf %>% filter(ca_sig, if (up) ca_log2FC > 0 else ca_log2FC < 0) %>%
      arrange(desc(abs(ca_log2FC))) %>% pull(gene)
    intersect(head(g, n_core), names(ranks))
  }
  sets <- list(PHF1_up = core(TRUE), PHF1_down = core(FALSE))
  sets <- sets[vapply(sets, length, 0L) >= 5]
  if (length(sets) == 0) return(NULL)
  cam <- limma::cameraPR(ranks, limma::ids2indices(sets, names(ranks)))
  cam$set <- rownames(cam)
  fg <- fgsea::fgsea(pathways = sets, stats = ranks, minSize = 5, maxSize = length(ranks),
                     eps = 0, nPermSimple = 10000)
  sc <- purrr::map_dfr(names(sets), function(nm) {
    v <- ranks[sets[[nm]]]
    tibble(set = nm, n_genes = length(v), mean_prox_t = mean(v),
           selfcontained_p = t.test(v, mu = 0)$p.value)
  })
  out <- sc %>%
    left_join(tibble(set = cam$set, camera_direction = cam$Direction,
                     camera_p = cam$PValue, camera_fdr = cam$FDR), by = "set") %>%
    left_join(tibble(set = fg$pathway, fgsea_NES = fg$NES, fgsea_p = fg$pval, fgsea_padj = fg$padj,
                     leading_edge = vapply(fg$leadingEdge, paste, "", collapse = ";")), by = "set")
  write.table(out %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))),
              file.path(out_dir, paste0("phf1_signature_convergence_", ct, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  list(table = out, ranks = ranks, sets = sets)
}

# ---- CAMERA visualisation: proximity-statistic distribution, signature vs background ----
# cameraPR is a COMPETITIVE test (is the signature's proximity statistic shifted relative to the
# rest of the transcriptome?). The faithful visual is the distribution of the proximity t for the
# PHF1+ up / down markers against the background, with a barcode rug of the individual genes.
make_camera_plot <- function(conv) {
  if (is.null(conv)) return(NULL)
  df <- tibble(gene = names(conv$ranks), prox_t = as.numeric(conv$ranks)) %>%
    mutate(group = dplyr::case_when(gene %in% conv$sets$PHF1_up   ~ "PHF1+ up markers",
                                    gene %in% conv$sets$PHF1_down ~ "PHF1+ down markers",
                                    TRUE ~ "background"),
           group = factor(group, levels = c("PHF1+ up markers", "PHF1+ down markers", "background")))
  pal <- c("PHF1+ up markers" = "#D7301F", "PHF1+ down markers" = "#2171B5")
  ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_density(data = dplyr::filter(df, group == "background"), aes(prox_t),
                 fill = "grey88", colour = "grey60", linewidth = 0.3) +
    geom_density(data = dplyr::filter(df, group != "background"), aes(prox_t, colour = group),
                 linewidth = 0.6) +
    geom_rug(data = dplyr::filter(df, group == "PHF1+ up markers"), aes(prox_t),
             colour = "#D7301F", sides = "t", alpha = 0.5, length = grid::unit(4, "pt")) +
    geom_rug(data = dplyr::filter(df, group == "PHF1+ down markers"), aes(prox_t),
             colour = "#2171B5", sides = "b", alpha = 0.5, length = grid::unit(4, "pt")) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Proximity t-statistic  (+ = up near tangles)", y = "Density") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8), legend.position = "right")
}

# ---- divergence figures: the two arms of the proximity response ------------------
# The up/down-near-tangles response splits into a PHF1-CONVERGENT arm (shared-concordant:
# PHF1+ genes that also follow the proximity trend) and a distinct PROXIMITY-UNIQUE arm.
# These two functions contrast them functionally (GO) and at gene level (volcano).
program_go_df <- function(go_list, n_top = 8) {
  grab <- function(nm, program, dir) {
    g <- go_list[[nm]]
    if (is.null(g) || nrow(g) == 0) return(tibble())
    g %>% transmute(description, odds_ratio, overlap, FDR, program = program,
                    signed_or = if (dir == "up") odds_ratio else -odds_ratio)
  }
  bind_rows(
    grab("shared_concordant_up",   "PHF1-convergent",  "up"),
    grab("shared_concordant_down", "PHF1-convergent",  "down"),
    grab("proximity_unique_up",    "proximity-unique", "up"),
    grab("proximity_unique_down",  "proximity-unique", "down")) %>%
    group_by(program) %>% slice_min(FDR, n = n_top, with_ties = FALSE) %>% ungroup() %>%
    mutate(program = factor(program, levels = c("PHF1-convergent", "proximity-unique")))
}

make_program_dotplot <- function(df) {
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(description = stringr::str_wrap(description, 34))
  ob <- c(3, 5, 10, 20, 30); sb <- ob[ob >= min(df$overlap) & ob <= max(df$overlap)]
  if (length(sb) < 2) sb <- ob[ob <= max(df$overlap)]
  ggplot(df, aes(signed_or, reorder(description, signed_or))) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    geom_point(aes(fill = FDR, size = overlap), shape = 21, stroke = 0.3, colour = "black") +
    facet_grid(rows = vars(program), scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "navy", high = "gold", name = "FDR", transform = "log10",
                        guide = guide_colourbar(reverse = TRUE, order = 2)) +
    scale_size(name = "Overlap", range = c(2, 6), limits = c(3, 40), breaks = sb,
               guide = guide_legend(order = 1)) +
    labs(x = "Odds ratio (UP +, DOWN -)", y = NULL) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          strip.text = element_text(size = 7, face = "bold"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          panel.spacing = grid::unit(4, "pt"),
          plot.margin = margin(6, 8, 6, 8), legend.position = "right")
}

make_proximity_volcano <- function(jf, label_genes) {
  d <- jf %>% filter(prox_sig, class %in% c("shared-concordant", "proximity-unique", "discordant")) %>%
    mutate(origin = dplyr::recode(as.character(class),
             "shared-concordant" = "PHF1-convergent", "proximity-unique" = "proximity-unique",
             "discordant" = "discordant"),
           origin = factor(origin, levels = c("PHF1-convergent", "proximity-unique", "discordant")),
           neglogp = -log10(prox_padj))
  pal <- c("PHF1-convergent" = "#D7301F", "proximity-unique" = "#2171B5", "discordant" = "#6A51A3")
  lab <- d %>% filter(gene %in% label_genes)
  ggplot(d, aes(prox_slope, neglogp)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = origin), size = 0.7, alpha = 0.8) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene, colour = origin), size = 2,
                             max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
                             show.legend = FALSE) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Non-cell-autonomous slope (+ = up near tangles)",
         y = expression(-log[10] * " proximity FDR")) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(6, 8, 6, 8), legend.position = "right") +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)))
}

# =============================================================================
# RUN: decomposition for every panel subtype
# =============================================================================
results <- lapply(PANEL, decompose_subtype)
names(results) <- PANEL

# ---- panel-level diagnostics + correlations (Task 6) ----
panel_diag <- purrr::map_dfr(results, function(r) {
  cc_all <- r$corr %>% filter(set == "all_shared_genes")
  cc_uni <- r$corr %>% filter(set == "union_sig_either")
  n_by <- r$class_counts %>% deframe()
  prox_genes <- r$joined %>% filter(class == "proximity-unique") %>% pull(gene)
  n_bg <- nrow(r$joined)
  n_prox_glial <- sum(prox_genes %in% glial_markers)
  tibble(
    celltype        = r$ct,
    is_focus        = r$ct == FOCUS,
    n_tested_ca     = r$prov$n_tested_ca,
    n_tested_px     = r$prov$n_tested_px,
    n_common        = r$prov$n_common,
    n_dropped_ca_only = r$prov$n_dropped_ca_only,
    n_dropped_px_only = r$prov$n_dropped_px_only,
    n_ca_sig        = sum(r$joined$ca_sig),
    n_px_sig        = sum(r$joined$prox_sig),
    n_shared_concordant = as.integer(n_by[["shared-concordant"]] %||% 0),
    n_proximity_unique  = as.integer(n_by[["proximity-unique"]]  %||% 0),
    n_tangle_unique     = as.integer(n_by[["tangle-unique"]]     %||% 0),
    n_discordant        = as.integer(n_by[["discordant"]]        %||% 0),
    frac_proximity_unique = round(as.integer(n_by[["proximity-unique"]] %||% 0) / n_bg, 4),
    spearman_all    = round(cc_all$spearman, 3),
    perm_p_all      = cc_all$perm_p,
    spearman_union  = round(cc_uni$spearman, 3),
    perm_p_union    = cc_uni$perm_p,
    prox_unique_glial_frac = ifelse(length(prox_genes) > 0,
                                    round(n_prox_glial / length(prox_genes), 4), NA_real_),
    prox_unique_glial_genes = paste(intersect(prox_genes, glial_markers), collapse = ";"))
})
write.table(panel_diag, file.path(OUT, "panel_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

panel_corr <- purrr::map_dfr(results, ~ .x$corr)
write.table(panel_corr, file.path(OUT, "panel_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- per-subtype scatter + quadrant dotplot (all panel members) ----
for (r in results) {
  d <- r$joined
  # Label genes: top 2 in each direction for every class/colour (concordant, proximity-unique,
  # tangle-unique, discordant) -> ~16 labels. Ranked by the effect that defines each class
  # (proximity slope for proximity/concordant/discordant; cell-autonomous log2FC for tangle-unique).
  lab_genes <- unique(c(
    d %>% filter(class == "shared-concordant", prox_slope > 0) %>% slice_max(ca_log2FC + prox_slope, n = 2) %>% pull(gene),
    d %>% filter(class == "shared-concordant", prox_slope < 0) %>% slice_min(ca_log2FC + prox_slope, n = 2) %>% pull(gene),
    d %>% filter(class == "proximity-unique",  prox_slope > 0) %>% slice_max(prox_slope, n = 2) %>% pull(gene),
    d %>% filter(class == "proximity-unique",  prox_slope < 0) %>% slice_min(prox_slope, n = 2) %>% pull(gene),
    d %>% filter(class == "tangle-unique",     ca_log2FC  > 0) %>% slice_max(ca_log2FC,  n = 2) %>% pull(gene),
    d %>% filter(class == "tangle-unique",     ca_log2FC  < 0) %>% slice_min(ca_log2FC,  n = 2) %>% pull(gene),
    d %>% filter(class == "discordant",        prox_slope > 0) %>% slice_max(prox_slope, n = 2) %>% pull(gene),
    d %>% filter(class == "discordant",        prox_slope < 0) %>% slice_min(prox_slope, n = 2) %>% pull(gene)))

  p_sc <- make_scatter(d, r$ct, lab_genes)
  sc_src <- d %>% transmute(gene, ca_log2FC, prox_slope, ca_padj, prox_padj,
                            class, labelled = gene %in% lab_genes)
  cc <- r$corr
  sc_stats <- c(
    paste0("Decomposition scatter -- ", r$ct),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "x = cell-autonomous log2FC (Fig2, PHF1+ vs PHF1-; + = up in PHF1+ cell)",
    "y = non-cell-autonomous slope (Fig3 dream coef on sd-scaled log-distance, NEGATED",
    "    so + = up near tangles). Units: logCPM per SD of log-distance.",
    paste0("Sig thresholds: Fig2 padj<", CA_PADJ, " & |log2FC|>", CA_LFC,
           "; Fig3 padj<", PX_PADJ, "."),
    "",
    "Class counts:",
    capture.output(print(as.data.frame(r$class_counts))),
    "",
    paste0("Correlations (gene-label permutation, n_perm = ", N_PERM, "):"),
    capture.output(print(as.data.frame(cc))),
    "",
    "sessionInfo() saved alongside outputs (sessionInfo.txt).")
  write_triple(p_sc, sc_src, sc_stats, r$ct_dir,
               paste0("decomposition_scatter_", r$ct), 9, 7)   # cm (compact, square panel + side legend)

  # ---- rank-rank map: gene concordance grid (shared = diagonal, unique = arms) ----
  p_rr <- make_rank_heatmap(d)
  rr_src <- d %>% mutate(ca_rank = rank(ca_log2FC, ties.method = "first"),
                         prox_rank = rank(prox_slope, ties.method = "first")) %>%
    transmute(gene, ca_log2FC, prox_slope, ca_rank, prox_rank, class)
  rr_stats <- c(
    paste0("Rank-rank map -- ", r$ct),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Each gene = one square at (rank by cell-autonomous log2FC, rank by proximity slope).",
    "Diagonal = concordant ranking (shared response); anti-diagonal = discordant;",
    "vertical arms = proximity-unique; horizontal arms = tangle-unique. Colour = class.",
    paste0("n genes: ", nrow(d)))
  write_triple(p_rr, rr_src, rr_stats, r$ct_dir,
               paste0("rank_rank_map_", r$ct), 8.5, 7)   # cm (compact, square)

  # ---- GENE Venn (area-proportional): PHF1+/- DEGs vs distance DEGs ----
  # The overlap (significant in BOTH) is subdivided into concordant (same direction = a
  # shared program) vs discordant (opposite direction), so "shared" is not overstated.
  n_ca_sig <- sum(d$ca_sig); n_px_sig <- sum(d$prox_sig)
  n_gene_shared <- sum(d$ca_sig & d$prox_sig)                     # concordant + discordant
  n_gene_conc <- sum(d$class == "shared-concordant")
  n_gene_disc <- sum(d$class == "discordant")
  p_gv <- euler2(n_ca_sig, n_px_sig, n_gene_shared,
                 "PHF1+/- DEGs", "distance DEGs",
                 shared_breakdown = c(n_gene_conc, n_gene_disc))
  gv_src <- data.frame(
    measure = c("PHF1_only", "shared_total", "shared_concordant", "shared_discordant",
                "distance_only", "PHF1_total", "distance_total"),
    n_genes = c(n_ca_sig - n_gene_shared, n_gene_shared, n_gene_conc, n_gene_disc,
                n_px_sig - n_gene_shared, n_ca_sig, n_px_sig))
  gv_stats <- c(
    paste0("Gene Venn (area-proportional) -- ", r$ct),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Sets: genes significant in the PHF1+ vs - (cell-autonomous) analysis vs the distance",
    "(non-cell-autonomous) analysis. Circle areas and overlap scaled to gene counts.",
    "The overlap is split: concordant (same direction, a shared program) vs discordant",
    "(opposite direction). 'Shared' is NOT the whole overlap -- only the concordant part is",
    "a shared program.",
    "",
    paste0("PHF1 DEGs total:       ", n_ca_sig),
    paste0("distance DEGs total:   ", n_px_sig),
    paste0("shared / both (total): ", n_gene_shared,
           "  (concordant ", n_gene_conc, " / discordant ", n_gene_disc, ")"),
    paste0("PHF1-only:             ", n_ca_sig - n_gene_shared),
    paste0("distance-only:         ", n_px_sig - n_gene_shared))
  write_triple(p_gv, gv_src, gv_stats, r$ct_dir,
               paste0("gene_venn_", r$ct), 7, 6.5)     # cm (small)

  # ---- PATHWAY Venn + shared-pathway table: which pathways are shared? ----
  pos <- r$pos
  n_ca_path <- sum(pos$x_or != 0)                      # enriched in the PHF1 measure
  n_px_path <- sum(pos$y_or != 0)                      # enriched in the distance measure
  n_path_shared <- sum(pos$x_or != 0 & pos$y_or != 0)  # enriched in BOTH
  shared_path <- pos %>% filter(x_or != 0 & y_or != 0) %>%
    transmute(geneset, description, database, region,
              phf1_signed_OR = x_or, distance_signed_OR = y_or,
              phf1_overlap = ca_overlap, distance_overlap = px_overlap, FDR) %>%
    arrange(desc(abs(phf1_signed_OR) + abs(distance_signed_OR)))
  write.table(shared_path, file.path(r$ct_dir, paste0("shared_pathways_", r$ct, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  n_path_conc <- sum(pos$region == "shared-concordant")
  n_path_disc <- sum(pos$region == "discordant")
  p_pv <- euler2(n_ca_path, n_px_path, n_path_shared,
                 "PHF1 pathways", "distance pathways",
                 shared_breakdown = c(n_path_conc, n_path_disc))
  pv_src <- data.frame(
    measure = c("PHF1_only", "shared_total", "shared_concordant", "shared_discordant",
                "distance_only", "PHF1_total", "distance_total"),
    n_pathways = c(n_ca_path - n_path_shared, n_path_shared, n_path_conc, n_path_disc,
                   n_px_path - n_path_shared, n_ca_path, n_px_path))
  pv_stats <- c(
    paste0("Pathway Venn (area-proportional) -- ", r$ct),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Sets: GO terms enriched in the PHF1+ vs - measure vs the distance measure (directional",
    "GO on the published signatures; a term counts as enriched in a measure if it is",
    "significant in either direction of that measure). Circle areas/overlap scaled to counts.",
    paste0("GO libs: ", paste(GO_DBS, collapse = ", "),
           "; size ", GO_MIN_SIZE, "-", GO_MAX_SIZE,
           ", overlap>=", GO_MIN_OVERLAP, ", odds_ratio>2, BH FDR<=0.05."),
    "",
    paste0("PHF1 pathways total:     ", n_ca_path),
    paste0("distance pathways total: ", n_px_path),
    paste0("shared (both):           ", n_path_shared,
           "  -> see shared_pathways_", r$ct, ".tsv"),
    "",
    "Shared pathways:",
    if (nrow(shared_path)) capture.output(print(as.data.frame(
        shared_path[, c("description", "region", "phf1_signed_OR", "distance_signed_OR")]),
        row.names = FALSE)) else "  (none)")
  write_triple(p_pv, pv_src, pv_stats, r$ct_dir,
               paste0("pathway_venn_", r$ct), 7, 6.5)  # cm (small)

  # ---- paired signed-OR plot of the shared (concordant + discordant) pathways ----
  p_sp <- make_shared_pathway_plot(shared_path)
  if (!is.null(p_sp)) {
    sp_stats <- c(
      paste0("Shared pathways -- paired signed-OR plot -- ", r$ct),
      paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "Pathways enriched in BOTH the PHF1+/- and distance measures. Each term = two points",
      "(its signed odds ratio in each measure) joined by a line. Concordant terms keep both",
      "points on the same side of 0; discordant terms straddle 0 (line crosses the zero axis).",
      paste0("n shared pathways: ", nrow(shared_path),
             " (concordant ", n_path_conc, " / discordant ", n_path_disc, ")"))
    n_rows <- nrow(shared_path)
    write_triple(p_sp, shared_path, sp_stats, r$ct_dir,
                 paste0("shared_pathways_paired_", r$ct), 9.48, max(5, 1 + 0.8 * n_rows))  # cm
  } else {
    message("  (no shared pathways to plot for ", r$ct, ")")
  }
}

# =============================================================================
# Task 5 -- directionality checks on the FOCUS subtype
# =============================================================================
rf <- results[[FOCUS]]
jf <- rf$joined
focus_dir <- rf$ct_dir

# ---- 5a. Neuronal apoptosis: leading edge from proximity GO + curated module ----
# Leading-edge genes behind apoptosis / cell-death terms in the proximity-unique GO.
prox_go <- bind_rows(rf$go$proximity_unique_up, rf$go$proximity_unique_down)
apo_terms <- prox_go %>% filter(grepl("apopto|cell death|programmed cell death|neuron death",
                                       description, ignore.case = TRUE))
apo_leading <- if (nrow(apo_terms))
  unique(str_trim(unlist(str_split(apo_terms$genes, ";")))) else character(0)

# Curated-module view over the whole proximity signature (robust to GO not flagging it).
apo_tbl <- jf %>%
  filter(gene %in% c(apoptosis_pro, apoptosis_anti)) %>%
  mutate(apo_role = ifelse(gene %in% apoptosis_pro, "pro-apoptotic", "anti-apoptotic"),
         in_go_leading_edge = gene %in% apo_leading,
         rises_near_tangles = prox_slope > 0) %>%
  transmute(gene, apo_role, ca_log2FC, ca_sig, prox_slope, prox_padj, prox_sig,
            in_go_leading_edge, rises_near_tangles, class) %>%
  arrange(apo_role, desc(prox_slope))
write.table(apo_tbl, file.path(focus_dir, "directionality_apoptosis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Direction-aware split: tests whether PRO-death genes RISE near tangles (i.e. whether the
# PHF1- NEIGHBOURS carry the pro-apoptotic signal). The verdict is computed from the data.
apo_pro_up    <- apo_tbl %>% filter(apo_role == "pro-apoptotic",  prox_sig, prox_slope > 0)
apo_pro_down  <- apo_tbl %>% filter(apo_role == "pro-apoptotic",  prox_sig, prox_slope < 0)
apo_anti_up   <- apo_tbl %>% filter(apo_role == "anti-apoptotic", prox_sig, prox_slope > 0)
apo_anti_down <- apo_tbl %>% filter(apo_role == "anti-apoptotic", prox_sig, prox_slope < 0)
apo_ca_sig    <- apo_tbl %>% filter(ca_sig)                 # apoptosis in the tangle-bearer itself
apo_pro_sig_n <- nrow(apo_pro_up) + nrow(apo_pro_down)
apo_supports_neighbour <- nrow(apo_pro_up) > nrow(apo_pro_down) && nrow(apo_pro_up) >= 3
apo_verdict <- if (apo_pro_sig_n == 0) {
  "INCONCLUSIVE: no pro-apoptotic gene reaches significance in the proximity signature."
} else if (apo_supports_neighbour) {
  sprintf(paste0("SUPPORTED: %d of %d significant pro-apoptotic genes RISE near tangles, ",
                 "consistent with dying cells being the PHF1- neighbours."),
          nrow(apo_pro_up), apo_pro_sig_n)
} else {
  sprintf(paste0("NOT SUPPORTED: %d of %d significant pro-apoptotic genes FALL near tangles ",
                 "(higher further away), and %d apoptosis-annotated genes are significant in the ",
                 "CELL-AUTONOMOUS (tangle-bearer) signature. ",
                 "Note: neighbours already committed to death may have dropped out of the section ",
                 "(survivorship), which would also depress pro-death signal near tangles."),
          nrow(apo_pro_down), apo_pro_sig_n, nrow(apo_ca_sig))
}

# ---- 5b. Cell-cycle re-entry / neurogenesis module ----
ccr_tbl <- jf %>%
  filter(gene %in% c(ccr_cellcycle, ccr_neurogenesis)) %>%
  mutate(ccr_role = ifelse(gene %in% ccr_cellcycle, "cell-cycle", "neurogenesis"),
         rises_near_tangles = prox_slope > 0) %>%
  transmute(gene, ccr_role, ca_log2FC, ca_sig, prox_slope, prox_padj, prox_sig,
            rises_near_tangles, class) %>%
  arrange(ccr_role, desc(prox_slope))
write.table(ccr_tbl, file.path(focus_dir, "directionality_ccr.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ccr_sig  <- ccr_tbl %>% filter(prox_sig)
ccr_up   <- ccr_tbl %>% filter(prox_sig, prox_slope > 0)   # rising near tangles = reentry with proximity
ccr_down <- ccr_tbl %>% filter(prox_sig, prox_slope < 0)

# ---- 5c. Does the PHF1+ signature converge over proximity? (set-level tests) ----
conv <- signature_convergence(jf, FOCUS, focus_dir)
if (!is.null(conv)) {
  p_up <- fgsea::plotEnrichment(conv$sets$PHF1_up, conv$ranks) +
    labs(title = "PHF1+ UP signature", x = NULL, y = "Enrichment score") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(size = 8))
  p_dn <- fgsea::plotEnrichment(conv$sets$PHF1_down, conv$ranks) +
    labs(title = "PHF1+ DOWN signature",
         x = "Genes ranked by proximity t  (left = up near tangles)", y = "Enrichment score") +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title = element_text(size = 8))
  p_conv <- cowplot::plot_grid(p_up, p_dn, ncol = 1, align = "v")
  conv_stats <- c(
    paste0("PHF1+ signature convergence over proximity -- ", FOCUS),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Genes ranked by proximity moderated t (+ = up near tangles). cameraPR = competitive",
    "(module-score analogue); fgsea = competitive enrichment + leading edge; t-test vs 0 =",
    "self-contained (absolute shift). See phf1_signature_convergence_", "",
    capture.output(print(as.data.frame(conv$table[, c("set","n_genes","mean_prox_t",
      "camera_direction","camera_p","fgsea_NES","fgsea_p","selfcontained_p")]), row.names = FALSE)))
  write_triple(p_conv, conv$table[, setdiff(names(conv$table), "leading_edge")], conv_stats,
               focus_dir, paste0("phf1_signature_convergence_", FOCUS), 9, 9)  # cm

  # CAMERA competitive-test visualisation: signature vs background proximity-stat distribution
  p_cam <- make_camera_plot(conv)
  cam_src <- tibble(gene = names(conv$ranks), prox_t = as.numeric(conv$ranks)) %>%
    mutate(group = ifelse(gene %in% conv$sets$PHF1_up, "PHF1_up",
                   ifelse(gene %in% conv$sets$PHF1_down, "PHF1_down", "background")))
  cam_stats <- c(
    paste0("CAMERA competitive test -- proximity-stat distribution -- ", FOCUS),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Density of the proximity t-statistic (+ = up near tangles) for the PHF1+ up / down marker",
    "sets vs background. cameraPR asks whether each set is shifted relative to background:",
    capture.output(print(as.data.frame(conv$table[, c("set","n_genes","mean_prox_t",
      "camera_direction","camera_p","camera_fdr")]), row.names = FALSE)))
  write_triple(p_cam, cam_src, cam_stats, focus_dir,
               paste0("phf1_signature_camera_", FOCUS), 9, 6.5)  # cm
}

# ---- 5d. Divergence of the two proximity arms (convergent vs proximity-unique) ----
prog_df <- program_go_df(rf$go)
p_prog  <- make_program_dotplot(prog_df)
if (!is.null(p_prog)) {
  prog_h <- max(6, 1.5 + 0.5 * nrow(prog_df))
  prog_stats <- c(
    paste0("Program divergence GO dotplot -- ", FOCUS),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "Top GO terms of the two arms of the proximity response, side by side:",
    "  PHF1-convergent  = shared-concordant genes (PHF1+ genes that also follow the proximity trend)",
    "  proximity-unique = genes changing with proximity but NOT in the PHF1 signature",
    "x = signed odds ratio (UP +, DOWN -); fill = FDR; size = overlap. enrichR GO, BH FDR<=0.05.",
    "The two arms enrich distinct biology -- this is the functional divergence.")
  write_triple(p_prog, prog_df, prog_stats, focus_dir,
               paste0("program_divergence_go_", FOCUS), 9.48, prog_h)  # cm
}

# gene-level companion: proximity volcano coloured by origin
vlab <- unique(c("CLU", "PRNP", "MAP1B", "YWHAH", "HOPX", "GPM6A",
                 jf %>% filter(class == "proximity-unique") %>% slice_max(abs(prox_slope), n = 10) %>% pull(gene),
                 jf %>% filter(class == "discordant")       %>% slice_max(abs(prox_slope), n = 5)  %>% pull(gene)))
p_volc <- make_proximity_volcano(jf, vlab)
volc_src <- jf %>% filter(prox_sig, class %in% c("shared-concordant", "proximity-unique", "discordant")) %>%
  transmute(gene, prox_slope, prox_padj, neglog10_prox_fdr = -log10(prox_padj),
            origin = dplyr::recode(as.character(class),
              "shared-concordant" = "PHF1-convergent", "proximity-unique" = "proximity-unique",
              "discordant" = "discordant"),
            labelled = gene %in% vlab)
volc_stats <- c(
  paste0("Proximity volcano coloured by origin -- ", FOCUS),
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Proximity-significant genes (BH padj<0.1). x = non-cell-autonomous slope (+ = up near tangles);",
  "y = -log10 proximity FDR. Colour: PHF1-convergent (shared-concordant), proximity-unique, or",
  "discordant -- i.e. how much of the neighbourhood response echoes PHF1 vs is genuinely new.")
write_triple(p_volc, volc_src, volc_stats, focus_dir,
             paste0("proximity_volcano_", FOCUS), 10, 8)  # cm

# =============================================================================
# Task 7 -- glial-marker spillover in the FOCUS proximity-unique set
# =============================================================================
prox_unique_focus <- jf %>% filter(class == "proximity-unique") %>% pull(gene)
spill <- jf %>% filter(class == "proximity-unique", gene %in% glial_markers) %>%
  transmute(gene, prox_slope, prox_padj,
            direction = ifelse(prox_slope > 0, "up-near", "down-near"),
            flag = "GLIAL_MARKER_review_ambient_RNA") %>%
  arrange(desc(prox_slope))
write.table(spill, file.path(focus_dir, "spillover_glial.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
spill_frac <- if (length(prox_unique_focus) > 0)
  nrow(spill) / length(prox_unique_focus) else NA_real_

# =============================================================================
# SUMMARY.md -- summary of findings
# =============================================================================
fnum <- function(x, d = 2) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fp   <- function(x)        ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 2))
fd   <- rf$class_counts %>% deframe()
n_focus_bg <- nrow(jf)
cc_all_f <- rf$corr %>% filter(set == "all_shared_genes")
cc_uni_f <- rf$corr %>% filter(set == "union_sig_either")

# Which component dominates the focus? (report, don't pass/fail)
dominant <- if ((fd[["proximity-unique"]] %||% 0) > (fd[["shared-concordant"]] %||% 0))
  "PROXIMITY-UNIQUE (non-cell-autonomous)" else "SHARED-CONCORDANT (propagating)"

summ <- c(
  "# Decomposition of the tangle response -- summary",
  paste0("_Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), "._"),
  "",
  paste0("**Focus subtype:** ", FOCUS, " (vulnerable calbindin population). ",
         "**Panel:** ", paste(PANEL, collapse = ", "), "."),
  "",
  "Sign convention: cell-autonomous log2FC positive = up in the PHF1+ (tangle-bearing) cell;",
  "proximity slope is the Fig3 distance coefficient NEGATED, so positive = up near tangles.",
  "Both axes therefore point the same way (positive = tangle-associated).",
  "",
  "## 1. Common background (Task 1)",
  paste0("Focus: ", rf$prov$n_tested_ca, " genes tested cell-autonomously, ",
         rf$prov$n_tested_px, " tested for proximity; intersection = ",
         rf$prov$n_common, " genes (dropped ", rf$prov$n_dropped_ca_only,
         " CA-only + ", rf$prov$n_dropped_px_only, " proximity-only). ",
         "Gene classification and correlations use this common set; GO enrichment uses each ",
         "analysis's own tested universe as background (cell-autonomous vs proximity)."),
  "",
  "## 2. Global relationship (Task 2)",
  paste0("- All shared genes: Spearman rho = ", fnum(cc_all_f$spearman),
         " (permutation p = ", fp(cc_all_f$perm_p), ", n = ", cc_all_f$n, ")."),
  paste0("- Union sig-in-either: Spearman rho = ", fnum(cc_uni_f$spearman),
         " (permutation p = ", fp(cc_uni_f$perm_p), ", n = ", cc_uni_f$n, ")."),
  "- Cross-subtype benchmark in `panel_diagnostics.tsv`;",
  "  a modest shared correlation with a large proximity-unique residual is the expected",
  "  signature of a partly-propagating program plus a distinct non-cell-autonomous response.",
  "",
  "## 2b. PHF1+ signature convergence over proximity (set-level)",
  "Set-level test (the analogue of the PHF1 module-score-vs-distance model): genes ranked by the",
  "proximity t-statistic (+ = up near tangles); the top-150 cell-autonomous PHF1+ markers per",
  "direction tested for coordinated enrichment. Complements the genome-wide gene-level",
  "correlation in section 2.",
  if (exists("conv") && !is.null(conv)) {
    u <- conv$table[conv$table$set == "PHF1_up", ]
    dn <- conv$table[conv$table$set == "PHF1_down", ]
    c(paste0("- PHF1+ UP markers (n=", u$n_genes, "): cameraPR ", u$camera_direction,
             " p = ", fp(u$camera_p), "; self-contained mean prox-t = ", fnum(u$mean_prox_t),
             ", p = ", fp(u$selfcontained_p), " -> up near tangles (competitive AND absolute)."),
      paste0("- PHF1+ DOWN markers (n=", dn$n_genes, "): self-contained mean prox-t = ",
             fnum(dn$mean_prox_t), ", p = ", fp(dn$selfcontained_p), " -> down near tangles."),
      paste0("Interpretation: the core PHF1+ signature is tested at set level in PHF1- ",
             "neighbours (cf. the module score); the genome-wide gene-level ",
             "correlation is rho ~ ", fnum(cc_all_f$spearman), ", with the broad ",
             "proximity response also containing a distinct proximity-unique program. See `", FOCUS,
             "/phf1_signature_convergence_", FOCUS, ".tsv` + enrichment plot."))
  } else "- (not computed)",
  "",
  "## 3. Quadrant decomposition (Task 3)",
  paste0("- shared-concordant : ", fd[["shared-concordant"]] %||% 0),
  paste0("- proximity-unique  : ", fd[["proximity-unique"]]  %||% 0,
         "   <- primary non-cell-autonomous deliverable"),
  paste0("- tangle-unique     : ", fd[["tangle-unique"]]     %||% 0),
  paste0("- discordant        : ", fd[["discordant"]]        %||% 0, "  (flagged, retained)"),
  paste0("- not significant   : ", fd[["ns"]]                %||% 0,
         "   (of ", n_focus_bg, " common)"),
  paste0("**Dominant component: ", dominant, ".**"),
  "",
  "## 4. GO programs + pathway sharing (Task 4)",
  paste0("All FOUR relationship classes are GO-assessed (proximity-unique, tangle-unique, ",
         "shared-concordant, discordant), each split by direction: `", FOCUS, "/go_*.tsv`."),
  paste0("**Pathway sharing:** of the enriched GO terms, ", sum(rf$pos$x_or != 0),
         " are enriched in the PHF1+/- measure and ", sum(rf$pos$y_or != 0),
         " in the distance measure, with only ", sum(rf$pos$x_or != 0 & rf$pos$y_or != 0),
         " enriched in BOTH -- i.e. the two measures pick out almost entirely different ",
         "pathways. Shown as an area-proportional Venn (`", FOCUS, "/plot_pathway_venn_", FOCUS,
         ".pdf`, overlap split concordant/discordant); the shared terms are named in `", FOCUS,
         "/shared_pathways_", FOCUS, ".tsv` and drawn as a paired signed-OR plot (`", FOCUS,
         "/plot_shared_pathways_paired_", FOCUS, ".pdf`) that shows each term's direction in both ",
         "measures (concordant = same side of 0; discordant = straddles 0)."),
  paste0("An area-proportional gene Venn (`", FOCUS, "/plot_gene_venn_", FOCUS,
         ".pdf`) shows the PHF1 vs distance DEG overlap. IMPORTANT: the overlap is split into ",
         "concordant (", sum(rf$joined$class == "shared-concordant"), ", same direction = a shared ",
         "program) vs discordant (", sum(rf$joined$class == "discordant"), ", opposite direction), ",
         "so 'shared' is not overstated -- only the concordant fraction is a genuinely shared ",
         "program. The gene-level decomposition scatter remains the per-gene view."),
  "",
  "## 5. Directionality (Task 5)",
  "**Apoptosis** (annotation source: GO:0043065 positive- vs GO:0043066 negative-regulation of",
  "apoptotic process, refined with the canonical BCL2-family / caspase / IAP literature).",
  paste0("Pro-apoptotic genes significant near tangles: ", apo_pro_sig_n,
         " (RISE near tangles: ",
         ifelse(nrow(apo_pro_up) > 0, paste(apo_pro_up$gene, collapse = ", "), "none"),
         " | FALL near tangles: ",
         ifelse(nrow(apo_pro_down) > 0, paste(apo_pro_down$gene, collapse = ", "), "none"), ")."),
  paste0("Anti-apoptotic significant near tangles: ",
         ifelse(nrow(apo_anti_up) + nrow(apo_anti_down) > 0,
                paste0("RISE: ", ifelse(nrow(apo_anti_up) > 0, paste(apo_anti_up$gene, collapse = ", "), "none"),
                       " | FALL: ", ifelse(nrow(apo_anti_down) > 0, paste(apo_anti_down$gene, collapse = ", "), "none")),
                "none"), "."),
  paste0("Apoptosis genes significant in the CELL-AUTONOMOUS (tangle-bearer) signature: ",
         ifelse(nrow(apo_ca_sig) > 0, paste(apo_ca_sig$gene, collapse = ", "), "none"), "."),
  paste0("GO apoptosis/cell-death terms surviving in the proximity signature: ",
         ifelse(nrow(apo_terms) > 0, paste(unique(apo_terms$description), collapse = "; "),
                "none at FDR<=0.05"), "."),
  paste0("**Verdict on 'dying cells are the neighbours': ", apo_verdict, "**"),
  "",
  paste0("**Cell-cycle re-entry / neurogenesis.** Module genes significant near tangles: ",
         ifelse(nrow(ccr_sig) > 0,
                paste0(nrow(ccr_sig), " (RISE near tangles: ",
                       ifelse(nrow(ccr_up) > 0, paste(ccr_up$gene, collapse = ", "), "none"),
                       " | FALL: ",
                       ifelse(nrow(ccr_down) > 0, paste(ccr_down$gene, collapse = ", "), "none"), ")"),
                "none"), ". MKI67 status: ",
         { m <- jf %>% filter(gene == "MKI67")
           if (nrow(m) == 0) "not in tested background" else
             paste0("prox_slope=", fnum(m$prox_slope), ", prox_padj=", fp(m$prox_padj),
                    " (", ifelse(m$prox_sig, "SIG", "ns"), ")") }, "."),
  "",
  "## 6. Cross-subtype panel (Task 6)",
  "Full side-by-side in `panel_diagnostics.tsv` / `panel_correlations.tsv`.",
  "",
  "## 7. Glial spillover (Task 7)",
  paste0("Of ", length(prox_unique_focus), " focus proximity-unique genes, ",
         nrow(spill), " are canonical glial/vascular markers (",
         ifelse(is.na(spill_frac), "NA", paste0(round(100 * spill_frac, 1), "%")),
         "): ", ifelse(nrow(spill) > 0, paste(spill$gene, collapse = ", "), "none"), "."),
  paste0("These are NOT dropped -- flagged in `", FOCUS, "/spillover_glial.tsv` for a downstream ",
         "ambient-RNA / neighbouring-glia analysis."),
  "",
  "## Notes",
  "- Proximity slope is per SD of *log* distance (Fig3 canonical run used log_um); comparable",
  "  across the panel (all log_um) but not a per-micron effect.",
  "- Cross-subtype DEG counts are not directly comparable (different tested backgrounds);",
  "  proximity-unique is also reported as a fraction of the common background.")
writeLines(summ, file.path(OUT, "SUMMARY.md"))

# ---- provenance dump ----
sink(file.path(OUT, "sessionInfo.txt"))
cat("decompose_tangle_response_phf1.r\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("set.seed(1); N_PERM =", N_PERM, "\n")
cat("Panel:", paste(PANEL, collapse = ", "), "\n\n")
print(sessionInfo())
sink()

message("Done. Outputs in ", OUT)
