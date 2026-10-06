# ---------------------------------------------------------------------------
# rrho2_decomposition_phf1.r
#
# Figure panels: 4C
# Also writes the quadrant statistics, gene lists and GO tables read by the 4D, 4E and S6
# scripts. Run after the Table S4 and S6 differential expression.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# =============================================================================
# rrho2_decomposition_phf1.r
#
# Threshold-free Rank-Rank Hypergeometric Overlap (RRHO2) between the two per-gene
# signatures used in the tangle-response decomposition:
#   list1 = CELL-AUTONOMOUS  (PHF1+ vs PHF1- pseudobulk, Fig 2)
#   list2 = PROXIMITY        (expression vs distance to nearest PHF1+ neuron, Fig 3)
#
# Improves on the threshold-based decomposition (R/decompose_tangle_response_phf1.r) by giving
# a statistically-grounded, threshold-free map of concordance (up-up / down-down) vs discordance
# (up-down / down-up), plus data-driven quadrant gene sets that feed the existing GO pipeline
# and cross-validate the padj-based classes.
#
# Ranking metric = sign(effect) * -log10(pval)  [RRHO2 standard signed significance].
# SIGN CONVENTION (matches the decomposition):
#   cell-autonomous:  + = up in PHF1+ (tangle-bearing) cell
#   proximity:        + = up NEAR tangles  -> we NEGATE the distance-model sign (whose coef is
#                     on DISTANCE, + = up FAR from tangles) to make + = up near tangles.
# So RRHO2 uu = up in both = shared-concordant-up; dd = down in both; ud/du = discordant.
#
# RRHO2 is installed from GitHub:
#   Rscript -e 'remotes::install_github("RRHO2/RRHO2")'
#
# Outputs: plots/decompose_tangle_response_phf1/{ct}/rrho2/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(ggplot2); library(purrr)
})

proj_hpc   <- "<PROJECT_ROOT>/phf1_v2"
proj_local <- "<PROJECT_ROOT>/phf1_v2"
setwd(if (dir.exists(proj_hpc)) proj_hpc else proj_local)

source("R/palettes.R")
source("R/pathway_enrichment_enrichR_with_background.r")   # pathway_analysis_enrichr()

set.seed(1)
options(timeout = 180)

if (!requireNamespace("RRHO2", quietly = TRUE))
  stop("RRHO2 not installed. Install it with:\n",
       "  Rscript -e 'if(!requireNamespace(\"remotes\",quietly=TRUE)) install.packages(\"remotes\"); remotes::install_github(\"RRHO2/RRHO2\")'")

# ------------------------------- config --------------------------------------
SUBTYPES   <- c("Exc-IT-L2-3-CBLN2-HOPX", "Exc-IT-L3-5-CHGA-IL1RAPL2")
GO_DBS     <- c("GO_Biological_Process_2025", "GO_Cellular_Component_2025", "GO_Molecular_Function_2025")
DECOMP_DIR <- "plots/decompose_tangle_response_phf1"        # write rrho2/ alongside the decomposition
PFLOOR     <- .Machine$double.xmin                          # floor for -log10(pval) (avoid Inf)
# Heatmap colour scale (-log10 BY-corrected hypergeometric p of overlap):
HM_MIN     <- -log10(0.05)   # floor: overlaps below this read as background (non-significant)
HM_MAX     <- 10             # cap: so mildly/moderately significant overlap uses the gradient
                             #      instead of collapsing into the background of one saturating corner
                             #      (lower HM_MAX to emphasise mild signal further)
TOP_DISCORDANT <- 100        # ORA on the strongest N genes of the ud direction (up PHF1+/down near),
                             #      ranked by ca_t - prox_t (~4% of the per-celltype background).

fig2_path <- function(ct) file.path("deg/de_cat/de_PHF1", ct, paste0(ct, "_PHF1TRUEVsPHF1FALSE.tsv"))
fig3_path <- function(ct) file.path("deg/de_linear_distance", ct, paste0(ct, "_dist_to_phf1_um_scaled.tsv"))

# ---- ranked signed-significance metrics over the common tested-gene set -----
load_ranks <- function(ct) {
  ca <- read.delim(fig2_path(ct), check.names = FALSE)
  px <- read.delim(fig3_path(ct), check.names = FALSE)
  ca2 <- ca %>% transmute(gene = as.character(gene),
                          ca_metric = sign(logFC) * -log10(pmax(pval, PFLOOR)),
                          ca_t = t)                                   # limma moderated t (+ = up in PHF1+)
  # NEGATE distance sign so + = up near tangles (proximity-signed), matching the decomposition
  px2 <- px %>% transmute(gene = as.character(gene),
                          prox_metric = -sign(logFC) * -log10(pmax(pval, PFLOOR)),
                          prox_t = -t)                                # dream t, negated (+ = up near tangles)
  inner_join(ca2, px2, by = "gene") %>%
    filter(is.finite(ca_metric), is.finite(prox_metric), is.finite(ca_t), is.finite(prox_t)) %>%
    distinct(gene, .keep_all = TRUE)
}

# ---- GO via the existing enrichR engine (returns a tidy combined table) -----
run_go <- function(genes, background) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) < 3) return(tibble())
  res <- tryCatch(
    pathway_analysis_enrichr(interest_gene = genes, enrichment_database = GO_DBS,
                             min_size = 15, max_size = 250, min_overlap = 3, background = background),  # matches linear-distance pathway script
    error = function(e) { message("  enrichR failed: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(tibble())
  res$plot <- NULL
  dfs <- purrr::keep(res, is.data.frame)
  if (length(dfs) == 0) return(tibble())
  bind_rows(lapply(dfs, as_tibble)) %>% dplyr::select(-dplyr::any_of("clusters"))
}

# ---- defensive extractor for a RRHO2 quadrant overlap gene list -------------
# Quadrants: uu (up both), dd (down both), ud (up list1 / down list2), du (down list1 / up list2).
get_quadrant <- function(obj, q) {
  gl <- obj[[paste0("genelist_", q)]]
  if (is.null(gl)) return(character(0))
  ov <- gl[[paste0("gene_list_overlap_", q)]]
  if (is.null(ov)) {                                   # fall back if slot names differ by version
    cand <- gl[grepl("overlap", names(gl), ignore.case = TRUE)]
    ov <- if (length(cand)) cand[[1]] else unlist(gl, use.names = FALSE)
  }
  unique(as.character(ov))
}

# ---- per-quadrant set sizes (list1, list2, overlap) for a clean Euler diagram -----
quad_counts <- function(obj, q) {
  gl <- obj[[paste0("genelist_", q)]]
  glen <- function(pat) { i <- grep(pat, names(gl)); if (length(i)) length(gl[[i[1]]]) else 0L }
  c(n1 = glen("list1"), n2 = glen("list2"), no = glen("overlap"))   # n1/n2 include the overlap
}

# ---- area-proportional 2-set Euler diagram (dependency-free; ported from -------
#      decompose_tangle_response_phf1.r) -- replaces RRHO2_vennDiagram for clean, paper-ready venns
euler2 <- function(n_a, n_b, n_shared, label_a, label_b,
                   fill_a = "#D94801", fill_b = "#2171B5") {
  n_a <- as.numeric(n_a); n_b <- as.numeric(n_b); n_shared <- as.numeric(n_shared)
  circ <- function(xc, r, grp) { t <- seq(0, 2 * pi, length.out = 240)
    data.frame(x = xc + r * cos(t), y = r * sin(t), grp = grp) }
  if (n_a <= 0 || n_b <= 0) {
    big <- if (n_a >= n_b) list(n = n_a, lab = label_a, fill = fill_a, other = label_b)
           else            list(n = n_b, lab = label_b, fill = fill_b, other = label_a)
    r <- sqrt(max(big$n, 1) / pi)
    return(ggplot() +
      geom_polygon(data = circ(0, r, "x"), aes(x, y), fill = big$fill, alpha = 0.45,
                   colour = "black", linewidth = 0.3) +
      annotate("text", x = 0, y = 0, label = format(big$n, big.mark = ","), size = 3) +
      annotate("text", x = 0, y = r * 1.18, label = paste0(big$lab, "  (", big$other, " = 0)"),
               size = 2.5, fontface = "bold") +
      coord_equal(clip = "off") + theme_void(base_size = 8) + theme(plot.margin = margin(12, 12, 12, 12)))
  }
  r_a <- sqrt(n_a / pi); r_b <- sqrt(n_b / pi)
  lens <- function(d) {
    if (d >= r_a + r_b) return(0)
    if (d <= abs(r_a - r_b)) return(pi * min(r_a, r_b)^2)
    r_a^2 * acos((d^2 + r_a^2 - r_b^2) / (2 * d * r_a)) +
      r_b^2 * acos((d^2 + r_b^2 - r_a^2) / (2 * d * r_b)) -
      0.5 * sqrt((-d + r_a + r_b) * (d + r_a - r_b) * (d - r_a + r_b) * (d + r_a + r_b))
  }
  d <- if (n_shared <= 0) r_a + r_b + 0.08 * (r_a + r_b)
       else if (n_shared >= min(n_a, n_b)) abs(r_a - r_b)
       else uniroot(function(x) lens(x) - n_shared, lower = abs(r_a - r_b) + 1e-6,
                    upper = r_a + r_b - 1e-6)$root
  polys   <- rbind(circ(0, r_a, "a"), circ(d, r_b, "b"))
  x_share <- (r_a - r_b + d) / 2
  count_df <- data.frame(x = c(-0.45 * r_a, x_share, d + 0.45 * r_b), y = 0,
    lab = format(c(n_a - n_shared, n_shared, n_b - n_shared), big.mark = ","))
  title_df <- data.frame(x = c(-0.4 * r_a, d + 0.4 * r_b), y = max(r_a, r_b) * 1.18,
                         lab = c(label_a, label_b))
  ggplot() +
    geom_polygon(data = polys, aes(x, y, fill = grp, group = grp), alpha = 0.45,
                 colour = "black", linewidth = 0.3) +
    geom_text(data = count_df, aes(x, y, label = lab), size = 3) +
    geom_text(data = title_df, aes(x, y, label = lab), size = 2.5, fontface = "bold") +
    scale_fill_manual(values = c(a = fill_a, b = fill_b), guide = "none") +
    coord_equal(clip = "off") + theme_void(base_size = 8) + theme(plot.margin = margin(12, 12, 12, 12))
}

# ---- compact, gap-free, SMOOTH heatmap from RRHO2's hypergeometric matrix ------------
# RRHO2_heatmap can't be shrunk to a subfigure ("figure margins too large") and its inter-block
# white gaps aren't tunable, so we render obj$hypermat directly (bilinear-smoothed; colour capped
# [hm_min, hm_max] so mild/moderate overlap stays visible).
#
# ORIENTATION is DETERMINISTIC -- it does NOT depend on RRHO2's (version-dependent) internal index
# order. RRHO2 steps each list with a fixed stepsize, so each up/down stratum occupies a section of
# the axis PROPORTIONAL to its gene count. We (1) locate the all-NA stratum seam, (2) read the two
# section sizes, (3) map the LARGER section to whichever direction has MORE genes (n_up vs n_dn),
# then (4) flip so the UP-regulated end sits at the HIGH index of each axis -- exactly matching the
# down->up axis titles (x: down->up in PHF1+; y: down->up near tangles). The four corner intensities
# (largest-overlap quadrant = hottest corner) are computed as an INDEPENDENT cross-check, and the
# observed vs expected up-section fraction is recorded so a wrong flip cannot pass silently.
make_rrho_raster <- function(obj, quad, hm_min, hm_max,
                             n_up1, n_dn1, n_up2, n_dn2,
                             lab_x = "Cell-autonomous: down -> up in PHF1+",
                             lab_y = "Proximity: down -> up near tangles") {
  M <- obj$hypermat; n <- nrow(M)
  seam <- function(is_na) if (any(is_na)) mean(which(is_na)) else NA_real_
  bx_raw <- seam(apply(M, 1, function(r)  all(is.na(r))))
  by_raw <- seam(apply(M, 2, function(cl) all(is.na(cl))))

  # deterministic up-end per axis: larger gene stratum -> larger axis section
  up_end <- function(b_raw, n_up, n_dn) {
    if (is.na(b_raw) || n_up == n_dn) return(NA_character_)
    big_side <- if ((n - b_raw) >= b_raw) "h" else "l"          # index end of the larger section
    if (n_up > n_dn) big_side else if (big_side == "h") "l" else "h"
  }
  u1 <- up_end(bx_raw, n_up1, n_dn1)
  u2 <- up_end(by_raw, n_up2, n_dn2)

  # independent cross-check: hottest corner should fall in the largest-overlap quadrant
  k  <- max(2L, floor(0.15 * n)); lo <- 1:k; hi <- (n - k + 1):n
  blk <- function(ri, ci) mean(M[ri, ci], na.rm = TRUE)
  cc  <- c(ll = blk(lo, lo), lh = blk(lo, hi), hl = blk(hi, lo), hh = blk(hi, hi))
  sz  <- c(uu = length(quad$uu), dd = length(quad$dd), ud = length(quad$ud), du = length(quad$du))
  cb  <- list(u1 = "h", u2 = "h", score = -Inf)
  for (a in c("l", "h")) for (b in c("l", "h")) {
    qof <- function(x, y) paste0(if (x == a) "u" else "d", if (y == b) "u" else "d")
    asg <- c(ll = qof("l", "l"), lh = qof("l", "h"), hl = qof("h", "l"), hh = qof("h", "h"))
    s   <- suppressWarnings(cor(cc, sz[asg]))
    if (is.finite(s) && s > cb$score) cb <- list(u1 = a, u2 = b, score = s)
  }
  if (is.na(u1)) u1 <- cb$u1                        # fall back to corner estimate if seam ambiguous
  if (is.na(u2)) u2 <- cb$u2

  if (identical(u1, "l")) M <- M[n:1, , drop = FALSE]   # up in PHF1+     -> right (high x)
  if (identical(u2, "l")) M <- M[, n:1, drop = FALSE]   # up near tangles -> top   (high y)

  # verification: fraction of each axis that is the UP stratum, observed vs expected from n genes
  up_frac <- function(b_raw, u) {
    if (is.na(b_raw)) return(NA_real_)
    lo_sz <- b_raw - 1; hi_sz <- n - b_raw
    (if (identical(u, "h")) hi_sz else lo_sz) / (lo_sz + hi_sz)
  }

  df <- data.frame(x = as.vector(row(M)), y = as.vector(col(M)), val = as.vector(M))
  p <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster(interpolate = TRUE) +               # bilinear-smoothed surface: no blocky squares, no seams
    scale_fill_gradientn(
      colours = c("#08306B", "#2171B5", "#41B6C4", "#FFFFBF", "#FDAE61", "#D73027", "#7F0000"),
      limits = c(hm_min, hm_max), oob = scales::squish,
      na.value = "grey80",                          # the up/down stratum divide, kept light + subtle
      name = expression(-log[10] ~ italic(p))) +
    coord_equal(expand = FALSE) + labs(x = lab_x, y = lab_y) +
    theme_classic(base_size = 8) + fig_theme +   # house standard
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.margin = margin(3, 4, 3, 3), legend.position = "right",
          legend.key.width = grid::unit(5, "pt"), legend.key.height = grid::unit(14, "pt"))
  attr(p, "orient") <- list(
    u1 = u1, u2 = u2, corner = cb, corners = cc, sizes = sz,
    counts = c(up_PHF1 = n_up1, down_PHF1 = n_dn1, up_near = n_up2, down_near = n_dn2),
    obs_up_frac = c(PHF1 = up_frac(bx_raw, u1), near = up_frac(by_raw, u2)),
    exp_up_frac = c(PHF1 = n_up1 / (n_up1 + n_dn1), near = n_up2 / (n_up2 + n_dn2)))
  p
}

# ---- per-gene t-statistic quadrant scatter (companion to the RRHO2 rank map) --------
# Places EACH gene by its two moderated t-statistics: x = cell-autonomous t (+ up in PHF1+),
# y = proximity t (+ up near tangles; the distance-model t is negated). Same sign convention and
# orientation as the heatmap, so top-right = up-up (concordant core). Unlike the rank-overlap map,
# this shows the actual per-gene effects and makes the four sign quadrants concrete. Triple output.
make_quadrant_tstat <- function(j, ct, n_label = 3L) {
  d <- j
  d$quadrant <- dplyr::case_when(
    d$ca_t > 0 & d$prox_t > 0 ~ "up PHF1+ / up near",
    d$ca_t < 0 & d$prox_t < 0 ~ "down PHF1+ / down near",
    d$ca_t > 0 & d$prox_t < 0 ~ "up PHF1+ / down near",
    d$ca_t < 0 & d$prox_t > 0 ~ "down PHF1+ / up near",
    TRUE ~ "on axis")
  d$concordance <- dplyr::case_when(
    d$quadrant %in% c("up PHF1+ / up near", "down PHF1+ / down near") ~ "concordant",
    d$quadrant == "on axis" ~ "on axis",
    TRUE ~ "discordant")
  lim <- max(abs(c(d$ca_t, d$prox_t)), na.rm = TRUE) * 1.05
  pal <- c(concordant = "#0072B2", discordant = "#D55E00", `on axis` = "grey70")   # Okabe-Ito, CB-safe

  d$r <- sqrt(d$ca_t^2 + d$prox_t^2)                          # distance from origin -> most extreme genes
  lab <- d %>% dplyr::filter(concordance != "on axis") %>%
    dplyr::group_by(quadrant) %>% dplyr::slice_max(r, n = n_label, with_ties = FALSE) %>% dplyr::ungroup()

  cnt <- d %>% dplyr::filter(quadrant != "on axis") %>% dplyr::count(quadrant)
  corner <- tibble::tibble(
    quadrant = c("up PHF1+ / up near", "down PHF1+ / up near", "up PHF1+ / down near", "down PHF1+ / down near"),
    x  = c( lim, -lim,  lim, -lim) * 0.98,
    y  = c( lim,  lim, -lim, -lim) * 0.98,
    hj = c(1, 0, 1, 0), vj = c(1, 1, 0, 0)) %>%
    dplyr::left_join(cnt, by = "quadrant") %>% dplyr::mutate(n = ifelse(is.na(n), 0L, n))

  p <- ggplot(d, aes(ca_t, prox_t, colour = concordance)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
    geom_text(data = corner, aes(x, y, label = paste0("n=", n), hjust = hj, vjust = vj),
              inherit.aes = FALSE, size = 2, colour = "grey30") +
    scale_colour_manual(values = pal, breaks = c("concordant", "discordant"), name = NULL) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = expression("Cell-autonomous " * italic(t) * "  (+ = up in PHF1+ vs PHF1-)"),
         y = expression("Proximity " * italic(t) * "  (+ = up near tangles)")) +
    theme_classic(base_size = 8) + fig_theme +   # house standard
    theme(legend.position = "top", legend.key.size = grid::unit(6, "pt"))
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab))
    p <- p + ggrepel::geom_text_repel(
      data = lab, aes(ca_t, prox_t, label = gene), inherit.aes = FALSE,
      size = 2, colour = "grey15", max.overlaps = Inf,
      min.segment.length = 0, segment.size = 0.2, box.padding = 0.2)
  attr(p, "data") <- d
  p
}

# =============================================================================
# RUN
# =============================================================================
for (ct in SUBTYPES) {
  message("== ", ct, " ==")
  out_dir <- file.path(DECOMP_DIR, ct, "rrho2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  j <- load_ranks(ct)
  message("  common genes: ", nrow(j))
  list1 <- data.frame(Genes = j$gene, DDE = j$ca_metric,   stringsAsFactors = FALSE)  # cell-autonomous
  list2 <- data.frame(Genes = j$gene, DDE = j$prox_metric, stringsAsFactors = FALSE)  # proximity

  obj <- RRHO2::RRHO2_initialize(
    list1, list2,
    labels = c("Cell-autonomous (PHF1+ vs -)", "Proximity (up near tangles)"),
    log10.ind = TRUE)

  # log the object structure once (RRHO2 slot names vary between versions; aids debugging)
  writeLines(capture.output(str(obj, max.level = 2, list.len = 20)),
             file.path(out_dir, paste0("rrho2_object_str_", ct, ".txt")))

  # ---- quadrant overlap gene sets (also needed to auto-orient the heatmap) ----
  quad <- list(uu = get_quadrant(obj, "uu"), dd = get_quadrant(obj, "dd"),
               ud = get_quadrant(obj, "ud"), du = get_quadrant(obj, "du"))
  for (q in names(quad))
    writeLines(quad[[q]], file.path(out_dir, paste0("rrho2_genes_", q, "_", ct, ".txt")))

  # ---- headline heatmap (compact, gap-free, deterministically oriented down->up on both axes) ----
  n_up1 <- sum(j$ca_metric > 0);   n_dn1 <- sum(j$ca_metric < 0)     # cell-autonomous up/down
  n_up2 <- sum(j$prox_metric > 0); n_dn2 <- sum(j$prox_metric < 0)   # proximity up/down (near tangles)
  p_hm <- make_rrho_raster(obj, quad, hm_min = HM_MIN, hm_max = HM_MAX,
                           n_up1 = n_up1, n_dn1 = n_dn1, n_up2 = n_up2, n_dn2 = n_dn2)
  ggsave(file.path(out_dir, paste0("plot_rrho2_heatmap_", ct, ".pdf")),
         p_hm, width = 6.5, height = 5.5, units = "cm", device = "pdf")
  ori <- attr(p_hm, "orient")

  # ---- per-gene t-statistic quadrant scatter (companion to the rank map) ----
  p_tq <- make_quadrant_tstat(j, ct)
  ggsave(file.path(out_dir, paste0("plot_rrho2_tstat_quadrant_", ct, ".pdf")),
         p_tq, width = 7, height = 7.5, units = "cm", device = "pdf")
  tq <- attr(p_tq, "data")
  write.table(tq[, c("gene", "ca_t", "prox_t", "quadrant", "concordance")],
              file.path(out_dir, paste0("source_data_rrho2_tstat_quadrant_", ct, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  tq_counts <- table(tq$quadrant)
  tq_spearman <- suppressWarnings(cor(tq$ca_t, tq$prox_t, method = "spearman"))

  # ---- clean area-proportional Euler per quadrant (own euler2, not RRHO2_vennDiagram) ----
  # RRHO2_vennDiagram produces overlapping/garbled labels and needs VennDiagram; euler2 is
  # dependency-free, area-proportional, and paper-styled. n1/n2 include the overlap.
  venn_labels <- list(uu = c("PHF1+ up",   "up near tangles"),
                      dd = c("PHF1+ down", "down near tangles"),
                      ud = c("PHF1+ up",   "down near tangles"),
                      du = c("PHF1+ down", "up near tangles"))
  for (q in names(venn_labels)) {
    cnt <- quad_counts(obj, q)
    if (cnt["n1"] == 0 && cnt["n2"] == 0) next
    p_v <- euler2(cnt["n1"], cnt["n2"], cnt["no"],
                  label_a = venn_labels[[q]][1], label_b = venn_labels[[q]][2])
    ggsave(file.path(out_dir, paste0("plot_rrho2_venn_", q, "_", ct, ".pdf")),
           p_v, width = 6.5, height = 5.5, units = "cm", device = "pdf")
    write.table(data.frame(quadrant = q,
                  set = c(venn_labels[[q]][1], "overlap", venn_labels[[q]][2]),
                  n   = c(cnt["n1"] - cnt["no"], cnt["no"], cnt["n2"] - cnt["no"])),
                file.path(out_dir, paste0("source_data_rrho2_venn_", q, "_", ct, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # ---- GO per quadrant (uu/dd = concordant/shared; ud/du = opposing) ----
  bg <- j$gene
  go_counts <- setNames(integer(4L), c("uu", "dd", "ud", "du"))
  for (q in names(quad)) {
    g <- run_go(quad[[q]], bg)
    go_counts[q] <- nrow(g)
    if (nrow(g)) { g$quadrant <- q
      write.table(g, file.path(out_dir, paste0("go_rrho2_", q, "_", ct, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE) }
  }

  # ---- ORA on the TOP-N strongest genes of the ud direction (up PHF1+ / down near) ----
  # Rank ud-quadrant genes by combined discordant strength (ca_t + (-prox_t) = ca_t - prox_t;
  # both large => strongly up in PHF1+ AND down near tangles) and test the top TOP_DISCORDANT.
  ud_ranked <- j %>% dplyr::filter(ca_t > 0, prox_t < 0) %>%
    dplyr::mutate(disc_score = ca_t - prox_t) %>% dplyr::arrange(dplyr::desc(disc_score))
  ud_top <- head(ud_ranked$gene, TOP_DISCORDANT)
  writeLines(ud_top, file.path(out_dir, paste0("rrho2_genes_ud_top", TOP_DISCORDANT, "_", ct, ".txt")))
  g_udtop <- run_go(ud_top, bg)
  n_udtop_go <- nrow(g_udtop)
  if (n_udtop_go) { g_udtop$quadrant <- paste0("ud_top", TOP_DISCORDANT)
    write.table(g_udtop, file.path(out_dir, paste0("go_rrho2_ud_top", TOP_DISCORDANT, "_", ct, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE) }
  message("  ud top-", TOP_DISCORDANT, " ORA: ", length(ud_top), " genes -> ", n_udtop_go, " GO terms")

  # ---- cross-validation vs the threshold-based decomposition classes ----
  jl_path <- file.path(DECOMP_DIR, ct, paste0("joined_gene_level_", ct, ".csv"))
  if (file.exists(jl_path)) {
    cl <- read.csv(jl_path)[, c("gene", "class")]
    qmap <- tibble(gene = unlist(quad, use.names = FALSE),
                   rrho2 = rep(names(quad), lengths(quad)))
    xt <- cl %>% left_join(qmap, by = "gene") %>% mutate(rrho2 = ifelse(is.na(rrho2), "none", rrho2))
    ctab <- as.data.frame.matrix(table(xt$class, xt$rrho2))
    ctab <- tibble::rownames_to_column(ctab, "class")
    write.table(ctab, file.path(out_dir, paste0("rrho2_vs_class_crosstab_", ct, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  # ---- source data + stats ----
  write.table(j, file.path(out_dir, paste0("source_data_rrho2_ranks_", ct, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c(
    paste0("RRHO2 decomposition -- ", ct),
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("n common genes: ", nrow(j)),
    "Ranking metric = sign(effect) * -log10(pval).",
    "list1 = cell-autonomous (+ up in PHF1+); list2 = proximity (+ up near tangles).",
    "Quadrants: uu (up both) = shared/convergent program [primary]; dd (down both) concordant;",
    "ud (up PHF1+/down near) and du (down PHF1+/up near) = opposing.",
    paste0("quadrant overlap sizes: uu=", length(quad$uu), " dd=", length(quad$dd),
           " ud=", length(quad$ud), " du=", length(quad$du)),
    "heatmap orientation (deterministic: larger gene stratum = larger axis section):",
    paste0("  up-end index -- cell-auto=", ori$u1, "  proximity=", ori$u2,
           "  [corner cross-check: cell-auto=", ori$corner$u1, " proximity=", ori$corner$u2,
           ", r=", round(ori$corner$score, 2), "]"),
    paste0("  up-section fraction observed vs expected -- PHF1+: ",
           sprintf("%.2f vs %.2f", ori$obs_up_frac[["PHF1"]], ori$exp_up_frac[["PHF1"]]),
           " ; near tangles: ",
           sprintf("%.2f vs %.2f", ori$obs_up_frac[["near"]], ori$exp_up_frac[["near"]])),
    paste0("  gene counts -- up-PHF1=", n_up1, " down-PHF1=", n_dn1,
           " ; up-near=", n_up2, " down-near=", n_dn2),
    paste0("  corner means ll/lh/hl/hh = ", paste(sprintf("%.1f", ori$corners), collapse = "/")),
    paste0("GO terms by quadrant: ", paste(sprintf("%s=%d", names(go_counts), go_counts), collapse = " ")),
    "",
    "t-statistic quadrant scatter (per gene; x = cell-auto t, y = proximity t with distance t negated):",
    paste0("  Spearman r(ca_t, prox_t) = ", round(tq_spearman, 3), " (n = ", nrow(j), " genes)"),
    paste0("  genes per quadrant: ",
           paste(sprintf("%s=%d", names(tq_counts), as.integer(tq_counts)), collapse = "; ")),
    "",
    paste0("ud discordant-direction ORA (up PHF1+/down near): top ", TOP_DISCORDANT, " of ",
           nrow(ud_ranked), " ud genes, ranked by ca_t - prox_t -> ", n_udtop_go, " GO terms",
           " (vs ", go_counts[["ud"]], " for the full ud set).")),
    file.path(out_dir, paste0("stats_rrho2_", ct, ".txt")))
}

sink(file.path(DECOMP_DIR, "rrho2_sessionInfo.txt"))
cat("rrho2_decomposition_phf1.r\nRun:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
message("Done. Outputs under ", DECOMP_DIR, "/{ct}/rrho2/")
