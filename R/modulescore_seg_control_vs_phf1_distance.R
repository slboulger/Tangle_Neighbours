#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# modulescore_seg_control_vs_phf1_distance.R
#
# Shared utility - sourced by the scripts above, no panel of its own
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# modulescore_seg_control_vs_phf1_distance.R
#
# SUPPLEMENTARY NEGATIVE CONTROL for the module-score-over-distance analyses.
#
# Sibling of plot_modulescore_vs_phf1_distance_modelp.R (glia: Mancuso Micro,
# Cameron Astro, Pandey Oligo) and plot_phf1_module_vs_phf1_distance_modelp.R
# (neurons: PHF1 signature, Otero-Garcia Exc1/2 UP/DOWN). Same helpers, same
# covariates, same cell-inclusion rules, same model, same plot style, same output
# triple - the only thing that changes is the gene set being scored.
#
# THE CONTROL  stably expressed genes (SEG) by the scSEGIndex method of
#   Lin et al. 2019 GigaScience 8(9):giz106.
#
#   DEFAULT (--seg_qs seg_control/Reference_SEG_panel.qs): PROJECT-SPECIFIC sets, one
#   per CosMx celltype, fitted on our own snRNA-seq reference by
#   R/prepare_seg_geneset_reference.R. Each celltype is scored with the SEG set derived
#   from that same celltype, so the control is matched to the tissue AND the cell type.
#
#   ALTERNATIVE: seg_control/Lin_SEG_2019.qs, the published human list shipped in
#   scMerge's `segList`, extracted by R/prepare_seg_geneset.R. A single shared set; it
#   is accepted here and reused for every celltype.
#
#   A set of genes that is by construction NOT expected to track pathology should
#   show no gradient with distance to the nearest PHF1+ neuron. A flat SEG trace
#   where the reported signatures are sloped shows the reported gradients are not
#   an artefact of scoring an arbitrary gene set on these cells.
#
# CELLTYPES  the four panels the reported gradients are drawn from:
#   Exc-IT-L2-3-CBLN2-HOPX, Exc-IT-L3-5-CHGA-IL1RAPL2 (PHF1 sig / Otero UP), Micro
#   (Mancuso states), Astro (Cameron pooled + subclusters).
#
# WHAT IS PLOTTED  the RAW (model-free) rolling mean, as in the sibling scripts' raw panels:
#   a sliding-window mean of the per-cell AddModuleScore over distance, window width
#   75 um (--window_um), on a 200-point grid, with a 95% CI ribbon (mean +/- 1.96*SEM).
#   NOT z-scored - native AddModuleScore scale, matching make_glia_raw_plot().
#   The fitted log-LMM curve is deliberately NOT plotted.
#
# MODEL  the log-distance LMM is still fitted, because it supplies the statistics and
#   the line styling (solid+thick = model FDR < 0.05, dashed+thin = n.s.) - exactly the
#   relationship the sibling scripts have between their raw panels and their log analysis.
#
#     score ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)
#
#   Distance transform is the canonical one: plain natural log(), divided by the
#   within-celltype SD, NOT centred. Never log1p (distances are strictly > 0 for
#   every modelled PHF1-negative cell; the script stop()s otherwise).
#   Only the log-distance LMM is fitted here.
#
# SCOPE  the whole panel-intersected SEG set, scored as one module. No null draws,
#   no equivalence/TOST machinery, no reported-module overlay, no intensity-adjusted
#   arm - the reported signatures stay in their own figures with their own colours.
#
# OUTPUT under --output_dir, group = "seg_control":
#   plot_modulescore_raw_dist_seg_control_<CT>.pdf
#   plot_modulescore_raw_dist_seg_control_composite.pdf          (4-panel, A4 width)
#   source_data_modulescore_raw_dist_seg_control_<CT>_rollmean.tsv   (the drawn rows)
#   stats_modulescore_raw_dist_seg_control_<CT>.txt   (rolling-mean spec + LMM stats)
#   stats_modulescore_seg_control_overlap.tsv, ..._members.tsv
#   stats_modulescore_lmm_coeffs.tsv
#
# RUN (interactive session, e.g. OpenOnDemand; conda activate dgenv):
#   Rscript R/prepare_seg_geneset.R                       # once, needs scMerge
#   Rscript R/modulescore_seg_control_vs_phf1_distance.R  # defaults = 1000 um cap

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(argparse)
})

# Project root: HPC path first (the project convention), local Mac mount as fallback,
# so this can be run either on the HPC or locally on <RDS_ROOT> without editing.
PROJECT_ROOTS <- c(
  "<PROJECT_ROOT>/phf1_v2",
  "<PROJECT_ROOT>/phf1_v2"
)
.root <- PROJECT_ROOTS[dir.exists(PROJECT_ROOTS)]
if (length(.root) == 0)
  stop("No project root found. Tried:\n  ", paste(PROJECT_ROOTS, collapse = "\n  "), call. = FALSE)
setwd(.root[1])
cat("Project root:", getwd(), "\n")

source("R/palettes.R")  # fig_theme, neuron_order, glia_order, ...

set.seed(42)

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--seu",
  default = "seu_PHF1.rds",
  help = "Path to seu_PHF1.rds [default: seu_PHF1.rds]")
parser$add_argument("--seg_qs", default = "seg_control/Reference_SEG_panel.qs",
  help = paste("SEG geneset .qs, relative to the project root. Either one element per",
               "celltype (R/prepare_seg_geneset_reference.R) or a single shared set",
               "(R/prepare_seg_geneset.R). [default: seg_control/Reference_SEG_panel.qs]"))
parser$add_argument("--output_dir",
  default = "plots/modulescore_seg_control_1000um",
  help = "Output directory")
parser$add_argument("--seed", type = "integer", default = 42,
  help = "Seed (AddModuleScore + provenance) [default: 42]")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "If > 0, restrict modelled cells to dist_to_phf1_um <= this (um); 0 = no cap")
parser$add_argument("--window_um", type = "double", default = 75,
  help = "Rolling-mean window WIDTH (um) for the plotted curve [default: 75]")
parser$add_argument("--neuron_coeffs_tsv",
  default = "plots/phf1_module_vs_phf1_distance_modelp_1000um/stats_modulescore_lmm_coeffs.tsv",
  help = "Reported neuron-signature coefficients, for the context block in the stats logs")
parser$add_argument("--glia_coeffs_tsv",
  default = "plots/module_score_vs_phf1_distance_modelp_1000um/stats_modulescore_lmm_coeffs.tsv",
  help = "Reported glia-module coefficients, for the context block in the stats logs")
args <- parser$parse_args()

MAX_DIST <- if (args$max_dist_um > 0) args$max_dist_um else Inf
if (is.finite(MAX_DIST))
  cat(sprintf("Distance cap: modelling cells with dist_to_phf1_um <= %g um\n", MAX_DIST))

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Constants (identical to the sibling scripts) ---------------------------------
MIN_GENES_PRESENT <- 10
CTRL              <- 100
NBIN              <- 24
SIG_ALPHA         <- 0.05
HALF_WINDOW       <- args$window_um / 2   # rolling-mean half-window (um), as the siblings

GROUP_TAG   <- "seg_control"   # filename + lmm_coeffs `group` token
SCALE_LABEL <- "log_um"        # lmm_coeffs `distance_scale`; log-distance LMM only

# The four panels the reported gradients come from. Plot/composite order follows
# the canonical c(neuron_order, glia_order) from R/palettes.R.
CELLTYPES  <- c("Exc-IT-L2-3-CBLN2-HOPX", "Exc-IT-L3-5-CHGA-IL1RAPL2", "Micro", "Astro")
CELLTYPES  <- CELLTYPES[order(match(CELLTYPES, c(neuron_order, glia_order)))]

# Neutral grey: the control must read as a control, not as an additional signature.
# Deliberately NOT from MODULE_COLOURS / MERGED_COLOURS - those stay pinned to the
# reported modules and are not touched here.
SEG_COLOUR <- "#666666"

XLAB <- expression("Distance to PHF1+ neuron (" * mu * "m)")  # mu via plotmath (survives default pdf device)

MODEL_LABEL <- "LMM (log distance)"
FORMULA_STR <- "score ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)"

##  ............................................................................
##  Helpers - copied VERBATIM from the sibling scripts (see header)         ####
# The sibling scripts run argparse + readRDS at top level and so cannot be source()d,
# so the helper block is copied. Only `extract_celltype_data` differs from the glia
# script, and only by carrying `phf1_pos` (as in the neuron script) so the
# PHF1-negative filter can be applied explicitly.

# First present column name from a set of candidates (defensive against naming).
pick_col <- function(meta, candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0)
    stop(sprintf("None of the %s columns (%s) found in seu meta.data",
                 what, paste(candidates, collapse = ", ")))
  hit[1]
}

# Intersect each set with the panel; keep those with >= min_genes; build overlap table.
prep_sets <- function(sets, panel, min_genes, collection, group) {
  present <- lapply(sets, function(g) intersect(g, panel))
  n_total <- vapply(sets, length, integer(1))
  n_pres  <- vapply(present, length, integer(1))
  scored  <- n_pres >= min_genes
  overlap <- tibble(
    collection   = collection,
    group        = group,
    set          = names(sets),
    n_total      = as.integer(n_total),
    n_present    = as.integer(n_pres),
    frac_present = round(n_pres / n_total, 4),
    scored       = scored
  )
  list(present = present[scored], overlap = overlap)
}

# AddModuleScore for one prefix; returns data.frame of scores renamed to set names.
# Retries with a smaller nbin if the control-gene binning fails.
score_group <- function(seu_ct, present_sets, prefix, seed) {
  set.seed(seed)
  seu_ct <- tryCatch(
    AddModuleScore(seu_ct, features = present_sets, name = prefix,
                   assay = "SCT", ctrl = CTRL, nbin = NBIN, seed = seed),
    error = function(e) {
      cat("    AddModuleScore failed at nbin=", NBIN, " (", conditionMessage(e),
          ") - retrying at nbin=15\n", sep = "")
      set.seed(seed)
      AddModuleScore(seu_ct, features = present_sets, name = prefix,
                     assay = "SCT", ctrl = CTRL, nbin = 15, seed = seed)
    })
  cols <- paste0(prefix, seq_along(present_sets))     # AddModuleScore appends list index
  stopifnot(
    "AddModuleScore column count != number of scored sets" =
      length(cols) == length(present_sets) && all(cols %in% colnames(seu_ct@meta.data))
  )
  sc <- seu_ct@meta.data[, cols, drop = FALSE]
  colnames(sc) <- names(present_sets)                 # load-bearing index->name mapping
  sc
}

# Subset Seurat to a celltype, score the supplied feature groups, return one
# data.frame with cell_id + PHF1 status + covariates + all module-score columns.
extract_celltype_data <- function(seu, ct, feature_groups, seed,
                                  sex_col, age_col, pmi_col) {
  seu_ct <- subset(seu, celltype == ct)
  DefaultAssay(seu_ct) <- "SCT"
  md <- seu_ct@meta.data
  cat(sprintf("  %s cells: %d\n", ct, nrow(md)))

  id_cands <- c("cell_ID", "cell_id", "CellID")
  id_col   <- id_cands[id_cands %in% colnames(md)]
  cell_id  <- if (length(id_col) > 0) as.character(md[[id_col[1]]]) else rownames(md)

  base <- data.frame(
    cell_id         = cell_id,
    phf1_pos        = as.logical(as.character(md$PHF1)),   # boolean post-stain
    dist_to_phf1_um = as.numeric(md$dist_to_phf1_um),
    percent_neg     = as.numeric(md$percent.neg),
    nUMI_log        = log2(as.numeric(md$nCount_RNA) + 1),
    sample_id       = as.character(md$sample_id),
    Sex             = as.character(md[[sex_col]]),
    Age             = as.numeric(md[[age_col]]),
    PMI             = as.numeric(md[[pmi_col]]),
    stringsAsFactors = FALSE
  )

  scores <- do.call(cbind, lapply(names(feature_groups), function(prefix)
    score_group(seu_ct, feature_groups[[prefix]], prefix, seed)))

  cbind(base, scores)   # row-aligned: all derived from seu_ct@meta.data order
}

# Complete-covariate filter + factor levels (verbatim from the sibling scripts).
prep_celltype_df <- function(df) {
  df <- df[!is.na(df$dist_to_phf1_um) & !is.na(df$percent_neg) &
             !is.na(df$nUMI_log) & !is.na(df$Age) & !is.na(df$PMI) & !is.na(df$Sex), ]
  df$Sex       <- droplevels(factor(df$Sex))
  df$sample_id <- factor(df$sample_id)
  df
}

# Observed LMM for one module. Returns fit stats + the fitted model.
fit_obs_lmer <- function(sub) {
  m <- tryCatch(
    lmerTest::lmer(score ~ dist_scaled + nUMI_log + percent_neg + Sex + Age_s + PMI_s +
                     (1 | sample_id), data = sub, REML = TRUE),
    error = function(e) { message("    obs lmer failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  cf <- coef(summary(m))
  if (!"dist_scaled" %in% rownames(cf)) return(NULL)
  list(model = m,
       estimate = cf["dist_scaled", "Estimate"],
       se       = cf["dist_scaled", "Std. Error"],
       t        = cf["dist_scaled", "t value"],
       df       = cf["dist_scaled", "df"],
       p        = cf["dist_scaled", "Pr(>|t|)"],
       singular = isSingular(m))
}

# (The sibling scripts' predict_grid() / predict_grid_gam() are NOT copied here:
#  the fitted curve is not plotted, so there is no grid to predict onto. The LMM is
#  used only for its coefficient, CI and p.)

# Model-free rolling (sliding-window) mean of y over x, evaluated on a grid.
# Returns mean, +/- 1 SEM, and window cell count at each grid point (95% CI = +/- 1.96 SEM).
# Verbatim from plot_modulescore_vs_phf1_distance_modelp.R.
roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}

# Fix the panel to a constant physical width so x-axes align across ALL plots
# regardless of legend size; total figure width adapts to the legend.
# Verbatim from plot_phf1_module_vs_phf1_distance_modelp.R.
save_fixed_panel <- function(p, path) {
  g <- ggplot2::ggplotGrob(p)
  pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
  g$widths[pcol] <- grid::unit(1.6, "in")   # wide enough for the x-axis title (avoids clip)
  ggsave(path, g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
         height = 2.1, units = "in", device = "pdf")
}

# Two significant figures, for the stats logs.
sf2 <- function(x) if (length(x) == 0 || all(is.na(x))) "NA" else formatC(signif(x, 2), format = "g")

##  ............................................................................
##  Load the SEG gene set                                                   ####
seg_path <- args$seg_qs
if (!file.exists(seg_path))
  stop(paste0("SEG geneset not found:\n  ", seg_path,
              "\nRun R/prepare_seg_geneset.R first (it needs scMerge installed)."),
       call. = FALSE)

cat("Loading SEG gene sets:", seg_path, "\n")
seg <- qs::qread(seg_path)
stopifnot("SEG geneset must be a named list of character vectors" =
            is.list(seg) && !is.null(names(seg)) &&
            all(vapply(seg, is.character, logical(1))))
cat("  sets:", paste(names(seg), collapse = ", "),
    "| sizes:", paste(vapply(seg, length, integer(1)), collapse = ", "), "\n")

# The sets are PER CELLTYPE: each CosMx celltype is scored with the SEG set derived
# from that same celltype in the snRNA-seq reference. A single shared set (the
# Lin_SEG_2019.qs shape, one element) is also accepted and reused for every celltype.
SEG_PER_CELLTYPE <- all(CELLTYPES %in% names(seg))
if (SEG_PER_CELLTYPE) {
  cat("  -> per-celltype SEG sets detected; each celltype uses its own set\n")
} else if (length(seg) == 1L) {
  cat("  -> single shared SEG set; the same genes are used for every celltype\n")
  seg <- setNames(rep(seg, length(CELLTYPES)), CELLTYPES)
} else {
  stop(paste0("SEG geneset must either contain one element per requested celltype\n",
              "  (", paste(CELLTYPES, collapse = ", "), ")\n",
              "  or be a single shared set. Found: ",
              paste(names(seg), collapse = ", ")), call. = FALSE)
}

##  ............................................................................
##  Load Seurat object, panel intersection                                  ####
cat("Loading Seurat object (this is large)...\n")
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "SCT"

stopifnot("celltype column missing"        = "celltype"        %in% colnames(seu@meta.data),
          "PHF1 column missing"            = "PHF1"            %in% colnames(seu@meta.data),
          "dist_to_phf1_um column missing" = "dist_to_phf1_um" %in% colnames(seu@meta.data),
          "percent.neg column missing"     = "percent.neg"     %in% colnames(seu@meta.data),
          "nCount_RNA column missing"      = "nCount_RNA"      %in% colnames(seu@meta.data),
          "sample_id column missing"       = "sample_id"       %in% colnames(seu@meta.data))
stopifnot("requested celltypes missing from seu$celltype" =
            all(CELLTYPES %in% unique(as.character(seu@meta.data$celltype))))

sex_col <- pick_col(seu@meta.data, c("Sex"), "Sex")
age_col <- pick_col(seu@meta.data, c("Age"), "Age")
pmi_col <- pick_col(seu@meta.data, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")
cat(sprintf("Donor covariate columns: Sex=%s Age=%s PMI=%s\n", sex_col, age_col, pmi_col))

panel_genes <- rownames(seu[["SCT"]])
cat("CosMx panel genes (SCT assay):", length(panel_genes), "\n")

# Panel-intersect each celltype's own set. `prep_sets` is applied once per celltype so
# the overlap table has one auditable row per celltype, in the sibling scripts' schema.
SEG_PRESENT <- list(); SEG_N <- integer(0)
overlap_rows <- list()
for (ct in CELLTYPES) {
  pr <- prep_sets(setNames(list(seg[[ct]]), ct), panel_genes, MIN_GENES_PRESENT,
                  "Reference_scSEGIndex", GROUP_TAG)
  overlap_rows[[ct]] <- pr$overlap
  if (length(pr$present) == 0)
    stop(sprintf(paste0("%s: only %d of %d SEG genes are on the panel, below ",
                        "MIN_GENES_PRESENT=%d; the control cannot be scored."),
                 ct, pr$overlap$n_present[1], pr$overlap$n_total[1], MIN_GENES_PRESENT),
         call. = FALSE)
  SEG_PRESENT[[ct]] <- pr$present[[1]]
  SEG_N[ct] <- length(SEG_PRESENT[[ct]])
}
overlap_tab <- bind_rows(overlap_rows)
write.table(overlap_tab,
            file.path(args$output_dir, "stats_modulescore_seg_control_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("SEG / panel overlap (per celltype):\n"); print(as.data.frame(overlap_tab))

# One constant legend label across panels, so the composite's collected legend stays a
# single key; the per-celltype n lives in the stats logs and the source data.
SEG_LABEL <- if (SEG_PER_CELLTYPE) "Reference SEG (per celltype)" else
  sprintf("SEG (n = %d)", SEG_N[[1]])
cat("Panel-intersected SEG sizes:",
    paste(sprintf("%s=%d", names(SEG_N), SEG_N), collapse = "  "), "\n")

# Every gene used, with on_panel = actually scored (same pattern as the sibling scripts).
members_df <- function(sets, group) bind_rows(lapply(names(sets), function(m)
  tibble(group = group, module = m, gene = sets[[m]], on_panel = sets[[m]] %in% panel_genes)))
seg_members <- bind_rows(lapply(CELLTYPES, function(ct)
  tibble(group = GROUP_TAG, module = ct, gene = seg[[ct]],
         on_panel = seg[[ct]] %in% panel_genes)))
write.table(seg_members,
            file.path(args$output_dir, "stats_modulescore_seg_control_members.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved SEG members (on_panel = scored):",
    sum(seg_members$on_panel), "of", nrow(seg_members), "genes on panel\n")

##  ............................................................................
##  Score + build model frames, per celltype                                ####
# Module scores are computed on the FULL celltype population (as the siblings), so
# the AddModuleScore expression-bin control background is identical to theirs; the
# PHF1-negative filter and the distance cap are applied afterwards.
CT_DATA  <- list()
CT_NOTES <- list()
for (ct in CELLTYPES) {
  cat(sprintf("\nScoring %s (SEG negative control)...\n", ct))
  d <- extract_celltype_data(
    seu, ct,
    feature_groups = list(SEG_ = setNames(list(SEG_PRESENT[[ct]]), "SEG")),
    seed = args$seed, sex_col = sex_col, age_col = age_col, pmi_col = pmi_col)

  n_all  <- nrow(d)
  n_pos  <- sum(d$phf1_pos, na.rm = TRUE)
  d      <- d[!is.na(d$phf1_pos) & !d$phf1_pos, ]      # PHF1-NEGATIVE only
  n_neg  <- nrow(d)
  d      <- prep_celltype_df(d)
  n_cc   <- nrow(d)
  cat(sprintf("  %s: %d cells, %d PHF1+ dropped, %d PHF1-negative, %d with complete covariates\n",
              ct, n_all, n_pos, n_neg, n_cc))
  CT_DATA[[ct]]  <- d
  CT_NOTES[[ct]] <- list(n_all = n_all, n_phf1_pos = n_pos, n_neg = n_neg, n_complete = n_cc)
}
rm(seu); invisible(gc())

##  ............................................................................
##  Reported-module context (read-only; no new statistics)                   ####
# The reported log-distance slopes for these celltypes, so each stats log states
# what the control is being compared against without cross-referencing files.
read_reported <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- utils::read.delim(path, stringsAsFactors = FALSE)
  if (!all(c("celltype", "module", "distance_scale", "slope_per_sd") %in% names(x))) return(NULL)
  x[x$distance_scale == "log_um" & x$celltype %in% CELLTYPES,
    intersect(c("celltype", "module", "slope_per_sd", "ci_lo_per_sd", "ci_hi_per_sd",
                "p_LMM", "lmm_padj"), names(x))]
}
REPORTED <- bind_rows(read_reported(args$neuron_coeffs_tsv),
                      read_reported(args$glia_coeffs_tsv))
if (is.null(REPORTED) || nrow(REPORTED) == 0) {
  REPORTED <- NULL
  cat("\nNOTE: no reported-module coefficient tables found; the stats logs will omit\n",
      "      the reported-slope context block.\n", sep = "")
} else {
  cat(sprintf("\nReported log-distance slopes loaded for context: %d rows, %d celltypes\n",
              nrow(REPORTED), dplyr::n_distinct(REPORTED$celltype)))
}

##  ............................................................................
##  Fit (log-distance LMM), pass 1                                          ####
# Every celltype is fitted first, so BH can be applied across the four celltypes
# (pass 2) before anything is plotted or written.
FITS <- list()
for (ct in CELLTYPES) {
  cat(sprintf("\n=== %s / %s  (model: %s%s) ===\n", GROUP_TAG, ct, SCALE_LABEL,
              if (is.finite(MAX_DIST)) sprintf(", <=%gum", MAX_DIST) else ""))
  df_full <- CT_DATA[[ct]]

  # Donor covariates scaled on the full (pre-cap) PHF1-negative set, as the siblings.
  df_full$Age_s <- as.numeric(scale(df_full$Age))
  df_full$PMI_s <- as.numeric(scale(df_full$PMI))

  df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full
  if (nrow(df) < 50) { cat("  <50 cells within cap; skipping.\n"); next }

  # Canonical distance transform: plain natural log(), divided by the within-celltype
  # SD, NOT centred. Distances are strictly > 0 for PHF1-negative cells, so there is
  # no offset and no log1p - stop() rather than absorb a bad value.
  if (any(df$dist_to_phf1_um <= 0, na.rm = TRUE))
    stop(sprintf(paste0("%s: %d modelled cell(s) have dist_to_phf1_um <= 0. ",
                        "Distance to the nearest PHF1+ neuron must be strictly > 0 for ",
                        "every PHF1-negative cell; log() is the canonical transform and ",
                        "log1p is not permitted here. Fix upstream in ",
                        "R/label_phf1_neighbours.r rather than adding an offset."),
                 ct, sum(df$dist_to_phf1_um <= 0, na.rm = TRUE)), call. = FALSE)

  dt      <- log(df$dist_to_phf1_um)
  dist_sd <- sd(dt)
  df$dist_scaled <- dt / dist_sd
  x_cap <- if (is.finite(MAX_DIST)) MAX_DIST else as.numeric(quantile(df$dist_to_phf1_um, 0.99))

  sub <- df[, c("cell_id", "dist_scaled", "nUMI_log", "percent_neg",
                "Sex", "Age_s", "PMI_s", "sample_id")]
  sub$score <- df[["SEG"]]
  obs <- fit_obs_lmer(sub)
  if (is.null(obs)) { cat("  LMM failed, skipped\n"); next }

  # PLOTTED CURVE: raw model-free rolling mean on the native AddModuleScore scale
  # (NOT z-scored), 200-point grid from 0 to the cap - as make_glia_raw_plot().
  # The LMM above supplies the statistics and the solid/dashed styling, not the curve.
  grid_um <- seq(0, x_cap, length.out = 200)
  y  <- df[["SEG"]]
  ok <- !is.na(y) & !is.na(df$dist_to_phf1_um)
  roll_df <- data.frame(module = SEG_LABEL,
                        roll_mean(df$dist_to_phf1_um[ok], y[ok], grid_um, HALF_WINDOW),
                        dist_to_phf1_um = grid_um)
  n_empty <- sum(roll_df$n_window == 0)
  if (n_empty > 0)
    cat(sprintf("  note: %d/200 grid points have no cells in the +/-%g um window\n",
                n_empty, HALF_WINDOW))

  # Standardised effect size: slope per SD of log-distance, divided by
  # the total outcome SD implied by the model (donor RE + residual).
  vc       <- as.data.frame(lme4::VarCorr(obs$model))
  sd_tot   <- sqrt(sum(vc$vcov, na.rm = TRUE))
  cohens_d <- obs$estimate / sd_tot

  FITS[[ct]] <- list(
    stats = tibble(
      group = GROUP_TAG, celltype = ct, module = SEG_LABEL, distance_scale = SCALE_LABEL,
      slope_per_sd = obs$estimate, se = obs$se,
      ci_lo_per_sd = obs$estimate - 1.96 * obs$se,
      ci_hi_per_sd = obs$estimate + 1.96 * obs$se,
      slope_per_unit = obs$estimate / dist_sd,
      t = obs$t, df = obs$df, p_LMM = obs$p,
      singular = obs$singular, dist_sd = dist_sd, cohens_d = cohens_d
    ),
    roll = roll_df,
    meta = list(n_pre_cap = nrow(df_full), n_modelled = nrow(df),
                n_donors = dplyr::n_distinct(df$sample_id),
                dist_sd = dist_sd, x_cap = x_cap, n_genes = SEG_N[[ct]],
                n_empty_windows = n_empty,
                vc_txt = paste(sprintf("%s=%.5g", vc$grp, vc$vcov), collapse = ", "))
  )
  cat(sprintf("  fitted: n=%d cells, %d donors, slope_per_sd=%s, p=%s\n",
              nrow(df), dplyr::n_distinct(df$sample_id), sf2(obs$estimate), sf2(obs$p)))
}

if (length(FITS) == 0) stop("No models were fitted; nothing to write.", call. = FALSE)

##  ............................................................................
##  BH adjustment (pass 2)                                                  ####
# BH family, stated in every stats log: the sibling scripts BH-adjust WITHIN
# (group x celltype x model), which here would be a family of ONE module and make
# lmm_padj identical to p_LMM. BH is therefore applied across the CELLTYPES, which is
# the actual multiplicity of this figure.
PADJ <- bind_rows(lapply(FITS, `[[`, "stats"))
PADJ$lmm_padj    <- p.adjust(PADJ$p_LMM, method = "BH")
PADJ$significant <- PADJ$lmm_padj < SIG_ALPHA
PADJ$verdict     <- ifelse(is.na(PADJ$significant), "not testable",
                    ifelse(PADJ$significant,
                           "significant distance effect (model FDR<0.05)",
                           "n.s. (model FDR)"))
cat(sprintf("\nBH across %d celltypes: %d significant\n",
            nrow(PADJ), sum(PADJ$significant, na.rm = TRUE)))

##  ............................................................................
##  Plot + write, per celltype                                              ####
SIG_BASIS     <- "padj<0.05"
lmm_coef_rows <- list()
PANELS        <- list()

for (ct in names(FITS)) {
  fo   <- FITS[[ct]]
  st   <- PADJ[PADJ$celltype == ct, ]
  meta <- fo$meta
  sig  <- isTRUE(st$significant[1])

  ct_tag <- gsub("[^A-Za-z0-9_-]", "_", ct)

  roll_df <- fo$roll %>%
    mutate(module = factor(SEG_LABEL, levels = SEG_LABEL),
           significant = sig,
           sig = factor(ifelse(sig, SIG_BASIS, "ns"), levels = c(SIG_BASIS, "ns")))

  # ---- figure: raw 75 um rolling mean, one neutral-grey control trace ----
  # Ribbon = rolling mean +/- 1.96*SEM (95% CI), as make_glia_raw_plot(). Line style
  # carries the LOG-model FDR, so the panel still says whether the control's distance
  # slope is significant. No plot title (panel label + figure legend do that work).
  p <- ggplot(roll_df, aes(dist_to_phf1_um, roll_mean,
                           colour = module, fill = module, group = module)) +
    geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
                alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = sig, linewidth = sig)) +
    scale_colour_manual(values = setNames(SEG_COLOUR, SEG_LABEL),
                        name = "Negative control", drop = FALSE) +
    scale_fill_manual(values = setNames(SEG_COLOUR, SEG_LABEL), guide = "none",
                      drop = FALSE) +
    # linetype + linewidth share one legend (same name + limits) so BOTH a
    # bold-solid (significant) and a thin-dashed (ns) key always show.
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(SIG_BASIS, "ns")),
                          name = "Log-model FDR", drop = FALSE,
                          limits = c(SIG_BASIS, "ns")) +
    scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(SIG_BASIS, "ns")),
                           name = "Log-model FDR", drop = FALSE,
                           limits = c(SIG_BASIS, "ns")) +
    guides(colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
    labs(x = XLAB, y = "Module score") +
    coord_cartesian(xlim = c(0, meta$x_cap)) +
    theme_classic(base_size = 8) + fig_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.key.width = grid::unit(18, "pt"),
          legend.key.height = grid::unit(8, "pt"),
          legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(t = 4, r = 6, b = 4, l = 5))

  PANELS[[ct]] <- p

  save_fixed_panel(p, file.path(args$output_dir,
    sprintf("plot_modulescore_raw_dist_%s_%s.pdf", GROUP_TAG, ct_tag)))

  # ---- source data: exact drawn rows (rolling mean + SEM + window n) ----
  roll_out <- roll_df %>%
    transmute(group = GROUP_TAG, celltype = ct, module = as.character(module),
              window_um = args$window_um, log_significant = significant,
              dist_to_phf1_um, roll_mean, sem, n_window)
  write.table(roll_out, file.path(args$output_dir,
    sprintf("source_data_modulescore_raw_dist_%s_%s_rollmean.tsv", GROUP_TAG, ct_tag)),
    sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- stats log (triple) ----
  sink(file.path(args$output_dir,
    sprintf("stats_modulescore_raw_dist_%s_%s.txt", GROUP_TAG, ct_tag)))
  cat("SEG negative control: module score vs distance to PHF1+ neuron\n")
  cat("Group:", GROUP_TAG, " | celltype:", ct, " | distance scale:", SCALE_LABEL, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  cat("PLOTTED CURVE (raw, model-free)\n")
  cat("  Rolling (sliding-window) mean of the per-cell AddModuleScore over distance,\n")
  cat("  +/- 95% CI (mean +/- 1.96*SEM). Window width:", args$window_um,
      "um (half-window", HALF_WINDOW, "um);\n")
  cat("  grid: 200 points from 0 to", meta$x_cap, "um.\n")
  cat("  NOT z-scored - native AddModuleScore scale, as the sibling scripts' raw panels.\n")
  cat("  Grid points with an empty window:", meta$n_empty_windows, "of 200.\n")
  cat("  The fitted log-LMM curve is NOT plotted; the LMM below supplies the statistics\n")
  cat("  and the line style (solid+thick = model FDR<0.05, dashed+thin = n.s.).\n")
  cat("  NO PHF1+ reference band on this figure (it is a PHF1-negative control).\n\n")

  cat("NEGATIVE CONTROL GENE SET\n")
  cat("  scSEGIndex stably expressed genes (Lin Y, Ghazanfar S, Strbenac D, et al.\n")
  cat("  Evaluating stably expressed genes in single cells. GigaScience 2019;8(9):giz106).\n")
  cat("  Source file:", args$seg_qs, "\n")
  if (SEG_PER_CELLTYPE) {
    cat("  PROJECT-SPECIFIC set, fitted on our own snRNA-seq reference for THIS celltype\n")
    cat("  by R/prepare_seg_geneset_reference.R (scMerge::scSEGIndex, cell_type = NULL,\n")
    cat("  panel-restricted universe, segIdx > 80th percentile). Per-gene stability\n")
    cat("  features (lambda, sigma, sigma^2, omega, omega*, mu, F_donor) are in\n")
    cat("  seg_control/seg_index_panel_<celltype>.tsv.\n")
    cat("  Genes in this celltype's set:", length(seg[[ct]]),
        " | on CosMx panel (scored):", meta$n_genes, "\n")
  } else {
    cat("  Published human list from scMerge's segList, extracted by\n")
    cat("  R/prepare_seg_geneset.R; the same shared set is used for every celltype.\n")
    cat("  Genes in published set:", length(seg[[ct]]),
        " | on CosMx panel (scored):", meta$n_genes, "\n")
  }
  cat("  Rationale: genes selected for stable expression should NOT track pathology,\n")
  cat("  so a flat trace here shows the reported signature gradients are not an\n")
  cat("  artefact of scoring an arbitrary gene set on these cells.\n\n")

  cat(MODEL_LABEL, ":\n", sep = "")
  cat("  ", FORMULA_STR, "\n", sep = "")
  cat("  Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
      " seed=", args$seed, ").\n", sep = "")
  cat("  Distance transform: plain natural log(), divided by the within-celltype SD,\n")
  cat("  NOT centred (dist_sd =", round(meta$dist_sd, 4), "), computed after the cap.\n")
  cat("  log1p is NOT used: distance is strictly > 0 for every modelled cell.\n")
  cat("  Only the log-distance LMM is fitted.\n")

  cat("\nCELLS\n")
  cat("  PHF1-negative only (PHF1 post-stain boolean), complete covariates.\n")
  cat("  ", ct, ": ", CT_NOTES[[ct]]$n_all, " cells -> ", CT_NOTES[[ct]]$n_phf1_pos,
      " PHF1+ dropped -> ", CT_NOTES[[ct]]$n_neg, " PHF1-negative -> ",
      CT_NOTES[[ct]]$n_complete, " with complete covariates\n", sep = "")
  if (is.finite(MAX_DIST)) {
    cat("  Distance cap: cells with dist_to_phf1_um <=", MAX_DIST, "um.\n")
    cat("  Cells modelled (within cap):", meta$n_modelled, " of ", meta$n_pre_cap,
        " PHF1-negative | donors:", meta$n_donors, "\n")
  } else {
    cat("  Cells (PHF1-negative):", meta$n_modelled, " | donors:", meta$n_donors, "\n")
  }
  cat("  Module scores were computed on the FULL celltype population before the PHF1\n")
  cat("  and distance filters, so the AddModuleScore expression-bin control background\n")
  cat("  matches the reported analyses exactly.\n")
  cat("  Donor random intercept (1|sample_id); n = ", meta$n_donors,
      " donors.\n", sep = "")

  cat("\nRESULT\n")
  print(as.data.frame(st %>% dplyr::select(
    celltype, module, distance_scale, slope_per_sd, se, ci_lo_per_sd, ci_hi_per_sd,
    slope_per_unit, t, df, p_LMM, lmm_padj, significant, singular, verdict)))

  cat("\n=== Effect sizes ===\n")
  cat("  Unstandardised: slope per SD of log-distance = ", sf2(st$slope_per_sd[1]),
      " (95% CI ", sf2(st$ci_lo_per_sd[1]), " to ", sf2(st$ci_hi_per_sd[1]), ")\n", sep = "")
  cat("  Per log-um: ", sf2(st$slope_per_unit[1]), " (dist_sd = ",
      round(meta$dist_sd, 4), ")\n", sep = "")
  cat("  Standardised: Cohen's d = slope / sqrt(donor RE var + residual var) = ",
      sf2(st$cohens_d[1]), "\n", sep = "")
  cat("  Variance components: ", meta$vc_txt, "\n", sep = "")

  if (!is.null(REPORTED)) {
    rep_ct <- REPORTED[REPORTED$celltype == ct, ]
    if (nrow(rep_ct) > 0) {
      cat("\n=== Reported modules for this celltype (log LMM), for context ===\n")
      cat("  Read from the published coefficient tables; no new statistics are computed here.\n")
      print(as.data.frame(rep_ct[order(-abs(rep_ct$slope_per_sd)), ]))
      steep   <- max(abs(rep_ct$slope_per_sd), na.rm = TRUE)
      steep_m <- rep_ct$module[which.max(abs(rep_ct$slope_per_sd))]
      cat("  Steepest reported |slope_per_sd|: ", sf2(steep), " (", steep_m, ")\n", sep = "")
      if (is.finite(st$slope_per_sd[1]))
        cat("  |SEG slope| / |steepest reported slope| = ",
            sf2(abs(st$slope_per_sd[1]) / steep), "\n", sep = "")
    }
  }

  cat("\nSignificance = BH-adjusted model p (lmm_padj); padj <", SIG_ALPHA,
      "-> significant.\n")
  cat("BH family: the sibling scripts BH-adjust WITHIN (group x celltype x model). With a\n")
  cat("  single control module that family has size 1 and lmm_padj would equal p_LMM.\n")
  cat("  BH is therefore applied across the", nrow(PADJ), "\n")
  cat("  CELLTYPES, which is the actual multiplicity of this figure.\n")

  cat("\nNOTES\n")
  cat("  1. AddModuleScore subtracts an expression-bin-matched control gene set, which\n")
  cat("     pins module scores near zero BY CONSTRUCTION. Any comparison with the\n")
  cat("     reported signatures is therefore between SLOPES, not baseline levels.\n")
  cat("  2. The SEG trace bounds the shared, panel-wide component of an apparent\n")
  cat("     distance trend.\n")
  cat("  3. This is a single gene set per celltype, not a resampled null distribution.\n")

  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()

  lmm_coef_rows[[ct]] <- st %>%
    dplyr::select(group, celltype, module, distance_scale, slope_per_sd, se,
                  ci_lo_per_sd, ci_hi_per_sd, slope_per_unit, t, df, p_LMM, lmm_padj,
                  significant, singular)

  cat(sprintf("  Wrote %s / %s (raw rolling mean)\n", GROUP_TAG, ct))
}

##  ............................................................................
##  Written-once table + composite figure                                   ####
write.table(bind_rows(lmm_coef_rows),
            file.path(args$output_dir, "stats_modulescore_lmm_coeffs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote stats_modulescore_lmm_coeffs.tsv\n")

# Four-panel supplementary composite (raw rolling mean), shared legend, A4 page width.
# 2 x 2 rather than 1 x 4: at 19.05 cm four panels in a row leave only
# ~4.8 cm each and the x-axis titles overlap between panels.
# The standalone panels stay title-free (house style: panel label + figure legend do
# that work); the composite copies get a small celltype title because otherwise the
# four panels are visually indistinguishable - every legend here is identical.
panels <- PANELS[CELLTYPES[CELLTYPES %in% names(PANELS)]]
if (length(panels) > 0) {
  titled <- lapply(names(panels), function(ct)
    panels[[ct]] + labs(title = ct) +
      theme(plot.title = element_text(size = 7, hjust = 0.5, face = "plain")))
  n_row <- if (length(titled) > 2) 2 else 1
  combo <- Reduce(`+`, titled) +
    patchwork::plot_layout(nrow = n_row, guides = "collect")
  ggsave(file.path(args$output_dir, "plot_modulescore_raw_dist_seg_control_composite.pdf"),
         combo, width = 19.05, height = if (n_row == 2) 10 else 5.5,
         units = "cm", device = "pdf")
  cat(sprintf("Wrote composite figure (A4 width, %d panels, %d row(s)): %s\n",
              length(titled), n_row, paste(names(panels), collapse = " | ")))
}

cat("\nAll outputs written to:", args$output_dir, "\n")
cat("Done.\n")
