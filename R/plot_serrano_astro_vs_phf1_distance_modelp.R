#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_serrano_astro_vs_phf1_distance_modelp.R
#
# Figure panels: 6B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# plot_serrano_astro_vs_phf1_distance_modelp.R
#
# The Cameron astrocyte-subcluster panel (Fig. 6b; cameron_subclusters x Astro x log in
# R/plot_modulescore_vs_phf1_distance_modelp.R) re-run with the Serrano-Pozo et al. 2024
# (Nat Neurosci) human astrocyte state markers in place of Cameron_2024.
#
# EVERYTHING EXCEPT THE GENE SETS IS COPIED VERBATIM from the parent script's Cameron path:
# same Seurat object, same Astro subset, same PHF1-negative filter, same AddModuleScore call
# (SCT, ctrl = 100, nbin = 24, seed = 42, nbin = 15 fallback), same 10-gene panel minimum,
# same 1000 um cap, same Set 1 model (docs/MODELS.md), same BH-within-group significance, same
# figure geometry (1.6 in panel, 2.1 in high), same colours (MODULE_COLOURS by position),
# same source-data and stats schema. Only these differ:
#   * gene sets  : h_astro_SerranoPozo.csv (columns state, Gene; one row per state x gene),
#                  used AS SUPPLIED -- no further filtering or ranking here -- except that
#                  astMic and astNeu are dropped (--exclude_states): the paper describes
#                  them as doublet artefacts, not astrocyte states.
#   * group name : serrano_states       (file names: ..._serrano_states_Astro_...)
#   * legend     : "Serrano-Pozo Astro"
#   * output dir : its own, so the Fig. 6b run directory and its shared tables are untouched.
#
# The parent script's model is Set 1:
#   score ~ log(dist)_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)
#   dist_scaled = log(dist_to_phf1_um) / sd(log(dist_to_phf1_um)), within Astro, within cap.
#
# Model types: "log" only (the reported fit, = Fig. 6b). "decay" is NOT run here: it is gated
# on the log model and needs the SEG baseline scoring. Also written, as in the parent: the raw
# model-free rolling mean overlay and the per-donor rolling-mean view.
#
# OUTPUT, under --output_dir:
#   plot_modulescore_dist_serrano_states_Astro_log.pdf                 -- the Fig. 6b analogue
#   source_data_modulescore_dist_serrano_states_Astro_log_fit.tsv
#   source_data_modulescore_dist_serrano_states_Astro_log_persample.tsv
#   stats_modulescore_dist_serrano_states_Astro_log.txt
#   plot_modulescore_raw_dist_serrano_states_Astro.pdf (+ _rollmean.tsv, stats .txt)
#   plot_modulescore_bydonor_dist_serrano_states_Astro.pdf (+ _rollmean.tsv, stats .txt)
#   stats_modulescore_geneset_overlap.tsv, stats_modulescore_geneset_members.tsv,
#   stats_modulescore_lmm_coeffs.tsv
#
# Run (interactive HPC session, from the project root):
#   Rscript R/plot_serrano_astro_vs_phf1_distance_modelp.R

##  ............................................................................
##  Packages + setup                                                        ####
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(argparse)
})

setwd("<PROJECT_ROOT>/phf1_v2")
source("R/palettes.R")             # fig_theme, celltype_palette, ...
source("R/rollmean_by_donor.R")    # per-donor rolling means, coloured by Braak

set.seed(42)

##  ............................................................................
##  Arguments                                                               ####
parser <- ArgumentParser()
parser$add_argument("--seu",
  default = "seu_PHF1.rds",
  help = "Path to seu_PHF1.rds [default: seu_PHF1.rds]")
parser$add_argument("--serrano_csv",
  default = "h_astro_SerranoPozo.csv",
  help = "Serrano-Pozo astrocyte state markers; columns 'state' and 'Gene', long format")
parser$add_argument("--exclude_states", default = "astMic,astNeu",
  help = paste("Comma-separated Serrano-Pozo states NOT analysed. Default astMic,astNeu:",
               "described in the paper as doublet artefacts, not astrocyte states"))
parser$add_argument("--output_dir",
  default = "plots/module_score_serrano_astro_vs_phf1_distance_modelp_1000um",
  help = "Output directory")
parser$add_argument("--seed",    type = "integer", default = 42,
  help = "Seed (AddModuleScore + provenance) [default: 42]")
parser$add_argument("--max_dist_um", type = "double", default = 1000,
  help = "If > 0, restrict modelled cells to dist_to_phf1_um <= this (um); 0 = no cap")
parser$add_argument("--window_um", type = "double", default = 75,
  help = "Rolling-mean window WIDTH (um) for the raw model-free overlay plots [default: 75]")
args <- parser$parse_args()

MAX_DIST <- if (args$max_dist_um > 0) args$max_dist_um else Inf
if (is.finite(MAX_DIST))
  cat(sprintf("Distance cap: modelling cells with dist_to_phf1_um <= %g um\n", MAX_DIST))

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Constants (identical to the parent script) -----------------------------------
MIN_GENES_PRESENT <- 10                         # drop signatures with fewer panel genes
CTRL              <- 100
NBIN              <- 24
RING_BREAKS       <- c(0, 50, 100, 200, 300, 500, 700)
SIG_ALPHA         <- 0.05                        # significance threshold (model FDR)
HALF_WINDOW       <- args$window_um / 2          # rolling-mean half-window (um) for raw plots

GROUP        <- "serrano_states"
CT           <- "Astro"
LEGEND_TITLE <- "Serrano-Pozo Astro"

# Same palette, same positional assignment as the parent's 8-module glia panels. With 9
# states the 9th module takes the pale yellow the parent reserves for position 9.
MODULE_COLOURS <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#000000",  # black
  "#999933",  # olive (Tol)
  "#F0E442",  # yellow
  "#332288"   # indigo (Tol)
)

# Ring labels / midpoints
.ring_lower <- head(RING_BREAKS, -1)
.ring_upper <- tail(RING_BREAKS, -1)
RING_LAB    <- sprintf("%g-%g", .ring_lower, .ring_upper)
RING_MID    <- setNames((.ring_lower + .ring_upper) / 2, RING_LAB)

##  ............................................................................
##  Helpers (verbatim from R/plot_modulescore_vs_phf1_distance_modelp.R)    ####

pick_col <- function(meta, candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0)
    stop(sprintf("None of the %s columns (%s) found in seu meta.data",
                 what, paste(candidates, collapse = ", ")))
  hit[1]
}

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

score_group <- function(seu_ct, present_sets, prefix, seed) {
  set.seed(seed)
  seu_ct <- tryCatch(
    AddModuleScore(seu_ct, features = present_sets, name = prefix,
                   assay = "SCT", ctrl = CTRL, nbin = NBIN, seed = seed),
    error = function(e) {
      cat("    AddModuleScore failed at nbin=", NBIN, " (", conditionMessage(e),
          ") — retrying at nbin=15\n", sep = "")
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
    dist_to_phf1_um = as.numeric(md$dist_to_phf1_um),
    percent_neg     = as.numeric(md$percent.neg),
    nUMI_log        = log2(as.numeric(md$nCount_RNA) + 1),
    sample_id       = as.character(md$sample_id),
    # Braak is carried only to COLOUR the per-donor rolling-mean lines; not a covariate.
    Braak           = if ("Braak" %in% colnames(md)) as.character(md$Braak) else NA_character_,
    Sex             = as.character(md[[sex_col]]),
    Age             = as.numeric(md[[age_col]]),
    PMI             = as.numeric(md[[pmi_col]]),
    stringsAsFactors = FALSE
  )

  scores <- do.call(cbind, lapply(names(feature_groups), function(prefix)
    score_group(seu_ct, feature_groups[[prefix]], prefix, seed)))

  cbind(base, scores)   # row-aligned: all derived from seu_ct@meta.data order
}

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

predict_grid <- function(model, df, dist_sd, log_distance, x_cap) {
  lo      <- if (log_distance) min(df$dist_to_phf1_um, na.rm = TRUE) else 0
  grid_um <- seq(lo, x_cap, length.out = 200)
  dt      <- if (log_distance) log(grid_um) else grid_um
  nd <- data.frame(
    dist_scaled = dt / dist_sd,
    nUMI_log    = mean(df$nUMI_log),
    percent_neg = mean(df$percent_neg),
    Sex         = factor(levels(df$Sex)[1], levels = levels(df$Sex)),
    Age_s       = mean(df$Age_s),
    PMI_s       = mean(df$PMI_s)
  )
  fe_form <- ~ dist_scaled + nUMI_log + percent_neg + Sex + Age_s + PMI_s
  X    <- model.matrix(fe_form, nd)
  beta <- lme4::fixef(model)
  X    <- X[, names(beta), drop = FALSE]
  fit  <- as.numeric(X %*% beta)
  V    <- as.matrix(vcov(model))
  se   <- sqrt(rowSums((X %*% V) * X))
  data.frame(dist_to_phf1_um = grid_um, fitted_score = fit,
             ci_lo = fit - 1.96 * se, ci_hi = fit + 1.96 * se)
}

roll_mean <- function(x, y, grid, half_window) {
  out <- lapply(grid, function(g) {
    idx <- which(x >= g - half_window & x <= g + half_window); n <- length(idx)
    if (n < 1) return(c(roll_mean = NA_real_, sem = NA_real_, n_window = 0))
    m <- mean(y[idx]); s <- if (n > 1) stats::sd(y[idx]) / sqrt(n) else NA_real_
    c(roll_mean = m, sem = s, n_window = n)
  })
  as.data.frame(do.call(rbind, out))
}

##  ............................................................................
##  Load Serrano-Pozo gene sets                                             ####
cat("Loading Serrano-Pozo astrocyte states from", args$serrano_csv, "...\n")
sp_raw <- read.csv(args$serrano_csv, stringsAsFactors = FALSE, check.names = FALSE,
                   fileEncoding = "UTF-8-BOM")
stopifnot("Serrano CSV must have columns 'state' and 'Gene'" =
            all(c("state", "Gene") %in% colnames(sp_raw)))
sp_raw$state <- trimws(sp_raw$state); sp_raw$Gene <- trimws(sp_raw$Gene)
sp_raw <- sp_raw[nzchar(sp_raw$state) & nzchar(sp_raw$Gene), ]
state_order <- unique(sp_raw$state)            # module order = order in the supplied file
serrano <- lapply(setNames(state_order, state_order),
                  function(s) unique(sp_raw$Gene[sp_raw$state == s]))
cat("  States:", paste(sprintf("%s (%d)", names(serrano), lengths(serrano)), collapse = ", "), "\n")

# Drop the states the paper itself calls doublet artefacts (astMic, astNeu by default).
EXCLUDED <- trimws(strsplit(args$exclude_states, ",", fixed = TRUE)[[1]])
EXCLUDED <- EXCLUDED[nzchar(EXCLUDED)]
if (length(setdiff(EXCLUDED, names(serrano))))
  stop("--exclude_states names states not in the CSV: ",
       paste(setdiff(EXCLUDED, names(serrano)), collapse = ", "))
serrano <- serrano[setdiff(names(serrano), EXCLUDED)]
cat("  Excluded (doublet artefacts per Serrano-Pozo 2024):",
    if (length(EXCLUDED)) paste(EXCLUDED, collapse = ", ") else "none", "\n")
cat("  Analysed:", paste(names(serrano), collapse = ", "), "\n")
stopifnot("More Serrano states than MODULE_COLOURS" = length(serrano) <= length(MODULE_COLOURS))

##  ............................................................................
##  Load Seurat object, determine panel, score Astro                        ####
cat("Loading Seurat object (this is large)...\n")
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "SCT"

stopifnot("celltype column missing"        = "celltype"        %in% colnames(seu@meta.data),
          "dist_to_phf1_um column missing" = "dist_to_phf1_um" %in% colnames(seu@meta.data),
          "percent.neg column missing"     = "percent.neg"     %in% colnames(seu@meta.data),
          "nCount_RNA column missing"      = "nCount_RNA"      %in% colnames(seu@meta.data),
          "sample_id column missing"       = "sample_id"       %in% colnames(seu@meta.data))

sex_col <- pick_col(seu@meta.data, c("Sex"), "Sex")
age_col <- pick_col(seu@meta.data, c("Age"), "Age")
pmi_col <- pick_col(seu@meta.data, c("PMI", "PostMortemInterval", "PostMortem_Interval"), "PMI")
cat(sprintf("Donor covariate columns: Sex=%s Age=%s PMI=%s\n", sex_col, age_col, pmi_col))

panel_genes <- rownames(seu[["SCT"]])
cat("CosMx panel genes (SCT assay):", length(panel_genes), "\n")

sp_prep <- prep_sets(serrano, panel_genes, MIN_GENES_PRESENT, "SerranoPozo_2024", GROUP)
write.table(sp_prep$overlap, file.path(args$output_dir, "stats_modulescore_geneset_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Gene-set / panel overlap:\n"); print(as.data.frame(sp_prep$overlap))

geneset_members <- bind_rows(lapply(names(serrano), function(m)
  tibble(group = GROUP, module = m, gene = serrano[[m]], on_panel = serrano[[m]] %in% panel_genes)))
write.table(geneset_members, file.path(args$output_dir, "stats_modulescore_geneset_members.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved gene-list members (on_panel = scored):",
    sum(geneset_members$on_panel), "of", nrow(geneset_members), "genes on panel\n")

cat("Scoring Astro (Serrano-Pozo states)...\n")
astro_df <- extract_celltype_data(
  seu, CT,
  feature_groups = list(SerranoPozo_ = sp_prep$present),
  seed = args$seed, sex_col = sex_col, age_col = age_col, pmi_col = pmi_col)

rm(seu); invisible(gc())

##  ............................................................................
##  Model frame (PHF1-negative, scaled covariates) -- as parent               ####
astro_df <- astro_df[!is.na(astro_df$dist_to_phf1_um) & !is.na(astro_df$percent_neg) &
                       !is.na(astro_df$nUMI_log) & !is.na(astro_df$Age) &
                       !is.na(astro_df$PMI) & !is.na(astro_df$Sex), ]
astro_df$Sex       <- droplevels(factor(astro_df$Sex))
astro_df$sample_id <- factor(astro_df$sample_id)
cat(sprintf("Modelled PHF1-negative cells (pre-cap): Astro=%d\n", nrow(astro_df)))

MODULES <- names(sp_prep$present)
COLS    <- setNames(MODULE_COLOURS[seq_along(MODULES)], MODULES)
x_cap_of <- function(df) if (is.finite(MAX_DIST)) MAX_DIST else
  as.numeric(quantile(df$dist_to_phf1_um, 0.99))

##  ............................................................................
##  Log-distance LMM per module (= Fig. 6b model)                           ####
cat(sprintf("\n=== %s / %s  (model: log_um%s) ===\n", GROUP, CT,
            if (is.finite(MAX_DIST)) sprintf(", <=%gum", MAX_DIST) else ""))
df_full <- astro_df
df_full$Age_s <- as.numeric(scale(df_full$Age))
df_full$PMI_s <- as.numeric(scale(df_full$PMI))
df <- if (is.finite(MAX_DIST)) df_full[df_full$dist_to_phf1_um <= MAX_DIST, ] else df_full
if (any(df$dist_to_phf1_um <= 0))
  stop("Non-positive dist_to_phf1_um among modelled cells: log() is undefined.")
dt      <- log(df$dist_to_phf1_um)
dist_sd <- sd(dt)
df$dist_scaled <- dt / dist_sd
x_cap <- x_cap_of(df)

fit_list <- list(); mod_stats <- list()
for (m in MODULES) {
  sub <- df[, c("cell_id", "dist_scaled", "nUMI_log", "percent_neg",
                "Sex", "Age_s", "PMI_s", "sample_id")]
  sub$score <- df[[m]]
  obs <- fit_obs_lmer(sub)
  if (is.null(obs)) { cat("  ", m, ": LMM failed, skipped\n"); next }
  grid <- predict_grid(obs$model, df, dist_sd, TRUE, x_cap)
  grid$module <- m
  fit_list[[m]] <- grid
  mod_stats[[m]] <- tibble(
    group = GROUP, celltype = CT, module = m, distance_scale = "log_um",
    slope_per_sd = obs$estimate, se = obs$se,
    ci_lo_per_sd = obs$estimate - 1.96 * obs$se,
    ci_hi_per_sd = obs$estimate + 1.96 * obs$se,
    slope_per_unit = obs$estimate / dist_sd,
    t = obs$t, df = obs$df, p_LMM = obs$p,
    singular = obs$singular, dist_sd = dist_sd
  )
}
if (length(mod_stats) == 0) stop("No Serrano-Pozo module could be modelled.")

mod_stats <- bind_rows(mod_stats)
mod_stats$lmm_padj    <- p.adjust(mod_stats$p_LMM, method = "BH")
mod_stats$significant <- mod_stats$lmm_padj < SIG_ALPHA
mod_stats$verdict     <- ifelse(is.na(mod_stats$significant), "not testable",
                         ifelse(mod_stats$significant,
                                "significant distance effect (model FDR<0.05)",
                                "n.s. (model FDR)"))
sig_basis <- "padj<0.05"

fit_df <- bind_rows(fit_list) %>%
  left_join(dplyr::select(mod_stats, module, significant, p_LMM, lmm_padj), by = "module") %>%
  mutate(module = factor(module, levels = MODULES),
         sig = factor(ifelse(significant, sig_basis, "ns"), levels = c(sig_basis, "ns")))

xlab <- expression("Distance to PHF1+ neuron (" * mu * "m)")   # mu via plotmath

p <- ggplot(fit_df, aes(dist_to_phf1_um, fitted_score, colour = module, group = module)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = module), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = sig, linewidth = sig)) +
  scale_colour_manual(values = COLS, name = LEGEND_TITLE, drop = FALSE) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c(sig_basis, "ns")),
                        name = NULL, drop = FALSE, limits = c(sig_basis, "ns")) +
  scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(sig_basis, "ns")),
                         name = NULL, drop = FALSE, limits = c(sig_basis, "ns")) +
  guides(colour = guide_legend(order = 1),
         linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
  labs(x = xlab, y = "Module score (fitted)", subtitle = NULL) +
  coord_cartesian(xlim = c(0, x_cap)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.text = element_text(size = 6),
        legend.key.width = grid::unit(20, "pt"),
        legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 5)) +
  scale_fill_manual(values = COLS, guide = "none")

g <- ggplot2::ggplotGrob(p)
pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
g$widths[pcol] <- grid::unit(1.6, "in")
ggsave(file.path(args$output_dir, sprintf("plot_modulescore_dist_%s_%s_log.pdf", GROUP, CT)),
       g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
       height = 2.1, units = "in", device = "pdf")

fit_out <- fit_df %>%
  transmute(celltype = CT, group = GROUP, distance_scale = "log_um",
            module = as.character(module), dist_to_phf1_um,
            fitted_score, ci_lo, ci_hi, significant)
write.table(fit_out, file.path(args$output_dir,
            sprintf("source_data_modulescore_dist_%s_%s_log_fit.tsv", GROUP, CT)),
            sep = "\t", quote = FALSE, row.names = FALSE)

ring_long <- lapply(MODULES, function(m) {
  d <- df[, c("sample_id", "dist_to_phf1_um")]
  d$score     <- df[[m]]
  d$dist_ring <- cut(d$dist_to_phf1_um, breaks = RING_BREAKS, labels = RING_LAB,
                     include.lowest = TRUE, right = TRUE)
  d <- d[!is.na(d$dist_ring), ]
  d %>% group_by(sample_id, dist_ring) %>%
    summarise(mean_score = mean(score), sd_score = sd(score), n_cells = dplyr::n(),
              .groups = "drop") %>%
    mutate(celltype = CT, group = GROUP, module = m,
           ring_mid_um = RING_MID[as.character(dist_ring)])
}) %>% bind_rows() %>%
  dplyr::select(celltype, group, module, sample_id, dist_ring, ring_mid_um,
                mean_score, sd_score, n_cells)
write.table(ring_long, file.path(args$output_dir,
            sprintf("source_data_modulescore_dist_%s_%s_log_persample.tsv", GROUP, CT)),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(args$output_dir, sprintf("stats_modulescore_dist_%s_%s_log.txt", GROUP, CT)))
cat("Module score vs distance to PHF1+ neuron (MODEL-FDR version, no permutation)\n")
cat("Group:", GROUP, " | celltype:", CT, " | distance scale: log_um\n")
cat("Gene sets: Serrano-Pozo et al. 2024 Nat Neurosci astrocyte states, from", args$serrano_csv,
    "(used as supplied).\n")
cat("States excluded before scoring (described in the paper as doublet artefacts):",
    if (length(EXCLUDED)) paste(EXCLUDED, collapse = ", ") else "none", "\n")
cat("  BH correction is therefore across the", length(MODULES), "analysed states only.\n")
cat("Everything else is identical to the Cameron panel (Fig. 6b) in\n")
cat("  R/plot_modulescore_vs_phf1_distance_modelp.R (Set 1, docs/MODELS.md).\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("LMM (log distance) (per module):\n")
cat("  score ~ dist_to_phf1_um_scaled + nUMI_log + percent_neg + Sex + Age + PMI + (1|sample_id)\n")
cat("  Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
    " seed=", args$seed, ").\n", sep = "")
cat("  Distance = log(um) / sd(log(um)), not centred (dist_sd =", round(dist_sd, 4), ").\n")
if (is.finite(MAX_DIST)) {
  cat("  Distance cap: cells with dist_to_phf1_um <=", MAX_DIST, "um.\n")
  cat("  Cells modelled (within cap):", nrow(df), " of ", nrow(df_full),
      " PHF1-negative | donors:", dplyr::n_distinct(df$sample_id), "\n")
} else {
  cat("  Cells (PHF1-negative):", nrow(df), " | donors:", dplyr::n_distinct(df$sample_id), "\n")
}
cat("  Donor random intercept: (1|sample_id).\n\n")
cat("Gene-set / panel overlap (sets with <", MIN_GENES_PRESENT, "panel genes are not scored):\n")
print(as.data.frame(sp_prep$overlap))
cat("\nSignificance = BH-adjusted model p (lmm_padj) within this group; padj <", SIG_ALPHA,
    " -> significant.\n\n")
cat("Effect size = slope per SD of log-distance, with 95% Wald CI (ci_lo/hi_per_sd).\n")
cat("  Negative slope = score higher NEAR PHF1+ neurons.\n\n")
cat("Per-module results:\n")
print(as.data.frame(mod_stats %>% dplyr::select(
  module, slope_per_sd, ci_lo_per_sd, ci_hi_per_sd, t, df, p_LMM, lmm_padj,
  significant, singular, verdict)))
cat("\nsessionInfo():\n"); print(sessionInfo())
sink()

write.table(mod_stats %>%
              dplyr::select(group, celltype, module, distance_scale, slope_per_sd, se,
                            ci_lo_per_sd, ci_hi_per_sd, slope_per_unit, t, df, p_LMM, lmm_padj,
                            significant, singular),
            file.path(args$output_dir, "stats_modulescore_lmm_coeffs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Done %s/%s_log: %d modules, %d significant (model FDR)\n",
            GROUP, CT, nrow(mod_stats), sum(mod_stats$significant, na.rm = TRUE)))

##  ............................................................................
##  Raw (model-free) rolling-mean overlay -- as parent make_glia_raw_plot()   ####
grid_r  <- seq(0, x_cap, length.out = 200)
sig_lut <- setNames(mod_stats$significant %in% TRUE, mod_stats$module)
sig_lut <- sig_lut[MODULES]; sig_lut[is.na(sig_lut)] <- FALSE; names(sig_lut) <- MODULES

roll_df <- bind_rows(lapply(MODULES, function(m) {
  y <- df[[m]]; ok <- !is.na(y) & !is.na(df$dist_to_phf1_um)
  data.frame(module = m, roll_mean(df$dist_to_phf1_um[ok], y[ok], grid_r, HALF_WINDOW),
             dist_to_phf1_um = grid_r)
}))
roll_df$module <- factor(roll_df$module, levels = MODULES)
roll_df$sig    <- factor(ifelse(sig_lut[as.character(roll_df$module)], sig_basis, "ns"),
                         levels = c(sig_basis, "ns"))

p_raw <- ggplot(roll_df, aes(dist_to_phf1_um, roll_mean, colour = module, fill = module, group = module)) +
  geom_ribbon(aes(ymin = roll_mean - 1.96 * sem, ymax = roll_mean + 1.96 * sem),
              alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = sig, linewidth = sig)) +
  scale_colour_manual(values = COLS, name = LEGEND_TITLE, drop = FALSE) +
  scale_fill_manual(values = COLS, guide = "none", drop = FALSE) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), c(sig_basis, "ns")),
                        name = "Log-model FDR", drop = FALSE, limits = c(sig_basis, "ns")) +
  scale_linewidth_manual(values = setNames(c(0.8, 0.4), c(sig_basis, "ns")),
                         name = "Log-model FDR", drop = FALSE, limits = c(sig_basis, "ns")) +
  guides(colour = guide_legend(order = 1),
         linetype = guide_legend(order = 2), linewidth = guide_legend(order = 2)) +
  labs(x = xlab, y = "Module score") +
  coord_cartesian(xlim = c(0, x_cap)) +
  theme_classic(base_size = 8) + fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.text = element_text(size = 6),
        legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(8, "pt"),
        legend.box.spacing = grid::unit(3, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 5))

g <- ggplot2::ggplotGrob(p_raw)
pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
g$widths[pcol] <- grid::unit(1.6, "in")
ggsave(file.path(args$output_dir, sprintf("plot_modulescore_raw_dist_%s_%s.pdf", GROUP, CT)),
       g, width = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
       height = 2.1, units = "in", device = "pdf")

write.table(roll_df %>% mutate(group = GROUP, celltype = CT, window_um = args$window_um,
                               log_significant = sig == sig_basis) %>%
              dplyr::select(group, celltype, module, window_um, log_significant,
                            dist_to_phf1_um, roll_mean, sem, n_window),
            file.path(args$output_dir,
                      sprintf("source_data_modulescore_raw_dist_%s_%s_rollmean.tsv", GROUP, CT)),
            sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(args$output_dir, sprintf("stats_modulescore_raw_dist_%s_%s.txt", GROUP, CT)))
cat("Raw (model-free) rolling-mean module scores vs distance to PHF1+ neuron\n")
cat("Group:", GROUP, " | celltype:", CT, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Rolling mean +/- 95% CI (mean +/- 1.96*SEM) of per-cell RAW AddModuleScore, per module\n")
cat("(NOT z-scored; native scale, shared axis).\n")
cat("Window:", args$window_um, "um; distance cap:",
    if (is.finite(MAX_DIST)) MAX_DIST else "none", "um; grid: 200 points. NO PHF1+ reference (glia).\n")
cat("Line style = LOG-model significance (solid+thick = model FDR<0.05, dashed+thin = n.s.),\n")
cat("  taken from stats_modulescore_dist_", GROUP, "_", CT, "_log.txt.\n", sep = "")
cat("Cells:", nrow(df), " | donors:", dplyr::n_distinct(df$sample_id), "\n")
cat("Modules (solid = log-significant):", paste(MODULES[sig_lut[MODULES]], collapse = ", "), "\n")
cat("Modules:", paste(MODULES, collapse = ", "), "\n\n")
cat("sessionInfo():\n"); print(sessionInfo())
sink()

##  ............................................................................
##  Per-donor rolling means -- as parent run_group_bydonor()                 ####
df_d <- df[is.finite(df$dist_to_phf1_um), ]
long <- rollmean_by_donor(df_d, MODULES, x_cap = x_cap, half_window = HALF_WINDOW)
if (is.null(long)) {
  cat("  [skip bydonor] no donor had enough cells\n")
} else {
  long$group <- GROUP; long$celltype <- CT
  ncol_facets <- if (length(MODULES) > 4) 4 else length(MODULES)
  p_d <- plot_rollmean_by_donor(long, xlab, "Module score (rolling mean)", x_cap,
                                facet_col = "module", facet_levels = MODULES,
                                ncol = ncol_facets)
  if (!is.null(p_d)) {
    g <- ggplot2::ggplotGrob(p_d)
    pcol <- unique(g$layout$l[grepl("^panel", g$layout$name)])
    g$widths[pcol] <- grid::unit(1.35, "in")
    nrow_facets <- ceiling(length(MODULES) / ncol_facets)
    ggsave(file.path(args$output_dir,
                     sprintf("plot_modulescore_bydonor_dist_%s_%s.pdf", GROUP, CT)), g,
           width  = grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE),
           height = 0.5 + 1.5 * nrow_facets, units = "in", device = "pdf")
  }
  write.table(long %>% dplyr::select(group, celltype, module, sample_id, Braak, window_um,
                                     dist_to_phf1_um, roll_mean, sem, n_window, n_cells_donor),
              file.path(args$output_dir,
                sprintf("source_data_modulescore_bydonor_dist_%s_%s_rollmean.tsv", GROUP, CT)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  sink(file.path(args$output_dir, sprintf("stats_modulescore_bydonor_dist_%s_%s.txt", GROUP, CT)))
  cat("Per-donor rolling-mean module scores vs distance to PHF1+ neuron\n")
  cat("Group:", GROUP, " | celltype:", CT, "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("Outcome = Seurat AddModuleScore (SCT data slot; ctrl=", CTRL, " nbin=", NBIN,
      " seed=", args$seed, ").\n", sep = "")
  cat("Modules:", paste(MODULES, collapse = ", "), "\n\n")
  rollmean_by_donor_stats(long, x_cap, HALF_WINDOW,
                          n_donors_total = dplyr::n_distinct(df_d$sample_id))
  cat("\nColour = Braak stage (braak_line_palette in R/rollmean_by_donor.R).\n")
  cat("\nsessionInfo():\n"); print(sessionInfo())
  sink()
}

cat("\nAll outputs written to:", args$output_dir, "\n")
cat("Done.\n")
