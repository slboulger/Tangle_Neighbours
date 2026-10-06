# ---------------------------------------------------------------------------
# build_aucell_annotation_tables.R
#
# Produces: Table S2, plus the supporting annotation tables behind the Methods counts
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
#
# WHAT THIS DOES
#   Assembles the cell-type annotation tables from the outputs of the two AUCell rounds
#   (AUCell_1.R: broad classes; AUCell_Neuron.R: neuronal subtypes). Cell types were
#   assigned by AUCell scoring of marker gene sets derived from an snRNA-seq entorhinal
#   cortex reference; no anchor-based label transfer is involved.
#
#   Table S2    - the marker gene sets of both rounds (TableS2_aucell_marker_gene_sets.tsv)
#   supporting  - set sizes, AUCell thresholds, assignment parameters, reference-label
#                 mapping, published-name crosswalk, and cells assigned per round, per
#                 final label and per donor (aucell_*.tsv). These document the Methods
#                 numbers (126,414 assigned; 62,742 neurons; 47,209 subtyped) and are not
#                 numbered supplementary tables.
#
#   The annotation is upstream of the manual tangle call, so none of these tables depend
#   on PHF1.
#
# INPUTS   AUCell_1/, AUCell_Neuron/ (gene sets, thresholds, per-cell calls),
#          PHF1/seu_coords.csv (the analysed cell set),
#          celltype_sce/Exc-ET-L5-SPON1-FGD4_sce.qs (rownames = the 6,175-gene panel)
# OUTPUTS  supp_tables/TableS2_aucell_marker_gene_sets.tsv, supp_tables/aucell_*.tsv
#
# Notes:
#   - The gene sets are saved as bare vectors, so no Wilcoxon statistics accompany them.
#   - The two rounds are not identically parameterised (marker cap 100 vs 60, min log2FC
#     0.25 vs 1, threshold and score floor 0.035 vs 0.027, gamma 0.40 vs 0.45, alpha 0.10
#     vs 0.11); see aucell_parameters.tsv.
#   - prior_offset_derived is back-calculated as (threshold - threshold_adj) / alpha. It is
#     exact where at_floor == FALSE and undefined where at_floor == TRUE.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr) })

setwd("<PROJECT_ROOT>/phf1_v2")
outd <- "supp_tables"; dir.create(outd, showWarnings = FALSE)
wr <- function(df, f) { write.table(df, file.path(outd, f), sep = "\t", quote = FALSE,
                                    row.names = FALSE, na = "NA"); cat("  wrote", f, "-", nrow(df), "rows\n") }

panel <- rownames(qs::qread("celltype_sce/Exc-ET-L5-SPON1-FGD4_sce.qs"))
stopifnot(length(panel) == 6175)

## published-name crosswalk (R/rename_celltypes.R)
exc_recode <- c("Exc-IT-L3-5"="Exc-IT-L3-5-CHGA-IL1RAPL2", "Exc-IT-L2-3"="Exc-IT-L2-3-CBLN2-HOPX",
                "Exc-CT-L6"="Exc-CT-L6-SYNPO2-SEMA3E", "Exc-ET-L5"="Exc-ET-L5-SPON1-FGD4",
                "Exc-IT-L6"="Exc-IT-L6-CTXN1-ERC2")
pub <- function(x) ifelse(x %in% names(exc_recode), exc_recode[as.character(x)], as.character(x))

## ---------------- Table S2: marker gene sets ----------------
cat("Table S2 + set sizes\n")
mk_set <- function(path, rnd, rlab) {
  gs <- readRDS(path)
  bind_rows(lapply(names(gs), function(l) data.frame(
    round = rnd, round_description = rlab, aucell_label = l, published_name = pub(l),
    rank_in_set = seq_along(gs[[l]]), gene = gs[[l]],
    on_cosmx_panel = gs[[l]] %in% panel, stringsAsFactors = FALSE)))
}
s2a <- bind_rows(
  mk_set("AUCell_1/gene_sets_from_SCT_snRNA.rds",      1L, "Broad cell classes"),
  mk_set("AUCell_Neuron/gene_sets_from_SCT_snRNA.rds", 2L, "Neuronal subtypes"))
wr(s2a %>% select(-on_cosmx_panel), "TableS2_aucell_marker_gene_sets.tsv")
cat("   on-panel:", sum(s2a$on_cosmx_panel), "/", nrow(s2a), "\n")

s2b <- s2a %>% group_by(round, round_description, aucell_label, published_name) %>%
  summarise(n_genes = n(), n_on_panel = sum(on_cosmx_panel), .groups = "drop") %>%
  arrange(round, aucell_label)
wr(as.data.frame(s2b), "aucell_marker_set_sizes.tsv")

## ---------------- thresholds + parameters ----------------
cat("thresholds + parameters\n")
a <- read.csv("AUCell_1/AUCell_label_thresholds_SCT_with_priors.csv", stringsAsFactors = FALSE)
b <- read.csv("AUCell_Neuron/AUCell_label_thresholds_SCT_with_priors.csv", stringsAsFactors = FALSE)
# threshold_prior_run comes from an earlier round-1 run without the class-prior adjustment.
# No script in this repository writes it, so the column is NA when the file is absent.
plain_f <- "AUCell_1/AUCell_label_thresholds_SCT.csv"
plain <- if (file.exists(plain_f)) {
  read.csv(plain_f, stringsAsFactors = FALSE) %>% select(label, threshold_prior_run = threshold)
} else data.frame(label = character(), threshold_prior_run = numeric())
alpha <- c(`1` = 0.10, `2` = 0.11); floor_v <- c(`1` = 0.035, `2` = 0.027)
mkthr <- function(d, rnd) d %>% mutate(
  round = rnd, aucell_label = label, published_name = pub(label),
  threshold_floored = threshold, threshold_adj = threshold_adj,
  prior_offset_derived = (threshold - threshold_adj) / alpha[[as.character(rnd)]],
  at_floor = abs(threshold_adj - floor_v[[as.character(rnd)]]) < 1e-12) %>%
  select(round, aucell_label, published_name, threshold_floored, threshold_adj,
         prior_offset_derived, at_floor)
s2c <- bind_rows(mkthr(a, 1L), mkthr(b, 2L)) %>%
  left_join(plain, by = c("aucell_label" = "label")) %>%
  select(round, aucell_label, published_name, threshold_prior_run,
         threshold_floored, threshold_adj, prior_offset_derived, at_floor)
wr(as.data.frame(s2c), "aucell_thresholds.tsv")

s2d <- data.frame(
  parameter = c("reference label column","marker test","min fraction detected (min.pct)",
    "min log2 fold change","FDR cutoff for markers","markers retained per label (top_n)",
    "min genes per set","AUCell threshold method","threshold floor","threshold ceiling",
    "prior softening exponent (gamma)","prior adjustment magnitude (alpha)",
    "per-cell AUC score floor","assignment rule","unassigned label"),
  round_1_broad = c("broad_celltype","Seurat FindAllMarkers, one-vs-rest Wilcoxon, positive only",
    "0.20","0.25","0.05","100","10","AUCell_exploreThresholds (fallback: 90th percentile)",
    "0.035","max(0.99, 99th pct of AUC)","0.40","0.10","0.035",
    "highest AUC/threshold_adj ratio, requiring ratio >= 1 and AUC >= score floor","Unassigned"),
  round_2_neuronal = c("broad_subcluster","Seurat FindAllMarkers, one-vs-rest Wilcoxon, positive only",
    "0.20","1","0.05","60","10","AUCell_exploreThresholds (fallback: 90th percentile)",
    "0.027","max(0.99, 99th pct of AUC)","0.45","0.11","0.027",
    "highest AUC/threshold_adj ratio, requiring ratio >= 1 and AUC >= score floor","Unassigned Neuron"),
  source = c("R/AUCell_1.R:55-63 / R/AUCell_Neuron.R:52-78","both scripts",
    "AUCell_1.R:72 / AUCell_Neuron.R:87","AUCell_1.R:73 / AUCell_Neuron.R:88",
    "AUCell_1.R:71 / AUCell_Neuron.R:86","AUCell_1.R:70 / AUCell_Neuron.R:85",
    "AUCell_1.R:69 / AUCell_Neuron.R:84","AUCell_1.R:160 / AUCell_Neuron.R:175",
    "AUCell_1.R:173 / AUCell_Neuron.R:188","AUCell_1.R:174 / AUCell_Neuron.R:189",
    "AUCell_1.R:191 / AUCell_Neuron.R:206","AUCell_1.R:195 / AUCell_Neuron.R:210",
    "AUCell_1.R:205 / AUCell_Neuron.R:220","AUCell_1.R:214 / AUCell_Neuron.R:229",
    "AUCell_1.R:224 / AUCell_Neuron.R:239"), stringsAsFactors = FALSE)
wr(s2d, "aucell_parameters.tsv")

## ---------------- reference mapping + crosswalk ----------------
cat("reference mapping + crosswalk\n")
r2map <- list(
  "Exc-IT-L2-3" = c("Exc-L2-3-CUX2-CALB1","Exc-L2-CUX2-PDGFD","Exc-L3-PCP4-CALB1","Exc-L2-RELN-BMPR1B"),
  "Exc-IT-L3-5" = c("Exc-L5-RORB-TLL1","Exc-L5-RORB-TPBG"),
  "Exc-ET-L5"   = c("Exc-L5-BCL11B-ADRA1A"),
  "Exc-CT-L6"   = c("Exc-L6-TLE4-SULF1","Exc-L5-6-TLE4-NXPH2"),
  "Exc-IT-L6"   = c("Exc-L6-THEMIS-CDH13"),
  "Inh-PVALB"   = c("Inh-PVALB-MYO5B","Inh-PVALB-UNC5B"),
  "Inh-SST"     = c("Inh-SST-NPY"),
  "Inh-VIP"     = c("Inh-VIP-RELN"),
  "Inh-LAMP5"   = c("Inh-LAMP5-RELN","Inh-RELN,SST"))
s2e <- bind_rows(lapply(names(r2map), function(l) data.frame(
  round = 2L, reference_subcluster = r2map[[l]], aucell_label = l, published_name = pub(l),
  mapping_type = "explicit", stringsAsFactors = FALSE)))
s2e_r1 <- data.frame(round = 1L,
  reference_subcluster = c("Astro","Oligo","OPC","Endo","VLMC","Micro",
                           "any cluster_celltype matching ^Exc","any cluster_celltype matching ^Inh","all others"),
  aucell_label = c("Astro","Oligo","OPC","Endo","VLMC","Micro","Glutamatergic","GABAergic","Other"),
  mapping_type = c(rep("pass-through", 6), rep("regex rule", 2), "fallback"), stringsAsFactors = FALSE) %>%
  mutate(published_name = pub(aucell_label)) %>%
  select(round, reference_subcluster, aucell_label, published_name, mapping_type)
s2e <- bind_rows(s2e_r1, s2e)
wr(as.data.frame(s2e), "aucell_reference_label_mapping.tsv")

s2f <- data.frame(aucell_label = names(exc_recode), published_name = unname(exc_recode),
                  renamed = TRUE, stringsAsFactors = FALSE)
wr(s2f, "aucell_published_name_crosswalk.tsv")

## ---------------- assignment outcome ----------------
cat("assignment outcome\n")
coords <- read.csv("PHF1/seu_coords.csv", stringsAsFactors = FALSE)
r1 <- read.csv("AUCell_1/cosmx_AUCell_labels_SCT_propAdj.csv", stringsAsFactors = FALSE)
r2 <- read.csv("AUCell_Neuron/cosmx_AUCell_labels_SCT_propAdj.csv", stringsAsFactors = FALSE)
an <- coords %>% select(cell_id, sample_id, celltype_final = celltype) %>%
  left_join(r1 %>% select(cell_id = cell, r1_label = final_label), by = "cell_id") %>%
  left_join(r2 %>% select(cell_id = cell, r2_label = final_label), by = "cell_id")
stopifnot(nrow(an) == 167176)
NEUR <- c("Glutamatergic","GABAergic")

s3a <- data.frame(
  round = c(1L, 2L),
  round_description = c("Broad cell classes","Neuronal subtypes"),
  cells_scored = c(nrow(an), sum(an$r1_label %in% NEUR)),
  cells_assigned = c(sum(an$r1_label != "Unassigned"),
                     sum(an$r1_label %in% NEUR & an$r2_label != "Unassigned Neuron")),
  cells_unassigned = c(sum(an$r1_label == "Unassigned"),
                       sum(an$r1_label %in% NEUR & an$r2_label == "Unassigned Neuron")),
  unassigned_label = c("Unassigned","Unassigned Neuron"), stringsAsFactors = FALSE) %>%
  mutate(pct_assigned = round(100 * cells_assigned / cells_scored, 2))
wr(s3a, "aucell_assignment_summary.tsv")

s3b <- an %>% count(celltype_final, name = "n_cells") %>%
  mutate(published_name = celltype_final,
         pct_of_analysed = round(100 * n_cells / sum(n_cells), 3)) %>%
  arrange(desc(n_cells)) %>% select(published_name, n_cells, pct_of_analysed)
wr(as.data.frame(s3b), "aucell_final_label_counts.tsv")

braak <- c("BBN00428931"="I","BBN_10231"="I","BBN00635914"="I",
           "BBN_9931"="IV","BBN00638047"="IV","BBN00628710"="IV",
           "BBN_9928"="VI","BBN_24895"="VI","BBN_9889"="VI")
stopifnot(all(unique(an$sample_id) %in% names(braak)))
s3c <- an %>% group_by(sample_id) %>% summarise(
  cells_analysed      = n(),
  r1_assigned         = sum(r1_label != "Unassigned"),
  r1_unassigned       = sum(r1_label == "Unassigned"),
  neurons             = sum(r1_label %in% NEUR),
  neurons_subtyped    = sum(r1_label %in% NEUR & r2_label != "Unassigned Neuron"),
  neurons_unassigned  = sum(r1_label %in% NEUR & r2_label == "Unassigned Neuron"),
  .groups = "drop") %>%
  mutate(Braak = braak[sample_id],
         pct_r1_assigned  = round(100 * r1_assigned / cells_analysed, 2),
         pct_subtyped     = round(100 * neurons_subtyped / neurons, 2)) %>%
  arrange(match(Braak, c("I","IV","VI")), sample_id) %>%
  select(sample_id, Braak, cells_analysed, r1_assigned, r1_unassigned, pct_r1_assigned,
         neurons, neurons_subtyped, neurons_unassigned, pct_subtyped)
wr(as.data.frame(s3c), "aucell_assignment_by_sample.tsv")

cat("\n--- checks against the manuscript ---\n")
chk <- function(l, g, w) cat(sprintf("  %-24s %8d  paper %8d  %s\n", l, g, w,
                                     ifelse(g == w, "MATCH", "*** DIFFERS ***")))
chk("round-1 assigned", s3a$cells_assigned[1], 126414)
chk("neurons annotated", s3a$cells_scored[2],  62742)
chk("neurons subtyped",  s3a$cells_assigned[2], 47209)
chk("panel genes",       length(panel), 6175)
cat("\nDONE\n")
