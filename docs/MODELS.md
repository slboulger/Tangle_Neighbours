# Model specifications

Which model backs which panel. Sets 1, 1p and 2 are the CosMx specifications; Set 3 is the
imaging mass cytometry (IMC) cohort, which is a different object keyed on its own donor column.

---

## The specifications

### Set 1 — cell-level, transcript-derived outcome

```
~ <predictor> + nUMI_log + percent_neg + Sex + Age + PMI + (1 | sample_id)
```

Fitted with `lmerTest::lmer(..., REML = TRUE)` for continuous outcomes (module scores) or
`variancePartition::dream` for gene-level differential expression.

- `nUMI_log = log2(nCount_RNA + 1)`
- `percent_neg` = `percent.neg` (negative-probe rate)
- `Age_s` / `PMI_s` = `as.numeric(scale(Age))` / `scale(PMI)`
- `Sex` is a factor.

**Panels:** Fig. 3A, 3C, 3D, 4B, 6B, 6E (and the CosMx series of 6C); Figs S1A-C, S2A-C,
S4A-E, S5A-C, S9A-F; Table S6.

### Set 1p — pseudobulk, transcript-derived outcome

```
~0 + <predictor> + nUMI_scaled + Sex + Age + PMI
```

with donor handled by `duplicateCorrelation` / `lmFit(block = sample_id)`, not a random effect.
`limma-trend` (`eBayes(trend = TRUE, robust = TRUE)`), CPM normalisation.

| Set 1 term | Set 1p | Why |
|---|---|---|
| `nUMI_log` | `nUMI_scaled` = `log2(nUMI)` | pseudocell-level aggregate |
| `(1 \| sample_id)` | `duplicateCorrelation` block | donor is the unit, not a level |
| `percent_neg` | absent | after aggregation it is a mean over ~10⁵ counts and carries almost no between-pseudocell information; it is already applied as a cell-removal criterion in QC |

`pseudocell_size_scale` is not a covariate: `nUMI_scaled` already carries it, since total
counts scale with the number of cells summed.

**Panels:** Fig. 2B; Table S4 (also consumed by Fig. 4C-E).

### Set 2 — image-derived outcome

```
~ <predictor> + Sex + Age + PMI + (1 | sample_id)
```

No transcript covariates. Conditioning an intensity or morphology outcome on transcript
abundance would adjust away part of the effect under test — PHF1+ cells carry ~50% more RNA,
so `nUMI_log` is a collider on the PHF1 contrast, not a nuisance term.

Image-derived covariates (`log_area`, `mask_px_log`, `dapi_z`, `depth_rel`, `edge_dist_um`,
`n_local`) are legitimate additions and are used as labelled sensitivities.

**Panels:** Fig. 4A, 5C; Figs S3 (distance fit, see *Documented exceptions*) and S7B.

### Set 3 — IMC (imaging mass cytometry, entorhinal cortex cohort)

```
~ <predictor> + Sex + Age_s + PMI_s + (1 | patient_id)
```

A different cohort and object (`<IMC_ROOT>/spe.rds`, 44 donors), so it has its own
specification. Three things differ from Set 2, all structural:

- **`patient_id` is the donor**; `sample_id` is the ROI and is not modelled. Set 2's
  `(1 | sample_id)` is the CosMx section, so `(1 | patient_id)` is its analogue here.
- **No transcript covariates exist**, so `percent_neg` / `nUMI_log` cannot arise.
- `Age_s` / `PMI_s` are `as.numeric(scale(...))`; `Sex` is a factor.

**Sensitivity terms, never in a headline:** `plaque` (`Matched_4G8_40`, nucleus within a 40 µm
dilated 4G8+ mask), `scaled_Y` / `depth_s` (relative cortical depth), `subcluster`, `log_area`,
`n50` (local neuron density). Random intercept, not random slope, throughout.

**Distance:** cap 300 µm, transform `log(dist_um) / sd(log(dist_um))` (see *Distance
transform*), with `dist_sd` recomputed per cell set, so "per s.d." magnitudes are not comparable
across panels.

**Panels:** Fig. 2E, 3E, 5A, 5B and the IMC series of 6C; Figs S7A and S8A-E.

| Panel | Script | n |
|---|---|---|
| Fig. 2E state markers, tangle-bearing vs tangle-free (all excitatory + subcluster) | `R/imc_phf1_state_markers_subsets.R` → `R/imc_phf1_state_markers_forest_panel.R` | 9,130 / 832 / 8,298 |
| Fig. 3E protein gradients | `R/imc_phf1_distance_markers_subsets.R` | 6,856 cells, 44 donors |
| Fig. 5A nucleus area vs distance (and S7A eccentricity) | `R/imc_phf1_morphology_dist_phf1_panel.R` | 26,653 cells, 44 donors |
| Fig. 5B total nuclear DNA vs distance | `R/imc_phf1_dna_dist_phf1_panel.R` | 26,653 cells |
| Fig. 6C GFAP vs distance (IMC series) | `R/imc_phf1_glia_distance.R`, drawn by `R/imc_cosmx_gfap_distance_panel.R` | 26,832 astrocytes, 44 donors |

**Within-donor z-scoring in the glia intensity panels.** The two intensity outcomes in
`R/imc_phf1_glia_distance.R` (including the Fig. 6C GFAP series) are z-scored within donor, so
that the plotted rolling curve and the fitted coefficient describe the same within-donor
gradient. Because the outcome then has no between-donor variance, the donor random-intercept
variance is estimated at zero and the fit is equivalent to pooled least squares on the
standardised outcome. The full covariate set is retained: z-scoring standardises the outcome,
not the predictor. Adding or removing Sex, Age and PMI moves the GFAP coefficient by 0.5%
(+0.05420 vs +0.05446).

#### Auxiliary regressions for the plaque-stratified sensitivity (Figs S8D, S8E)

`R/imc_amyloid_tangle_coupling.R` decomposes the change in the tangle-distance coefficient when
`plaque` is added as **γ × δ**, where γ is the plaque term of the plaque-adjusted model and δ is
an auxiliary regression of plaque status on the included predictor. δ is specified by three
rules:

1. **Identity link.** `plaque` enters the outcome models as a 0/1 term with an identity link,
   and the omitted-variable identity is a property of linear projection, so δ is a linear
   probability fit. A logistic fit is also reported, as an odds ratio describing the coupling,
   and is never multiplied by γ.
2. **Same right-hand side as the outcome model**, minus the omitted term (including any
   panel-specific terms such as `subcluster`), built programmatically from the panel metadata.
3. **Same weight matrix as the outcome fit**: δ is computed by whitened GLS at the outcome fit's
   own λ = τ²/σ², which makes the decomposition exact. Where a donor random effect is estimated at
   zero, λ ≈ 0 and δ reduces to pooled least squares, matching its γ.

The identity is asserted numerically (single-λ form to 1e-8; REML decomposition to 0.05 ×
SE(b1)). All decomposition arithmetic is on the raw estimate scale, and the reporting currency is
shift / SE(b1).

---

## Classification rule: module-anchored analyses take Set 1

Fits that use `phf1_intensity_p95_z` as the regression outcome but are anchored on a gene-set
module score take Set 1, because the question is about the module. In
`R/phf1_module_distance_intensity_adjust.R` (Fig. 4B, S4A-E) this covers the a-path
(`intensity ~ dist_scaled + covariates`) and the non-neuronal spillover fit. Classification
follows the predictor and the scientific question, not the outcome's modality; Fig. 4A, which
asks about intensity itself, takes Set 2.

## Classification rule: technical and image outcomes take Set 2

When the outcome is itself an image or technical quantity — channel intensity, segmentation
geometry, nuclear morphology, transcript density or library size — the specification is Set 2,
with no transcript covariates, whatever the outcome's modality:

1. Set 1 is unidentifiable when a transcript covariate is the outcome; it cannot appear on both
   sides.
2. Conditioning on transcript abundance would adjust away the effect under test (the collider
   argument in Set 2).

This covers Fig. 5C (rRNA channel) and Fig. S7B (nucleus area).

---

## Panel → script → specification

| Panel | Script | Specification |
|---|---|---|
| Fig. 1F | `R/plot_phf1_intensity_validation_neurons.R` (per-cell table from `R/plot_phf1_intensity_validation.R`) | ROC AUC; two-sided paired Wilcoxon signed-rank on donor means (n = 9) |
| Fig. 2A | `R/plot_phf1_neuron_density_by_celltype.R` | negative-binomial GLMM, `n_tangle ~ celltype + (1 \| sample_id) + offset(log(area))` (see *Documented exceptions*) |
| Fig. 2B | `R/deg_pb_category_phf1.r` (Table S4) → `R/plot_dotplot_pathway_cat_phf1.r` | Set 1p; enrichR over-representation |
| Fig. 2E | `R/imc_phf1_state_markers_subsets.R` → `R/imc_phf1_state_markers_forest_panel.R` | Set 3 (+ `subcluster`) |
| Fig. 3A | `R/deg_dream_linear_distance_phf1.r` (Table S6) → `R/nDEG_by_celltype_linear.R` | Set 1 |
| Fig. 3B | `R/pathway_enrichment_enrichR_linear.r` → `R/plot_dotplot_pathway_linear_distance_phf1.r` | enrichR over-representation |
| Fig. 3C, 3D | `R/plot_reactome_stress_death_vs_phf1_distance.R` | Set 1 |
| Fig. 3E | `R/imc_phf1_distance_markers_subsets.R` | Set 3 |
| Fig. 4A | `R/phf1_intensity_gradient_neurons.R` | Set 2 |
| Fig. 4B | `R/phf1_module_distance_intensity_adjust.R` | Set 1 (+ per-cell PHF1 intensity) |
| Fig. 4C-E | `R/rrho2_decomposition_phf1.r` → `R/plot_rrho2_*` | RRHO2 on the Table S4 and S6 statistics; no new model |
| Fig. 5A | `R/imc_phf1_morphology_dist_phf1_panel.R` | Set 3 |
| Fig. 5B | `R/imc_phf1_dna_dist_phf1_panel.R` | Set 3 |
| Fig. 5C | `R/phf1_channel_intensity_exc.R` | Set 2 |
| Fig. 6A | `R/plot_dist_to_phf1_by_glia_vs_null.R` + `R/null_glia_distance_within_celltype.R` | label-permutation null within cell type and donor |
| Fig. 6B | `R/plot_serrano_astro_vs_phf1_distance_modelp.R` | Set 1 |
| Fig. 6C | `R/imc_phf1_glia_distance.R` (IMC) and the GFAP row of `R/deg_dream_linear_distance_phf1.r` (CosMx), drawn by `R/imc_cosmx_gfap_distance_panel.R` | Set 3 (IMC) / Set 1 (CosMx) |
| Fig. 6E | `R/plot_modulescore_vs_phf1_distance_modelp.R` | Set 1 |
| Fig. S1A-C | `R/run_deg_linear_distance.sh` → `R/collate_exclusion_radius.r`; `R/deg_dream_linear_distance_depth_phf1.r` → `R/collate_depth_covariate.r`; `R/deg_dream_linear_distance_3d_phf1.r` → `R/collate_3d_distance_deg.r` | Set 1 (+ `depth_s` in S1B; recomputed 3D distance in S1C) |
| Fig. S2A-C | `R/plot_reactome_stress_death_vs_phf1_distance.R`; `R/plot_reactome_gene_contribution.R` | Set 1 |
| Fig. S3 | `R/run_nn3.sh` → `R/nn3_dist_phf1_panel_null.R`, `R/nn3_null_deviation.R` | Set 2 + `edge_dist_um`, calibrated against label permutations (see *Documented exceptions*) |
| Fig. S4A-E | `R/phf1_module_distance_intensity_adjust.R` | Set 1 (+ per-cell PHF1 intensity) |
| Fig. S5A-C | `R/modulescore_seg_control_vs_phf1_distance.R` → `R/plot_seg_control_composite_matched.R` | Set 1 |
| Fig. S6 | `R/rrho2_decomposition_phf1.r` → `R/plot_dotplot_pathway_rrho2_phf1.r` | enrichR over-representation |
| Fig. S7A | `R/imc_phf1_morphology_dist_phf1_panel.R` | Set 3 |
| Fig. S7B | `R/phf1_morphology_dist_phf1_panel.R` | Set 2 |
| Fig. S8A-C | `R/imc_covariate_sensitivity.R` | Set 3 with each sensitivity term added in turn |
| Fig. S8D, S8E | `R/imc_amyloid_tangle_coupling.R` | Set 3 within plaque strata, plus the auxiliary δ regressions above |
| Fig. S9A-F | `R/plot_reactome_stress_death_vs_phf1_distance.R` | Set 1 |

Panels 1B-1E, 2C, 2D, 2F, 5D and 6D are descriptive or image panels and fit no model.

---

## Documented exceptions

### Donor-level covariates aliased with a donor term

Sex, Age and PMI are exactly aliased with a per-donor intercept when donor enters as a fixed
effect. Including them yields a rank-deficient design, not better adjustment, so they are
dropped in:

| Script | Design | What is dropped |
|---|---|---|
| `R/decay_length_utils.R` | per-donor profile least squares, `0 + donor` | Sex/Age/PMI (aliased with the per-donor plateau). Keeps `nUMI_log + percent_neg` |
| `R/null_replica_nn3.R` (fast projection for the Fig. S3 null) | donor as fixed dummies | Sex/Age/PMI (aliased with the donor dummies). The observed coefficient from the same projection is checked against the mixed-model coefficient in the stats log |

### Within-donor contrasts fitted without Sex, Age and PMI

Two group contrasts compare cells within each donor and are fitted with a donor random intercept
but without the donor-level covariates. Because the predictor varies within donor and Sex, Age
and PMI are constant within donor, the covariates can only explain between-donor differences
that the random intercept already absorbs. Refitting both with the covariates added confirms
this:

- **Fig. S3A, tangle-bearing strip** (`R/nn3_phf1_vs_neg.R`, per subtype:
  `log(nn3) ~ grp + edge_dist_um + (1 | sample_id)`). Adding Sex, Age and PMI changes the
  Exc-IT-L2-3-CBLN2-HOPX estimate by 0.1% (+5.63% vs +5.62% spacing) and leaves the BH-adjusted
  p (0.13) and the plotted call unchanged; Exc-IT-L3-5-CHGA-IL1RAPL2 likewise.
- **Fig. 2A** (`R/plot_phf1_neuron_density_by_celltype.R`, negative-binomial GLMM across cell
  types). Adding Sex, Age and PMI leaves every leave-one-out call unchanged.

The companion contrasts in `R/phf1_morphology_exc.R` (run by `R/run_morphology_observed.sh`) are
also fitted without the donor-level covariates; they are not drawn in any panel. The Fig. S7B
contrast drawn in the panel is fitted by `R/phf1_morphology_dist_phf1_panel.R` under Set 2.

### Edge distance in the 3-NN spacing models (Fig. S3)

The 3-NN spacing outcome is a nearest-neighbour distance, so a neuron near a field-of-view edge
has some of its true neighbours unimaged and its spacing is inflated. The spacing models
(`R/nn3_phf1_vs_neg.R`, `R/nn3_over_phf1_distance.R`) therefore add `edge_dist_um` (distance to
the FOV edge, capped) to the Set 2 right-hand side. It is an image-derived covariate of the kind
Set 2 allows.

### The PHF1 module is a gene-set definition, not an inferential model

The PHF1 module (Fig. 4B, S4, Table S9) is defined by `R/find_phf1_markers_by_celltype.R`:
within each cell type, tangle-bearing vs tangle-free cells by Wilcoxon rank-sum test on the SCT
data slot (`Seurat::FindMarkers`, min.pct 0.2, log2FC threshold 0.25), keeping genes higher in
tangle-bearing cells (log2FC > 0.25) at adjusted p < 0.05. This step selects genes; it carries no covariates or
donor term because its p-values are not reported as results. All inference about the module is
the Set 1 model of module score over distance in tangle-free cells. The donor-aware differential
expression of tangle-bearing vs tangle-free neurons is the Set 1p analysis in Table S4.

---

## Permutation nulls

A null refits the specification of the observed model it calibrates, on the same cells. Label
permutations are generated once (`R/generate_phf1_null_labels.r`, 1,000 labellings, seed 42) and
shared by the Fig. 6A and Fig. S3 nulls; each permutation preserves the number of tangle-bearing
neurons within each subtype and donor. The Fig. S3 null refits the distance coefficient with the
same covariates as the observed model (`R/null_replica_nn3.R`).

---

## Distance transform

Uniform across every linear-in-distance model:

```r
dist_scaled <- log(dist_to_phf1_um) / sd(log(dist_to_phf1_um))
```

Natural log, divided by SD only, never centred, with `dist_sd` computed within each cell type
after the distance cap. Distance is strictly positive for every modelled cell, and scripts
`stop()` on a non-positive distance rather than absorbing it with an offset.

Two documented variants estimate a quantity in microns and so use the raw axis:

| Variant | Where | Why |
|---|---|---|
| raw µm inside `s(dist, k = 5)` | GAM paths in the module-score drivers | the smooth is estimated on the micron axis |
| raw µm inside `exp(-d / λ)` | `R/decay_length_utils.R` | λ has units of microns; a log-rescaled axis has no length constant |
