# Analysis code for Boulger et al. - Neuronal and glial transcriptional states vary with proximity to tangle-bearing neurons in the human entorhinal cortex

Code for the spatial analysis of the neuronal and glial microenvironment around
tau-tangle-bearing neurons in human entorhinal cortex, combining 6,000-plex CosMx
spatial transcriptomics with post-run PHF1 (phospho-tau) immunofluorescence, and
imaging mass cytometry (IMC) from a separate, larger cohort.

This repository contains **only the scripts behind the figures and supplementary
tables of the final manuscript**. Exploratory and superseded analyses are not included.

## Before you start

- **All absolute paths are placeholders.** Nothing runs until you substitute them —
  see [PLACEHOLDERS.md](PLACEHOLDERS.md).
- **Donor IDs are UK Brain Banks Network (BBN) IDs.** Brain-bank case IDs have been
  removed and appear nowhere in this repository.
- **Nine CosMx donors are analysed, but ten sections were run.** The tenth donor was
  excluded because the section was found to be white matter rather than grey matter. It
  still appears in [`merge_seu.R`](R/merge_seu.R) (slide 2, FOVs 42-100) because its FOVs
  are part of the AtoMx export. It is removed in [`add_PHF1.R`](R/add_PHF1.R) by the filter
  that keeps donors with both PHF1+ and PHF1- cells, since it has no PHF1+ cells. No object
  or result after that step includes it.
- **No data is tracked here.** The Seurat/SCE objects, images and result tables live
  with the data release, not in git.

## Two cohorts, two conventions

| | CosMx | IMC |
|---|---|---|
| Donors | 9 (3 each Braak I / IV / VI) | 44 AD (Braak III–VI); Fig. 4 adds 18 non-diseased controls (Braak 0–I) |
| Tangle call | manual annotation of PHF1 post-stain | IMC PHF1 channel (Otsu threshold) |
| Distance cap | 1,000 µm | 300 µm |
| Model | `~ log-distance + nUMI + percent_neg + Sex + Age + PMI + (1\|donor)` | `~ predictor + Sex + Age + PMI + (1\|donor)` |

Distance to the nearest tangle-bearing neuron is always transformed as
`log(d) / sd(log(d))` — natural log, scaled by SD only, **not** centred, with `sd`
computed within each cell type. See [docs/MODELS.md](docs/MODELS.md) for the full
model specifications and the documented exceptions.

### Where the IMC object comes from

The IMC scripts read one SpatialExperiment, `<IMC_ROOT>/spe.rds` (67 donors, 195 ROIs,
215,739 cells, 35 channels). The distance and state-marker panels subset it to the 44 Braak
III-VI donors; the neuronal-loss panels (Fig. 4) use the 62 donors with three ROIs (44 AD,
18 controls). It is the entorhinal cortex IMC dataset described previously (Boulger et al.,
2026, *Intrinsic molecular susceptibility underlies selective neuronal vulnerability in the
Alzheimer's disease entorhinal cortex*, bioRxiv,
[doi:10.64898/2026.06.22.733753](https://doi.org/10.64898/2026.06.22.733753)), segmented with `steinbock` v0.16.3 and built with
the IMC data analysis workflow of Windhager et al. (*Nat. Protoc.* 18, 3565-3613, 2023).
Tissue, antibody panel, acquisition, segmentation and phenotyping are described there.

The public, de-identified release of that object carries the donor as `BBN_ID` and has no
`patient_id` column. Every IMC model here keys its donor random effect on `patient_id`, so
set `spe$patient_id <- spe$BBN_ID` after loading it.

## Figure panels → scripts

Machine-readable version: [figure_panel_map.tsv](figure_panel_map.tsv).

### Main figures

| Panel | Panel shows | Data | Script(s) |
|---|---|---|---|
| **1A** | Cohort / workflow schematic | schematic | *not archived* |
| **1B** | Marker z-score heatmap, broad classes | CosMx | [`plot_aucell_heatmaps_renamed.R`](R/plot_aucell_heatmaps_renamed.R) |
| **1C** | Marker z-score heatmap, neuronal subtypes | CosMx | [`plot_aucell_heatmaps_renamed.R`](R/plot_aucell_heatmaps_renamed.R) |
| **1D** | Per-donor cell-type composition | CosMx | [`plot_celltype_composition_by_sample.R`](R/plot_celltype_composition_by_sample.R) |
| **1E** | Segmented FOV, DAPI + PHF1, arrowheads | CosMx | [`plot_phf1_celltype_overlay.py`](python/plot_phf1_celltype_overlay.py) |
| **1F** | PHF1 p95 intensity recovers manual calls in neurons | CosMx | [`plot_phf1_intensity_validation_neurons.R`](R/plot_phf1_intensity_validation_neurons.R)<br>[`plot_phf1_intensity_validation.R`](R/plot_phf1_intensity_validation.R) (per-cell table) |
| **2A** | Tangle-bearing neurons per mm2 by subtype | CosMx | [`plot_phf1_neuron_density_by_celltype.R`](R/plot_phf1_neuron_density_by_celltype.R) |
| **2B** | UMAP of Exc-IT-L2-3-CBLN2-HOPX coloured by pseudotime, tangle-bearing neurons marked | CosMx | [`plot_traj_umap_publication.R`](R/plot_traj_umap_publication.R)<br>[`run_traj_pipeline.sh`](R/run_traj_pipeline.sh) (embedding, trajectory search and fit)<br>[`cbln2_trajectory_selection_null.R`](R/cbln2_trajectory_selection_null.R) (permutation p) |
| **2C** | % tangle-bearing per pseudotime bin | CosMx | [`plot_cbln2_traj_phf1_by_bin_mean.R`](R/plot_cbln2_traj_phf1_by_bin_mean.R) |
| **2D** | Otero-Garcia Exc1/2 UP score, tangle-bearing vs tangle-free and across pseudotime bins | CosMx | [`plot_cbln2_traj_otero_by_bin.R`](R/plot_cbln2_traj_otero_by_bin.R)<br>[`prepare_fig2_contrast_genesets.R`](R/prepare_fig2_contrast_genesets.R) (gene list) |
| **2E** | Tangle-bearing CBLN2 neuron crop, Otero-Garcia Exc1/2 UP transcripts | CosMx image | [`plot_tangle_cell_crops.py`](python/plot_tangle_cell_crops.py) (see *Crop panels* below) |
| **2F** | Tangle-bearing CBLN2 neuron crop, microtubule (de)polymerisation transcripts | CosMx image | [`plot_tangle_cell_crops.py`](python/plot_tangle_cell_crops.py) |
| **2G** | GO enrichment of genes rising along the trajectory | CosMx | [`plot_dotplot_pathway_linear_pseudotime_phf1.r`](R/plot_dotplot_pathway_linear_pseudotime_phf1.r)<br>(inputs from [`deg_dream_linear_pseudotime_phf1.r`](R/deg_dream_linear_pseudotime_phf1.r) and [`pathway_enrichment_cbln2_trajectory.R`](R/pathway_enrichment_cbln2_trajectory.R)) |
| **2H** | IMC state markers, tangle-bearing vs tangle-free (donor-interval forest) | IMC | [`imc_phf1_state_markers_forest_panel.R`](R/imc_phf1_state_markers_forest_panel.R)<br>[`imc_phf1_state_markers_subsets.R`](R/imc_phf1_state_markers_subsets.R) (the model estimates) |
| **2I** | Representative IMC layer 2/3 region (PHF1/p62/ubK48/DNA) | IMC image | [`imc_l23_roi_phf1_p62_ubk48.py`](python/imc_l23_roi_phf1_p62_ubk48.py) |
| **3A** | PHF1 intensity in tangle-free CBLN2 vs distance | CosMx | [`phf1_intensity_gradient_neurons.R`](R/phf1_intensity_gradient_neurons.R) |
| **3B** | Otero-Garcia module score over distance, PHF1-intensity-adjusted | CosMx | [`plot_otero_distance_adjustedonly_panel.R`](R/plot_otero_distance_adjustedonly_panel.R)<br>[`phf1_module_distance_intensity_adjust.R`](R/phf1_module_distance_intensity_adjust.R) (the fits) |
| **3C** | Pseudotime of tangle-free CBLN2 neurons vs distance, with hinge breakpoint | CosMx | [`pseudotime_by_phf1_distance.R`](R/pseudotime_by_phf1_distance.R) |
| **3D** | nDEG vs distance, per neuronal subtype | CosMx | [`nDEG_by_celltype_linear.R`](R/nDEG_by_celltype_linear.R) |
| **3E** | GO dotplot, distance DEG (proximity-signed) | CosMx | [`plot_dotplot_pathway_linear_distance_phf1.r`](R/plot_dotplot_pathway_linear_distance_phf1.r) |
| **3F** | Reactome 'cellular responses to stimuli' over distance, CBLN2 | CosMx | [`plot_reactome_stress_death_vs_phf1_distance.R`](R/plot_reactome_stress_death_vs_phf1_distance.R)<br>[`run_reactome_stress_death_distance.sh`](R/run_reactome_stress_death_distance.sh) |
| **3G** | Reactome 'autophagy' over distance, CBLN2 | CosMx | [`plot_reactome_stress_death_vs_phf1_distance.R`](R/plot_reactome_stress_death_vs_phf1_distance.R) |
| **3H** | IMC protein gradients vs tangle distance | IMC | [`imc_phf1_distance_markers_subsets.R`](R/imc_phf1_distance_markers_subsets.R) |
| **4A** | Cluster loss (AD vs control) vs distance to nearest PHF1high neuron | IMC | [`imc_cluster_loss_vs_tangle_distance.R`](R/imc_cluster_loss_vs_tangle_distance.R) |
| **4B** | Cluster loss vs PHF1 intensity | IMC | [`imc_cluster_loss_vs_phf1_intensity.R`](R/imc_cluster_loss_vs_phf1_intensity.R) |
| **4C** | Excitatory loss by cortical-depth decile vs PHF1 intensity | IMC | [`imc_tangle_depth_density_shift.R`](R/imc_tangle_depth_density_shift.R) |
| **4D** | Inhibitory loss by cortical-depth decile vs PHF1 intensity | IMC | [`imc_tangle_depth_density_shift.R`](R/imc_tangle_depth_density_shift.R) `inh` |
| **5A** | IMC nucleus area vs tangle distance | IMC | [`imc_phf1_morphology_dist_phf1_panel.R`](R/imc_phf1_morphology_dist_phf1_panel.R) |
| **5B** | IMC total nuclear DNA vs tangle distance | IMC | [`imc_phf1_dna_dist_phf1_panel.R`](R/imc_phf1_dna_dist_phf1_panel.R) |
| **5C** | CosMx rRNA channel vs tangle distance, CBLN2 | CosMx | [`phf1_channel_intensity_exc.R`](R/phf1_channel_intensity_exc.R)<br>[`run_morphology_observed.sh`](R/run_morphology_observed.sh)<br>[`nn3_phf1_vs_neg.R`](R/nn3_phf1_vs_neg.R) (covariate cache) |
| **5D** | Representative IMC nuclei, tangle-bearing vs tangle-free, paired by subcluster | IMC image | [`imc_phf1_nucleus_crops.py`](python/imc_phf1_nucleus_crops.py)<br>[`imc_phf1_nucleus_crops_figure.py`](python/imc_phf1_nucleus_crops_figure.py)<br>[`imc_phf1_nucleus_crops_pairs.py`](python/imc_phf1_nucleus_crops_pairs.py)<br>[`imc_phf1_nucleus_crops_singles.py`](python/imc_phf1_nucleus_crops_singles.py) |
| **6A** | Tangle-to-glia distance vs label-shuffled null | CosMx | [`plot_dist_to_phf1_by_glia_vs_null.R`](R/plot_dist_to_phf1_by_glia_vs_null.R)<br>[`null_glia_distance_within_celltype.R`](R/null_glia_distance_within_celltype.R)<br>[`plot_dist_to_phf1_braak_celltype.R`](R/plot_dist_to_phf1_braak_celltype.R) (per-cell distances) |
| **6B** | Astrocyte state module scores vs distance (Serrano-Pozo states) | CosMx | [`plot_serrano_astro_vs_phf1_distance_modelp.R`](R/plot_serrano_astro_vs_phf1_distance_modelp.R) |
| **6C** | GFAP vs tangle distance | IMC (+ CosMx) | [`imc_cosmx_gfap_distance_panel.R`](R/imc_cosmx_gfap_distance_panel.R)<br>[`imc_phf1_glia_distance.R`](R/imc_phf1_glia_distance.R) (IMC curve and fit) |
| **6D** | Representative IMC ROI (PHF1/GFAP/S100B/DNA) | IMC image | *not archived* |
| **6E** | Microglial reactive module scores vs distance | CosMx | [`plot_modulescore_vs_phf1_distance_modelp.R`](R/plot_modulescore_vs_phf1_distance_modelp.R) |

### Supplementary figures

| Panel | Panel shows | Data | Script(s) |
|---|---|---|---|
| **S1A** | Otero-Garcia distance slope before vs after PHF1-intensity adjustment | CosMx | [`replot_intensity_adjust_forest_otero.R`](R/replot_intensity_adjust_forest_otero.R)<br>[`phf1_module_distance_intensity_adjust.R`](R/phf1_module_distance_intensity_adjust.R) (the fits) |
| **S1B** | Intensity-adjusted slope after 0-200 µm exclusion | CosMx | [`phf1_module_distance_intensity_adjust.R`](R/phf1_module_distance_intensity_adjust.R) |
| **S1C** | Slope within PHF1-intensity deciles | CosMx | [`phf1_module_distance_intensity_adjust.R`](R/phf1_module_distance_intensity_adjust.R) |
| **S2A-C** | Stably expressed gene (SEG) negative control over distance | CosMx | [`plot_seg_control_composite_matched.R`](R/plot_seg_control_composite_matched.R) `--astro_source serrano --cbln2_series otero`<br>[`modulescore_seg_control_vs_phf1_distance.R`](R/modulescore_seg_control_vs_phf1_distance.R) (the fits)<br>[`prepare_seg_geneset_reference.R`](R/prepare_seg_geneset_reference.R) (the SEG reference) |
| **S3A** | Exclusion-radius concordance of distance coefficients | CosMx | [`collate_exclusion_radius.r`](R/collate_exclusion_radius.r)<br>[`run_deg_linear_distance.sh`](R/run_deg_linear_distance.sh) |
| **S3B** | Distance coefficients with scaled cortical depth added as a covariate | CosMx | [`collate_depth_covariate.r`](R/collate_depth_covariate.r)<br>[`deg_dream_linear_distance_depth_phf1.r`](R/deg_dream_linear_distance_depth_phf1.r) (sources [`cortical_depth_utils.R`](R/cortical_depth_utils.R)) |
| **S3C** | Distance coefficients after simulating unseen tangles in the third dimension | CosMx | [`collate_3d_distance_deg.r`](R/collate_3d_distance_deg.r)<br>[`deg_dream_linear_distance_3d_phf1.r`](R/deg_dream_linear_distance_3d_phf1.r)<br>[`run_deg_linear_distance_3d.sh`](R/run_deg_linear_distance_3d.sh) |
| **S4A** | Reactome programmed cell death over distance, CBLN2 | CosMx | [`plot_reactome_stress_death_vs_phf1_distance.R`](R/plot_reactome_stress_death_vs_phf1_distance.R) |
| **S4B** | Per-gene decomposition, intrinsic apoptosis - tangle contrast | CosMx | [`plot_reactome_gene_contribution.R`](R/plot_reactome_gene_contribution.R)<br>[`run_reactome_gene_contribution.sh`](R/run_reactome_gene_contribution.sh) |
| **S4C** | Per-gene decomposition, intrinsic apoptosis - proximity | CosMx | [`plot_reactome_gene_contribution.R`](R/plot_reactome_gene_contribution.R) |
| **S5A** | IMC nucleus eccentricity vs tangle distance | IMC | [`imc_phf1_morphology_dist_phf1_panel.R`](R/imc_phf1_morphology_dist_phf1_panel.R) |
| **S5B** | CosMx nucleus area vs tangle distance, CBLN2 | CosMx | [`phf1_morphology_dist_phf1_panel.R`](R/phf1_morphology_dist_phf1_panel.R)<br>[`run_morphology_observed.sh`](R/run_morphology_observed.sh), which also runs [`phf1_morphology_exc.R`](R/phf1_morphology_exc.R) and [`phf1_distance_morphology_exc.R`](R/phf1_distance_morphology_exc.R) |
| **S6A** | Covariate ladder: IMC GFAP vs distance | IMC | [`imc_covariate_sensitivity.R`](R/imc_covariate_sensitivity.R) |
| **S6B** | Covariate ladder: IMC state markers, tangle contrast | IMC | [`imc_covariate_sensitivity.R`](R/imc_covariate_sensitivity.R)<br>[`imc_phf1_state_markers_subsets.R`](R/imc_phf1_state_markers_subsets.R) (reference fits) |
| **S6C** | Covariate ladder: IMC state markers vs distance | IMC | [`imc_covariate_sensitivity.R`](R/imc_covariate_sensitivity.R) |
| **S6D** | IMC state markers vs tangle distance, stratified by plaque proximity | IMC | [`imc_amyloid_tangle_coupling.R`](R/imc_amyloid_tangle_coupling.R) |
| **S6E** | IMC GFAP in astrocytes vs tangle distance, stratified by plaque proximity | IMC | [`imc_amyloid_tangle_coupling.R`](R/imc_amyloid_tangle_coupling.R) |
| **S7A-C** | Reactome cell death / stress / autophagy modules over distance, Astro | CosMx | [`plot_reactome_stress_death_vs_phf1_distance.R`](R/plot_reactome_stress_death_vs_phf1_distance.R) |
| **S7D-F** | Reactome cell death / stress / autophagy modules over distance, Micro | CosMx | [`plot_reactome_stress_death_vs_phf1_distance.R`](R/plot_reactome_stress_death_vs_phf1_distance.R) |

### Crop panels (2E, 2F)

Both are rendered by `plot_tangle_cell_crops.py` with `--celltypes Exc-IT-L2-3-CBLN2-HOPX
--phf1-status pos --window-um 60 --size-cm 4 --scalebar-um 10`, highlighting:

- **2E**: the Otero-Garcia Exc1/2 UP panel genes (`genesets/fig2_contrast/otero_l23_up.txt`,
  319 genes) via `--highlight-genes`;
- **2F**: `--highlight-genes "MAP1B;STMN1;STMN2;STMN3;STMN4;MAPT"`, the microtubule
  (de)polymerisation genes of the Fig. 2G term.

The figure cells were chosen by eye from the per-donor galleries and re-rendered at 600 dpi
with `--only-cells`.

## Supplementary tables → scripts

| Table | Contents | Script(s) |
|---|---|---|
| **S1** | Donor characteristics | *cohort metadata, not script-generated* |
| **S2** | Marker gene sets for the two rounds of AUCell marker-set annotation ("label-transfer marker gene sets" in the manuscript) | [`build_aucell_annotation_tables.R`](R/build_aucell_annotation_tables.R), which also writes the set sizes, thresholds, assignment parameters and per-round counts quoted in the Methods |
| **S3** | Cell quality control by donor | *counts of the objects written by [`add_metadata.R`](R/add_metadata.R), [`QC.R`](R/QC.R) and [`doublet_detection_scDblFinder.R`](R/doublet_detection_scDblFinder.R); no separate script* |
| **S4** | Genes vs pseudotime along the tangle-associated trajectory (CBLN2) | [`deg_dream_linear_pseudotime_phf1.r`](R/deg_dream_linear_pseudotime_phf1.r) (step 5 of [`run_traj_pipeline.sh`](R/run_traj_pipeline.sh)) |
| **S5** | GO overrepresentation for S4 | [`pathway_enrichment_cbln2_trajectory.R`](R/pathway_enrichment_cbln2_trajectory.R) |
| **S6** | Genes vs distance to nearest tangle-bearing neuron (CBLN2) | [`deg_dream_linear_distance_phf1.r`](R/deg_dream_linear_distance_phf1.r) |
| **S7** | GO overrepresentation for S6 | [`pathway_enrichment_enrichR_linear.r`](R/pathway_enrichment_enrichR_linear.r) |
| **S8** | Exclusion-radius sensitivity | [`collate_exclusion_radius.r`](R/collate_exclusion_radius.r), [`run_deg_linear_distance.sh`](R/run_deg_linear_distance.sh) |
| **S9** | Phospho-tau marker genes and the Otero-Garcia layer 2/3 signature | [`find_phf1_markers_by_celltype.R`](R/find_phf1_markers_by_celltype.R), [`otero_geneset.R`](R/otero_geneset.R) |

## Running order

The figure scripts all read objects built by the upstream pipeline, so that has to run
first. Within a stage, order is as listed.

1. **From the AtoMx export to annotated cell types** —
   [`merge_seu.R`](R/merge_seu.R) (sources
   [`assign_sample_id_from_fov.R`](R/assign_sample_id_from_fov.R)) →
   [`add_metadata.R`](R/add_metadata.R) →
   [`QC.R`](R/QC.R) →
   [`doublet_detection_scDblFinder.R`](R/doublet_detection_scDblFinder.R) →
   [`SCTransform.R`](R/SCTransform.R) →
   [`AUCell_1.R`](R/AUCell_1.R) →
   [`AUCell_Neuron.R`](R/AUCell_Neuron.R) →
   [`reintegrate_seu.R`](R/reintegrate_seu.R)

   Cell types come from **AUCell marker-set scoring against an snRNA-seq entorhinal
   cortex reference**, in two rounds: broad classes, then neuronal subtypes among the
   cells called neuronal. Figures 1B and 1C are the marker heatmaps from these rounds.

   Unsupervised clustering and Seurat anchor-based label transfer were both explored but
   **not used**, so those scripts are not published here. `AUCell_1.R` reads
   `seu_transformed.RDS` directly. The original run read a later object that also carried
   those unused annotations; AUCell_1.R uses none of them, so the input change does not
   alter its result.

2. **Objects and the distance field** — [`add_PHF1.R`](R/add_PHF1.R) →
   [`rename_celltypes.R`](R/rename_celltypes.R) →
   [`export_seu_coords_for_phf1.R`](R/export_seu_coords_for_phf1.R) →
   [`extract_phf1_intensity.py`](python/extract_phf1_intensity.py) →
   [`add_phf1_intensity.R`](R/add_phf1_intensity.R) →
   [`add_phf1_threshold_calls.R`](R/add_phf1_threshold_calls.R) →
   [`splitting_sce_cluster_celltype.r`](R/splitting_sce_cluster_celltype.r) →
   [`run_label_phf1_neighbours.sh`](R/run_label_phf1_neighbours.sh) →
   [`add_dist_to_phf1_seu.R`](R/add_dist_to_phf1_seu.R) →
   [`add_dist_to_full_sce.R`](R/add_dist_to_full_sce.R) →
   [`build_group_sce_neighbours.r`](R/build_group_sce_neighbours.r)

   `add_PHF1.R` is where the analysis is defined: it sets the `PHF1` boolean from the
   manual tangle annotation. Every distance in the paper is a nearest-neighbour distance
   to that anchor set, so changing it changes `dist_to_phf1_um` for **every** cell, not
   just those near an added tangle.

   `export_seu_coords_for_phf1.R` writes `PHF1/seu_coords.csv`, the per-cell coordinate
   table. The Python extractor reads it to measure each cell's PHF1 mask intensity, and
   `add_phf1_intensity.R` joins the result onto the objects. The same table is read by
   the image scripts behind 1E, 2E and 2F, by the S3B/S3C sensitivity arms and by the
   Table S2/S3 builder.

   Do not skip `rename_celltypes.R`: the upstream object carries short cell-type labels
   (`Exc-IT-L2-3`) while every downstream script expects the long ones
   (`Exc-IT-L2-3-CBLN2-HOPX`), and the neighbour step would otherwise find zero anchors.

3. **Derived sets** — [`find_phf1_markers_by_celltype.R`](R/find_phf1_markers_by_celltype.R)
   (the PHF1 marker set, Table S9),
   [`otero_geneset.R`](R/otero_geneset.R) then
   [`prepare_fig2_contrast_genesets.R`](R/prepare_fig2_contrast_genesets.R) (the Otero-Garcia
   layer 2/3 set as panel probes),
   [`run_generate_phf1_null_labels.sh`](R/run_generate_phf1_null_labels.sh) (the 1,000
   label permutations used by Figure 6A) and
   [`prepare_seg_geneset_reference.R`](R/prepare_seg_geneset_reference.R) (the stably
   expressed gene reference for Fig. S2).

4. **Tangle-associated trajectory (Exc-IT-L2-3-CBLN2-HOPX)** —
   `bash R/run_traj_pipeline.sh Exc-IT-L2-3-CBLN2-HOPX <cores>` (the cell type must be given:
   the runner's default is another subtype) runs, in order,
   [`plot_cbln2_umap_phf1_p95.R`](R/plot_cbln2_umap_phf1_p95.R) (embedding, cached via
   [`cbln2_embedding_utils.R`](R/cbln2_embedding_utils.R)),
   [`cbln2_trajectory_optimise.R`](R/cbln2_trajectory_optimise.R) (settings search),
   [`cbln2_slingshot_trajectory.R`](R/cbln2_slingshot_trajectory.R) (fit; writes the per-cell
   pseudotime table),
   [`plot_cbln2_traj_phf1_by_bin_mean.R`](R/plot_cbln2_traj_phf1_by_bin_mean.R),
   [`deg_dream_linear_pseudotime_phf1.r`](R/deg_dream_linear_pseudotime_phf1.r) (Table S4),
   [`pathway_enrichment_cbln2_trajectory.R`](R/pathway_enrichment_cbln2_trajectory.R) (Table S5)
   and [`cbln2_trajectory_selection_null.R`](R/cbln2_trajectory_selection_null.R) (the
   permutation p of Fig. 2B). Then
   [`plot_traj_umap_publication.R`](R/plot_traj_umap_publication.R),
   [`plot_cbln2_traj_otero_by_bin.R`](R/plot_cbln2_traj_otero_by_bin.R),
   [`plot_dotplot_pathway_linear_pseudotime_phf1.r`](R/plot_dotplot_pathway_linear_pseudotime_phf1.r)
   and, after stage 2, [`pseudotime_by_phf1_distance.R`](R/pseudotime_by_phf1_distance.R).

5. **Differential expression over distance** —
   [`deg_dream_linear_distance_phf1.r`](R/deg_dream_linear_distance_phf1.r)
   (cell-level, vs distance), plus the sensitivity arms
   [`run_deg_linear_distance.sh`](R/run_deg_linear_distance.sh) (exclusion radii),
   [`deg_dream_linear_distance_depth_phf1.r`](R/deg_dream_linear_distance_depth_phf1.r)
   and [`run_deg_linear_distance_3d.sh`](R/run_deg_linear_distance_3d.sh).

6. **Enrichment** — [`pathway_enrichment_enrichR_linear.r`](R/pathway_enrichment_enrichR_linear.r)
   (Table S7). The enrichment steps call the enrichR web API and need outbound internet.

7. **Figures** — everything in the panel tables above. The module-score scripts load the
   full Seurat object and are the memory-hungry step. Some panels read tables that
   another script writes, so run that script first:
   - **1F**: [`plot_phf1_intensity_validation.R`](R/plot_phf1_intensity_validation.R)
     before `plot_phf1_intensity_validation_neurons.R`. It writes the per-cell p95 table in
     `plots/phf1_intensity_validation/`; the 1F script restricts it to neurons.
   - **3B, S1**: [`phf1_module_distance_intensity_adjust.R`](R/phf1_module_distance_intensity_adjust.R)
     before `plot_otero_distance_adjustedonly_panel.R` and
     `replot_intensity_adjust_forest_otero.R`, which redraw its source data without refitting.
   - **S2**: [`run_phf1_module_distance_modelp.sh`](R/run_phf1_module_distance_modelp.sh),
     [`run_modulescore_distance_modelp.sh`](R/run_modulescore_distance_modelp.sh),
     `plot_serrano_astro_vs_phf1_distance_modelp.R` and
     `phf1_module_distance_intensity_adjust.R` before
     `modulescore_seg_control_vs_phf1_distance.R` and `plot_seg_control_composite_matched.R`.
   - **5C, S5B**: [`nn3_phf1_vs_neg.R`](R/nn3_phf1_vs_neg.R) before
     [`run_morphology_observed.sh`](R/run_morphology_observed.sh). It writes the per-neuron
     covariate cache (`plots/nn3_neuron_spacing/source_data_nn3_cells.tsv`) that the
     morphology and channel-intensity scripts read.
   - **6A**: [`plot_dist_to_phf1_braak_celltype.R`](R/plot_dist_to_phf1_braak_celltype.R)
     before `plot_dist_to_phf1_by_glia_vs_null.R`. It writes the per-cell distances,
     glia summary and leave-one-out contrasts in `plots/dist_to_phf1/`.
   - **2H, S6B**: [`imc_phf1_state_markers_subsets.R`](R/imc_phf1_state_markers_subsets.R)
     before `imc_phf1_state_markers_forest_panel.R` and `imc_covariate_sensitivity.R`.
     It writes the per-marker estimates in `plots/imc_phf1_state_markers_subsets/`.
   - **6C**: [`imc_phf1_glia_distance.R`](R/imc_phf1_glia_distance.R) before
     `imc_cosmx_gfap_distance_panel.R`. It writes the IMC rolling mean and fit in
     `plots/imc_phf1_glia_distance/`.

8. **Supplementary Table S2** —
   [`build_aucell_annotation_tables.R`](R/build_aucell_annotation_tables.R), after the
   AUCell rounds and `export_seu_coords_for_phf1.R`.

## Inputs these scripts do not create

Besides the AtoMx export and the IMC object, these inputs are curated by hand or come
from elsewhere:

- `PHF1/PHF1_Score.xlsx`: the manual tangle annotation (one `cell_id` per PHF1+ cell),
  read by `add_PHF1.R`, plus the PHF1 images and the per-sample alignment matrices in
  `PHF1/Transformed/` used by the image scripts.
- `CosMx_Samples.xlsx`: donor metadata (Braak stage, sex, age, PMI), read by
  `add_metadata.R`.
- The snRNA-seq entorhinal cortex reference objects that the AUCell marker sets are
  derived from (`AUCell_1.R`, `AUCell_Neuron.R`).
- The MAP snRNA-seq reference (`sce_ref.qs`) that `prepare_seg_geneset_reference.R` fits
  the stably expressed gene index on.
- Published gene sets (Otero-Garcia, Serrano-Pozo, Cameron, Mancuso, Pandey). The
  Otero-Garcia supplementary workbooks are read by `otero_geneset.R` (Table S6 of the source
  publication) and `prepare_fig2_contrast_genesets.R` (Tables S2 and S5, as
  `OteroGarcia_S2.xlsx` / `OteroGarcia_S5.xlsx` in `<PROJECT_ROOT>`).
- Reactome: the Enrichr `Reactome_Pathways_2024` library, cached as
  `genesets/enrichr_gmt/Reactome_Pathways_2024.gmt` (download from
  `https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_Pathways_2024`;
  the published run used the copy fetched 2026-08-10, md5
  `c30ddaf5ce1a1e8fd7f7850b39ebe4f6`), plus the Reactome hierarchy files, which
  `reactome_pcd_family.R` downloads itself.
- IMC, Fig. 4: `<IMC_ROOT>/Figures/SupplementaryTable4_patient_percent_change.csv`, the
  per-donor cluster changes published with the IMC dataset. The two cluster-loss scripts
  recompute the same quantity from `spe.rds` and check it against this table.
- IMC, Fig. 2I: `<IMC_ROOT>/imc_cell_table.csv`, a CSV export of `colData(spe)` with at least
  `sample_id`, `ObjectNumber`, `layer` and `PHF1_Otsu`, and the ROI image and
  `regionprops/` files under `<USER_DATA_ROOT>/EC_IMC_Project/`.

## GEO flat files

AtoMx exported only the transcript, metadata and polygon files for each slide. The per-slide
`exprMat_file` and `fov_positions_file` deposited on GEO were rebuilt from the transcript file
by [`reconstruct_atomx_flatfiles.R`](R/reconstruct_atomx_flatfiles.R), which checks the rebuilt
counts against the AtoMx Seurat export (identical for every cell and target) and the FOV
positions against the instrument run summary, and stops on any disagreement.

## Repository layout

```
R/                    119 R and shell scripts
python/                9 python scripts (PHF1 intensity extraction, CosMx overlays
                        and crops, IMC images and nucleus crops)
docs/MODELS.md          model specifications - which model backs which panel
figure_panel_map.tsv    panel -> script map, machine-readable
PLACEHOLDERS.md         every path token and what to set it to
CITATION.cff            citation metadata for this code
LICENSE                 MIT
```

## Additional points

- **`percent_neg` is in the cell-level transcript models only.** It is absent from every
  model whose outcome is an image measure.
- **enrichR output folders are named by logFC sign.** In the distance analysis `DOWN/` is
  the *up-near-tangles* direction; in the pseudotime analysis `UP/` is the direction rising
  towards the tangle-rich end of the trajectory.
- **Plots are written as a triple**: `plot_<name>.pdf`, `source_data_<name>.tsv` (the exact
  rows drawn) and `stats_<name>.txt` (model, summary, effect size, `sessionInfo()`).

## Software versions

As recorded in the manuscript Key Resources Table, plus the packages used by the trajectory
analysis. These are the versions the analyses **ran under**, which are not necessarily what is
installed today — the `sessionInfo()` block at the foot of every stats log in the data release
is the primary record.

| Package | Version | | Package | Version |
|---|---|---|---|---|
| R (CosMx) | 4.4.3 | | enrichR | 3.4 |
| R (IMC) | 4.4.2 | | AUCell | 1.28.0 |
| Seurat | 5.3.1 | | spatstat.geom | 3.6-0 |
| SeuratObject | 5.4.0 | | spatstat.explore | 3.5-3 |
| sctransform | 0.4.2 | | spatstat.random | 3.4-2 |
| SingleCellExperiment | 1.28.1 | | ComplexHeatmap | 2.22.0 |
| SummarizedExperiment | 1.36.0 | | ggplot2 | 4.0.3 |
| limma | 3.62.2 | | dplyr | 1.2.1 |
| statmod | 1.5.2 | | qs | 0.27.3 |
| edgeR | 4.4.2 | | Matrix | 1.7-6 |
| variancePartition | 1.36.2 | | argparse | 2.3.1 |
| lme4 | 1.1-38 | | zoo | 1.9-0 |
| lmerTest | 3.2-0 | | glmGamPoi | 1.18.0 |
| emmeans | 2.0.3 | | scDblFinder | 1.23.4 |
| multcomp | 1.4-29 | | scMerge | 1.22.0 |
| mgcv | 1.9-4 | | RANN | 2.6.2 |
| slingshot | 2.14.0 | | TrajectoryUtils | 1.14.0 |
| princurve | 2.1.6 | | harmony | see stats log |

Imaging and instrument software:

| Software | Version |
|---|---|
| AtoMx Spatial Informatics Platform | 1.3.2 |
| napari | 0.4.19.post1 |
| napari-cosmx plugin | 0.4.17.3 |
| affinder | 0.5.0 |
| LAS X | 4.7.0.28176 |

Python 3.11 for the imaging scripts, with numpy 2.2.0, tifffile 2025.10.16 and
scikit-image 0.25.2.

Two version splits are deliberate and worth knowing if you compare logs: the CosMx and
IMC analyses ran on different stacks (R 4.4.3 vs 4.4.2), and the later sensitivity
analyses — the plaque stratification (Figs S6D/S6E) and the 3D tangle simulation
(Fig S3C) — ran on lme4 2.0-6 / lmerTest 3.2-1 / RANN 2.6.3. The neuron-only
Fig. 1F (`plot_phf1_intensity_validation_neurons.R`; dplyr/ggplot2 only, no models) and the
IMC neuronal-loss panels (Fig. 4) ran locally on R 4.4.2.

## Citation

Please cite the manuscript when using this code. Citation metadata for the code itself
is in [CITATION.cff](CITATION.cff).
