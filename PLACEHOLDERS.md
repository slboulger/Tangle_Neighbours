# Placeholders

Every absolute path has been removed from these scripts and replaced with an
angle-bracket token. **No script will run until you substitute them.**

| Placeholder | What to set it to |
|---|---|
| `<PROJECT_ROOT>` | The analysis project directory — the folder holding `R/`, `plots/`, `deg/`, `phf1_v2/` and the Seurat/SCE objects. Most scripts call `setwd("<PROJECT_ROOT>/phf1_v2")`. |
| `<RDS_ROOT>` | The storage root above the project (institutional RDS mount or its local equivalent). |
| `<USER_DATA_ROOT>` | The per-user data directory under `<RDS_ROOT>` that contains the CosMx run folder and `EC_IMC_Project/`, the IMC image folder (`img/`, `DNA_Otsu_Masks/`, `regionprops/`) read by the Fig. 5D nucleus-crop scripts. |
| `<IMC_ROOT>` | The directory holding the IMC SpatialExperiment `spe.rds`, read by every IMC script (Figs 2E, 3E, 5A, 5B, 5D, 6C, S7A, S8). See the README for where the object comes from. |
| `<EXTERNAL_GENESET_DIR>` | The shared snRNA-seq gene-set directory used by the glial module-score scripts. Read-only in the original environment. |
| `<CONDA_ROOT>` | Conda/miniforge installation root; the analysis env is `dgenv`, the imaging env `spatial_env`. |

To substitute them all at once, export the ones your copy needs and run the `perl` line
(it works with both GNU and BSD/macOS tools). A token you leave unset is left in place.

```bash
export PROJECT_ROOT=/path/to/your/project IMC_ROOT=/path/to/imc   # ...and any others
grep -rlE '<(PROJECT_ROOT|RDS_ROOT|USER_DATA_ROOT|EXTERNAL_GENESET_DIR|CONDA_ROOT|IMC_ROOT)>' R python |
  xargs perl -pi -e 's{<(PROJECT_ROOT|RDS_ROOT|USER_DATA_ROOT|EXTERNAL_GENESET_DIR|CONDA_ROOT|IMC_ROOT)>}{$ENV{$1} // "<$1>"}ge'
```

Check what is still unset with:

```bash
grep -rhoE '<(PROJECT_ROOT|RDS_ROOT|USER_DATA_ROOT|EXTERNAL_GENESET_DIR|CONDA_ROOT|IMC_ROOT)>' R python | sort | uniq -c
```

Other angle-bracket strings in the scripts (`<CT>`, `<CAP_UM>`, `<K>` ...) are templates
in comments and file-name descriptions, not paths.

## Donor identifiers

Donor IDs throughout are **UK Brain Banks Network (BBN) IDs**. Brain-bank case IDs have been
removed and do not appear anywhere in this repository.

A tenth CosMx ID appears only in `R/merge_seu.R` and `R/add_PHF1.R`. It is
the **excluded** donor: ten sections were run, and this one was found to be white matter
rather than entorhinal cortex. Its FOVs (slide 2, FOVs 42-100) are labelled in
`merge_seu.R` because they are in the AtoMx export, and `add_PHF1.R` removes the donor
(it has no PHF1+ cells). It is in no analysed object or result.

IMC donors and ROIs use the IDs of the public de-identified IMC object: the donor is its
`BBN_ID` and an ROI is `<BBN_ID>_<n>`. The Fig. 5D crop scripts
select ROIs by that name, so the image and mask files they read must be named the same
way.

One consequence worth knowing: `R/nn3_utils.r` carries a lookup table
(`l1_dir`) keyed by donor ID, used by the depth-adjusted sensitivity variants. Its
keys are BBN IDs, so the `sample_id` column of any object you run it against must
also be BBN-relabelled. Relabel at export time; do not edit working objects.
