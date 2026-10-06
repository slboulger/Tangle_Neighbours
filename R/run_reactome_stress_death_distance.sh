#!/bin/bash
# ---------------------------------------------------------------------------
# run_reactome_stress_death_distance.sh
#
# Figure panels: 3C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_reactome_stress_death_distance.sh
#
# Curated Reactome DEATH-COMMITMENT, STRESS-ADAPTATION and PROTEOSTASIS module scores vs
# distance to the nearest PHF1+ (ptau) neuron, in Exc-IT-L2-3-CBLN2-HOPX, Astro and Micro.
# Runs plot_reactome_stress_death_vs_phf1_distance.R: AddModuleScore (SCT) per celltype,
# canonical log-distance LMM, model-FDR significance (BH within celltype x group), one
# rolling-mean overlay figure per celltype x group (9 total) with a PHF1+ mean/95% CI
# marginal strip at the left edge of the neuron panels (plus a _nostrip copy of each),
# and PHF1+ vs PHF1- forest + paired-donor figures for the neuron celltype.
#
# Distance cap is set by CAP_UM below. The output/log dirs are cap-specific
# (..._<CAP_UM>um), so changing the cap NEVER overwrites a previous run.
#
# Runs either way:
#   qsub run_reactome_stress_death_distance.sh     # batch
#   bash R/run_reactome_stress_death_distance.sh   # interactive session (preferred; ~64 GB is enough,
#                                                  # the Seurat object is freed right after scoring)

#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=4:mem=120gb
#PBS -N reactome_stress_death_distance
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

# PBS_O_WORKDIR is only set under qsub; fall back to the project root when run directly.
project_root=<PROJECT_ROOT>/phf1_v2
cd "${PBS_O_WORKDIR:-$project_root}"

# ---- distance cap (um): the ONLY input to change ----
CAP_UM=1000

resource_dir=$project_root/R
out_dir=plots/reactome_stress_death_vs_phf1_distance_${CAP_UM}um
log_dir=$project_root/R/logs/reactome_stress_death_vs_phf1_distance_${CAP_UM}um

mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START (cap ${CAP_UM} um)"

Rscript $resource_dir/plot_reactome_stress_death_vs_phf1_distance.R \
  --seu          seu_PHF1.rds \
  --output_dir   "$out_dir" \
  --max_dist_um  ${CAP_UM} \
  --window_um    75 \
  --seed         42 \
  --overwrite    yes \
  > "$log_dir/reactome_stress_death_vs_phf1_distance_${CAP_UM}um.log" 2>&1

END=$(date)
echo "Job ended at $END" >> "$log_dir/reactome_stress_death_vs_phf1_distance_${CAP_UM}um.log"
