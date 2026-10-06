#!/bin/bash
# ---------------------------------------------------------------------------
# run_modulescore_distance_modelp.sh
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
# run_modulescore_distance_modelp.sh
#
# Runs plot_modulescore_vs_phf1_distance_modelp.R: significance comes from the model's
# BH-adjusted p (lmm_padj, from the log-distance LMM). Coloured-line figures with
# 95% CI shading.
#
# Distance cap is set by CAP_UM below. The output/log dirs are cap-specific
# (..._<CAP_UM>um), so changing the cap NEVER overwrites a previous run.
#
#   qsub run_modulescore_distance_modelp.sh

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=120gb
#PBS -N modulescore_distance_modelp
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

cd $PBS_O_WORKDIR

# ---- distance cap (um): the ONLY input to change ----
CAP_UM=1000

resource_dir=<PROJECT_ROOT>/phf1_v2/R
out_dir=plots/module_score_vs_phf1_distance_modelp_${CAP_UM}um
log_dir=<PROJECT_ROOT>/phf1_v2/R/logs/module_score_vs_phf1_distance_modelp_${CAP_UM}um

mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START (cap ${CAP_UM} um)"

Rscript $resource_dir/plot_modulescore_vs_phf1_distance_modelp.R \
  --seu          seu_PHF1.rds \
  --geneset_dir  <EXTERNAL_GENESET_DIR> \
  --output_dir   "$out_dir" \
  --max_dist_um  ${CAP_UM} \
  --seed         42 \
  > "$log_dir/module_score_vs_phf1_distance_modelp_${CAP_UM}um.log" 2>&1

END=$(date)
echo "Job ended at $END" >> "$log_dir/module_score_vs_phf1_distance_modelp_${CAP_UM}um.log"
