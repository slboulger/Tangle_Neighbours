#!/bin/bash
# ---------------------------------------------------------------------------
# run_phf1_module_distance_modelp.sh
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
# run_phf1_module_distance_modelp.sh
#
# Runs plot_phf1_module_vs_phf1_distance_modelp.R: scores THREE ptau signatures per
# eligible neuron subtype -- the CELLTYPE-SPECIFIC PHF1 marker geneset (from
# find_phf1_markers_by_celltype.R) and the Otero-Garcia AT8 UP / DOWN sets (L2-3 by
# default) -- over distance to the nearest PHF1+ neuron, with the same model-FDR
# machinery as run_modulescore_distance_modelp.sh (log-distance LMM, plus the decay
# fit). Each signature gets its individual plots (fitted log-distance curve + actual-data
# raw, PHF1+ reference band on each), plus a per-celltype INTEGRATED plot overlaying the three z-scored
# rolling-mean trends with their PHF1+ references.
#
# Distance cap is set by CAP_UM below; output dir is cap-specific (..._<CAP_UM>um).
#
#   qsub run_phf1_module_distance_modelp.sh

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=120gb
#PBS -N phf1_module_distance_modelp
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

cd $PBS_O_WORKDIR

# ---- distance cap (um): the ONLY input to change ----
CAP_UM=1000

resource_dir=<PROJECT_ROOT>/phf1_v2/R
out_dir=plots/phf1_module_vs_phf1_distance_modelp_${CAP_UM}um
log_dir=<PROJECT_ROOT>/phf1_v2/R/logs/phf1_module_vs_phf1_distance_modelp_${CAP_UM}um

mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START (cap ${CAP_UM} um)"

Rscript $resource_dir/plot_phf1_module_vs_phf1_distance_modelp.R \
  --seu            seu_PHF1.rds \
  --markers_dir    phf1_markers \
  --otero_rds      otero_signatures/otero_at8_signatures.rds \
  --otero_up_set   otero_L23_up \
  --otero_down_set otero_L23_down \
  --output_dir     "$out_dir" \
  --max_dist_um    ${CAP_UM} \
  --window_um      75 \
  --seed           42 \
  > "$log_dir/phf1_module_vs_phf1_distance_modelp_${CAP_UM}um.log" 2>&1

END=$(date)
echo "Job ended at $END" >> "$log_dir/phf1_module_vs_phf1_distance_modelp_${CAP_UM}um.log"
