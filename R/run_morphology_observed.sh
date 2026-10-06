#!/bin/bash
# ---------------------------------------------------------------------------
# run_morphology_observed.sh
#
# Figure panels: 5C, S7B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_morphology_observed.sh
#
# All OBSERVED analyses of the CosMx morphology / channel-intensity family, in one job.
# Fast (~10 min total).
#
#   qsub R/run_morphology_observed.sh
#
# Produces, per distance cap (300 and 1000 um):
#   plots/phf1_morphology/                            PHF1+ vs PHF1- (cap-independent)
#   plots/phf1_distance_morphology_<cap>um/           morphology vs distance
#   plots/phf1_morphology_dist_phf1_panel_<cap>um/    the combined one-panel figures
#   plots/phf1_channel_intensity/                     DAPI/Histone/rRNA (cap-independent)
#   plots/phf1_channel_intensity_distance_<cap>um/
#   plots/phf1_channel_intensity_combined_<cap>um/

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -N morphology_observed
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

cd $PBS_O_WORKDIR

PROJ=<PROJECT_ROOT>/phf1_v2
log_dir="$PROJ/R/logs/morphology_observed"
mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START"

# --- cell-autonomous halves: cap-independent, run once ----------------------
Rscript "$PROJ/R/phf1_morphology_exc.R" \
  > "$log_dir/phf1_morphology_exc.log" 2>&1

# --- everything with a distance cap ----------------------------------------
for CAP_UM in 1000 300; do
  Rscript "$PROJ/R/phf1_distance_morphology_exc.R" ${CAP_UM} \
    > "$log_dir/phf1_distance_morphology_exc_${CAP_UM}um.log" 2>&1

  Rscript "$PROJ/R/phf1_morphology_dist_phf1_panel.R" ${CAP_UM} \
    > "$log_dir/phf1_morphology_dist_phf1_panel_${CAP_UM}um.log" 2>&1

  Rscript "$PROJ/R/phf1_channel_intensity_exc.R" ${CAP_UM} \
    > "$log_dir/phf1_channel_intensity_exc_${CAP_UM}um.log" 2>&1
done

END=$(date)
echo "Job ended at $END"
echo "Logs in $log_dir"
