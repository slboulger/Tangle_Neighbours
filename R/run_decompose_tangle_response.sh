#!/bin/bash
# ---------------------------------------------------------------------------
# run_decompose_tangle_response.sh
#
# Figure panels: 4C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_decompose_tangle_response.sh
#
# Single PBS job. Runs decompose_tangle_response_phf1.r, which loops internally
# over the neuronal panel (focus + comparison subtypes). Needs internet for the
# enrichR GO calls. Log written to logs/decompose_tangle_response.log.

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -N decompose_tangle
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

cd $PBS_O_WORKDIR

resource_dir=<PROJECT_ROOT>/phf1_v2/R
log_dir=<PROJECT_ROOT>/phf1_v2/R/logs
mkdir -p "$log_dir"

echo "Job started at $(date)"
Rscript "$resource_dir/decompose_tangle_response_phf1.r" \
  > "$log_dir/decompose_tangle_response.log" 2>&1
echo "Job ended at $(date)"
