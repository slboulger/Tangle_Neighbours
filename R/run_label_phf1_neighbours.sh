#!/bin/bash
# ---------------------------------------------------------------------------
# run_label_phf1_neighbours.sh
#
# Upstream pipeline - builds the objects every panel reads
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_label_phf1_neighbours.sh
#
# ONE-TIME job: labels all cells with near/distal proximity to PHF1+ neurons.
# Run this ONCE before the DEG array job below.

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=120gb
#PBS -N label_phf1_neighbours
#PBS -o <PROJECT_ROOT>/phf1_v2/R/logs/label_neighbours.o
#PBS -e <PROJECT_ROOT>/phf1_v2/R/logs/label_neighbours.e

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

resource_dir=<PROJECT_ROOT>/phf1_v2/R
data_dir=<PROJECT_ROOT>/phf1_v2
output_dir=${data_dir}/celltype_sce_neighbours

START=$(date)
echo "Job started at $START"

Rscript $resource_dir/label_phf1_neighbours.r \
  --all_sce          ${data_dir}/sce.qs \
  --mode             single \
  --celltype_sce_dir ${data_dir}/celltype_sce \
  --output_dir       ${output_dir} \
  --neuron_celltypes "Exc-IT-L2-3-CBLN2-HOPX,Exc-IT-L3-5-CHGA-IL1RAPL2,Exc-ET-L5-SPON1-FGD4,Exc-CT-L6-SYNPO2-SEMA3E,Exc-IT-L6-CTXN1-ERC2,Inh-SST,Inh-PVALB,Inh-VIP,Inh-LAMP5,Unassigned Neuron" \
  --distance_um      100 \
  --coord_x          x_slide_mm \
  --coord_y          y_slide_mm \
  --sample_col       sample_id

END=$(date)
echo "Job ended at $END"
