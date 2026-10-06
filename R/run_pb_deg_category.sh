#!/bin/bash
# ---------------------------------------------------------------------------
# run_pb_deg_category.sh
#
# Produces: Table S4
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=4:mem=100gb
#PBS -N DEG_category
#PBS -J 1-17
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

cd $PBS_O_WORKDIR

resource_dir=<PROJECT_ROOT>/phf1_v2/R
input_dir=<PROJECT_ROOT>/phf1_v2/celltype_sce
output_dir=<PROJECT_ROOT>/phf1_v2/deg/de_cat
log_dir=<PROJECT_ROOT>/phf1_v2/R/logs/de_cat

# Quote everything: celltype names can contain spaces (e.g. "Unassigned Neuron").
input_sce=$(ls "$input_dir" | head -n $PBS_ARRAY_INDEX | tail -n 1)
celltype=$(basename "$input_sce" _sce.qs)

mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START"
echo "Input SCE: $input_sce"
echo "Celltype: $celltype"

Rscript $resource_dir/deg_pb_category_phf1.r \
  --sce              "$input_dir/$input_sce" \
  --dependent_var    PHF1 \
  --ref_class        FALSE \
  --confounding_vars Sex,Age,PMI \
  --output_dir       "$output_dir" \
  > "$log_dir/${celltype}.log" 2>&1

END=$(date)
echo "Job ended at $END" >> "$log_dir/${celltype}.log"
