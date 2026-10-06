#!/bin/bash
# ---------------------------------------------------------------------------
# run_generate_phf1_null_labels.sh
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
# run_generate_phf1_null_labels.sh
#
# ONE-TIME global engine for the PHF1 label-permutation nulls (docs/MODELS.md,
# Permutation nulls). Produces phf1_label_permutations.qs +
# phf1_null_distance_matrix.qs under deg/null_label_deg/engine/. Run this ONCE
# before any consumer job, and submit consumers with a dependency on this job:
#   jid=$(qsub run_generate_phf1_null_labels.sh)
#   qsub -W depend=afterok:$jid <consumer job script>

#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=8:mem=120gb
#PBS -N phf1_null_labels
#PBS -o <PROJECT_ROOT>/phf1_v2/R/logs/phf1_null_labels.o
#PBS -e <PROJECT_ROOT>/phf1_v2/R/logs/phf1_null_labels.e

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

resource_dir=<PROJECT_ROOT>/phf1_v2/R
data_dir=<PROJECT_ROOT>/phf1_v2
output_dir=${data_dir}/deg/null_label_deg/engine

mkdir -p "${data_dir}/R/logs"

START=$(date)
echo "Job started at $START"

# Distance-source neuron celltypes MUST match label_phf1_neighbours.r exactly.
Rscript $resource_dir/generate_phf1_null_labels.r \
  --all_sce          ${data_dir}/sce.qs \
  --neuron_celltypes "Exc-IT-L2-3-CBLN2-HOPX,Exc-IT-L3-5-CHGA-IL1RAPL2,Exc-ET-L5-SPON1-FGD4,Exc-CT-L6-SYNPO2-SEMA3E,Exc-IT-L6-CTXN1-ERC2,Inh-SST,Inh-PVALB,Inh-VIP,Inh-LAMP5,Unassigned Neuron" \
  --coord_x          x_slide_mm \
  --coord_y          y_slide_mm \
  --sample_col       sample_id \
  --n_perm           1000 \
  --seed             42 \
  --ncores           8 \
  --output_dir       "${output_dir}"

END=$(date)
echo "Job ended at $END"
