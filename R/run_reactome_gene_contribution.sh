#!/bin/bash
# ---------------------------------------------------------------------------
# run_reactome_gene_contribution.sh
#
# Figure panels: S2B
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_reactome_gene_contribution.sh
#
# Decompose ONE Reactome module score's distance gradient into per-gene contributions,
# in ONE celltype. Runs plot_reactome_gene_contribution.R.
#
# Defaults: intrinsic apoptosis (R-HSA-109606) in Exc-IT-L2-3-CBLN2-HOPX.
# Change MODULE / CELLTYPE below to decompose a different one; MODULE must be a key from
# REACTOME_PRIMARY in plot_reactome_stress_death_vs_phf1_distance.R, which this script
# reads its configuration from so the two can never disagree.
#
# CAP_UM and SEED must match the run being explained, or the module score will differ
# from the one in the figure.
#
# Runs either way:
#   qsub run_reactome_gene_contribution.sh      # batch
#   bash R/run_reactome_gene_contribution.sh    # interactive session (preferred)

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=120gb
#PBS -N reactome_gene_contribution
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

project_root=<PROJECT_ROOT>/phf1_v2
cd "${PBS_O_WORKDIR:-$project_root}"

# ---- the inputs to change ----
# MODULE: 'significant' (every module the log-distance model called significant for this
#         celltype, read from the main run's coefficient table), 'all', or a
#         comma-separated list of REACTOME_PRIMARY keys e.g. "intrinsic_apop,hsf1".
CAP_UM=1000
SEED=42
MODULE=significant
CELLTYPE="Exc-IT-L2-3-CBLN2-HOPX"
LMM_COEFFS=plots/reactome_stress_death_vs_phf1_distance_${CAP_UM}um/stats_reactome_lmm_coeffs.tsv

resource_dir=$project_root/R
out_dir=plots/reactome_gene_contribution_${CAP_UM}um
log_dir=$project_root/R/logs/reactome_gene_contribution_${CAP_UM}um

mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START (module ${MODULE}, celltype ${CELLTYPE}, cap ${CAP_UM} um)"

Rscript $resource_dir/plot_reactome_gene_contribution.R \
  --seu          seu_PHF1.rds \
  --celltype     "$CELLTYPE" \
  --module       "$MODULE" \
  --lmm_coeffs   "$LMM_COEFFS" \
  --output_dir   "$out_dir" \
  --max_dist_um  ${CAP_UM} \
  --seed         ${SEED} \
  --overwrite    yes \
  > "$log_dir/reactome_gene_contribution_${MODULE}_${CAP_UM}um.log" 2>&1

END=$(date)
echo "Job ended at $END" >> "$log_dir/reactome_gene_contribution_${MODULE}_${CAP_UM}um.log"
