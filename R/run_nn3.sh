#!/bin/bash
# ---------------------------------------------------------------------------
# run_nn3.sh
#
# Figure panels: S3 (all)
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_nn3.sh
#
# 3-NEAREST-NEURONAL-NEIGHBOUR SPACING analyses (Zwang et al. 2024 metric).
# Runs the three scripts IN ORDER -- the dependency chain is real, not stylistic:
#
#   1. nn3_phf1_vs_neg.R      Analysis A (cell-autonomous, PHF1+ vs PHF1-).
#                             ALSO writes the shared per-neuron cache
#                             plots/nn3_neuron_spacing/source_data_nn3_cells.tsv,
#                             which the next two scripts read instead of the
#                             327 MB SCE. Must run first.
#   2. nn3_over_phf1_distance.R  Analysis B (spacing vs distance to nearest PHF1+
#                             neuron, per subtype). Reads the cache; also writes
#                             the observed coefficients the null replica needs.
#   3. null_replica_nn3.R     SUPPLEMENTARY permutation null bounding the
#                             geometric confound. Reads the cache, the observed
#                             coefficients, and the 1.14 GB null distance matrix
#                             in deg/null_label_deg/engine/. This is the step that
#                             needs the memory allocation below.
#
# Steps 1 and 2 are cheap (~1 min each, dominated by I/O) and run fine on a login
# node or locally. Step 3 holds a 167176 x 1001 double matrix in memory, hence
# 120 gb requested.
#
# Prerequisite: R/generate_phf1_null_labels.r must already have been run
# (deg/null_label_deg/engine/ populated; seed 42, 1000 labellings).

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=4:mem=120gb
#PBS -N nn3_neuron_spacing
#PBS -o <PROJECT_ROOT>/phf1_v2/R/logs/nn3_neuron_spacing.o
#PBS -e <PROJECT_ROOT>/phf1_v2/R/logs/nn3_neuron_spacing.e

set -x
set -e   # the chain is order-dependent: if a step fails, do not run the next one

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

resource_dir=<PROJECT_ROOT>/phf1_v2/R
data_dir=<PROJECT_ROOT>/phf1_v2

mkdir -p "${resource_dir}/logs"
cd "${data_dir}"

# 1. Analysis A + shared cache
Rscript "${resource_dir}/nn3_phf1_vs_neg.R"

# 2. Analysis B
Rscript "${resource_dir}/nn3_over_phf1_distance.R"

# 3. SUPPLEMENTARY null. Drop --n_perm to subsample if memory is tight; note that
#    p_emp is bounded below by 2 / (1 + n_perm), so fewer labellings means a floor
#    on how small any empirical p-value can be.
Rscript "${resource_dir}/null_replica_nn3.R" --n_perm 1000

echo "Done. Outputs in ${data_dir}/plots/nn3_neuron_spacing"
