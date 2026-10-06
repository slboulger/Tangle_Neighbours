#!/bin/bash
# ---------------------------------------------------------------------------
# run_deg_linear_distance.sh
#
# Figure panels: S1A; Table S8
# Runs the exclusion-radius sweep (10, 20, 30, 50 um) of deg_dream_linear_distance_phf1.r
# that collate_exclusion_radius.r collates.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_deg_linear_distance.sh
#
# PBS job array: one job per (celltype x exclusion radius).
# Runs deg_dream_linear_distance_phf1.r for each celltype SCE file at each radius
# of the exclusion-radius sensitivity sweep.
# Logs (stdout + stderr) written to logs/de_linear_distance_min<radius>um/<celltype>.log
# — a separate tree per radius, so no run overwrites another and the existing
# logs/de_linear_distance/*.log files from the canonical run are untouched.
#
# ARRAY RANGE DERIVATION
#   17 celltype SCE files in celltype_sce_neighbours/  x  4 radii (10 20 30 50)
#   = 68 jobs  ->  #PBS -J 1-68
#
#   Radius 0 is deliberately NOT in the array: it is the canonical run already on
#   disk in deg/de_linear_distance/, which this sweep must not overwrite. Run it
#   explicitly (with --min_distance_um 0) only if you intend to regenerate it.
#
#   Index arithmetic — celltype varies fastest, so each contiguous block of 17 is
#   one radius:
#     idx0   = PBS_ARRAY_INDEX - 1
#     ct_idx = idx0 % 17 + 1        (1-based index into `ls` of the SCE dir)
#     radius = RADII[idx0 / 17]     (integer division)
#
#   Submit everything:      qsub R/run_deg_linear_distance.sh
#   Submit one radius only: qsub -J 1-17  R/run_deg_linear_distance.sh   # 10 um
#                           qsub -J 18-34 R/run_deg_linear_distance.sh   # 20 um
#                           qsub -J 35-51 R/run_deg_linear_distance.sh   # 30 um
#                           qsub -J 52-68 R/run_deg_linear_distance.sh   # 50 um

#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=4:mem=100gb
#PBS -N DEG_linear_distance
#PBS -J 1-68
#PBS -o /dev/null
#PBS -e /dev/null

set -x

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv

cd $PBS_O_WORKDIR

resource_dir=<PROJECT_ROOT>/phf1_v2/R
input_dir=<PROJECT_ROOT>/phf1_v2/celltype_sce_neighbours
output_dir=<PROJECT_ROOT>/phf1_v2/deg
log_root=<PROJECT_ROOT>/phf1_v2/R/logs

# Exclusion radii for the sweep (um). Radius 0 = the canonical run, excluded.
RADII=(10 20 30 50)
N_CT=17

# The array range is baked into '#PBS -J' above and assumes exactly N_CT files in
# input_dir. If a celltype is ever added or removed the mapping silently shifts,
# so fail loudly rather than analysing the wrong file.
n_files=$(ls -1 "$input_dir" | wc -l)
if [ "$n_files" -ne "$N_CT" ]; then
  echo "ERROR: expected $N_CT SCE files in $input_dir, found $n_files."
  echo "Update N_CT and the '#PBS -J' range (N_CT * ${#RADII[@]}) and resubmit."
  exit 1
fi

idx0=$(( PBS_ARRAY_INDEX - 1 ))
ct_idx=$(( idx0 % N_CT + 1 ))
r_idx=$(( idx0 / N_CT ))
radius=${RADII[$r_idx]}

if [ -z "$radius" ]; then
  echo "ERROR: PBS_ARRAY_INDEX=$PBS_ARRAY_INDEX maps to radius index $r_idx,"
  echo "which is outside RADII (${RADII[*]}). Check the '#PBS -J' range."
  exit 1
fi

# Quote everything: celltype names can contain spaces (e.g. "Unassigned Neuron").
input_sce=$(ls "$input_dir" | head -n $ct_idx | tail -n 1)
celltype=$(basename "$input_sce" _sce_neighbours.qs)

log_dir="$log_root/de_linear_distance_min${radius}um"
mkdir -p "$log_dir"

START=$(date)
echo "Job started at $START"
echo "Array index: $PBS_ARRAY_INDEX  (celltype index $ct_idx of $N_CT, radius index $r_idx)"
echo "Input SCE: $input_sce"
echo "Celltype: $celltype"
echo "Exclusion radius (um): $radius"

Rscript $resource_dir/deg_dream_linear_distance_phf1.r \
  --sce              "$input_dir/$input_sce" \
  --confounding_vars Sex,Age,PMI \
  --output_dir       "$output_dir" \
  --ncores           4 \
  --max_dist_um      1000 \
  --min_distance_um  "$radius" \
  --log_distance \
  > "$log_dir/${celltype}.log" 2>&1

END=$(date)
echo "Job ended at $END" >> "$log_dir/${celltype}.log"
