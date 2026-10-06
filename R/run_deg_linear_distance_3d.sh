#!/bin/bash
# ---------------------------------------------------------------------------
# run_deg_linear_distance_3d.sh
#
# Figure panels: S1C
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_deg_linear_distance_3d.sh
#
# STAGE 2 of the out-of-plane tangle correction: the DEG sensitivity fits.
# One PBS array job per (distance-table configuration x replicate draw).
#
#   qsub R/run_deg_linear_distance_3d.sh            # all 27 jobs
#   qsub -J 1-1 R/run_deg_linear_distance_3d.sh     # the IDENTITY CHECK only
#
# ============================================================================
# RUN INDEX 1 FIRST AND CHECK IT BEFORE QUEUEING THE REST
# ============================================================================
# Index 1 is the identity arm: no simulated tangles at all, so the distance
# column is the canonical one, taken from the SCE's own stored values.
# Its output MUST match
#   deg/de_linear_distance/<ct>/<ct>_dist_to_phf1_um_scaled.tsv
# to floating-point noise (logFC to ~1e-12; identical gene ranking and
# significance calls). That is the proof the fork changed nothing but the
# distance. The check is scripted -- see the tail of this file.
#
# ============================================================================
# WHAT THESE FITS TEST
# ============================================================================
# Tangles are simulated ONLY in the planes the section never cut (|z| > 9 um,
# out to +/-500 um) at the volumetric density implied by the observed in-plane
# density, and each cell's distance becomes the 3D minimum over the observed
# in-plane anchors and those. NOTHING in the imaged plane is assumed missed and
# nothing is relabelled.
#
# The correction is concentrated in the TAIL -- median shortening runs 7.8% in
# the nearest d2D decile to 80.7% in the farthest -- and it DESTROYS DISTANCE
# ORDERING there (Spearman rho 0.79 clustered, 0.47 uniform). Once cells at
# 400-800 um all collapse to 100-200 um, which is nearest is decided by a
# simulated tangle. So:
#
#     survives  -> the effect is carried by the NEAR cells, where d3D ~= d2D.
#     collapses -> it rested on ordering in the far tail.
#
# --dilution_tsv is still passed because the DEG script requires it, but the
# 0.906 dilution constant was calibrated for the hidden-anchor arm
# (R/phf1_hidden_anchors.R) and does NOT apply here. Ignore
# attenuation_vs_dilution for these configs.
#
# ============================================================================
# ARRAY LAYOUT -- derived from dist3d/MANIFEST.tsv, NOT hardcoded here
# ============================================================================
#   idx 1        aug_p1p000_identity    rep_01           (identity check)
#   idx 2-11     unseen_clustered       rep_01..rep_10   (HEADLINE)
#   idx 12-21    unseen_uniform         rep_01..rep_10   (maximal-scrambling bound)
#   idx 22-24    unseen_clustered_h16   rep_01..rep_03   (thinner sampled slice)
#   idx 25-27    unseen_clustered_h20   rep_01..rep_03   (thicker sampled slice)
#                                                        = 27 jobs
#
# The arms simulate tangles in the planes the section never cut and take the 3D
# minimum against them. They do NOT assume any tangle in the imaged plane was
# missed. Note this test DESTROYS DISTANCE ORDERING IN THE TAIL (Spearman rho
# 0.79 clustered, 0.47 uniform), so survival means the effect is carried by the
# near cells.
# The job list is REBUILT from the manifest at submit time and the array range
# is asserted against it, so adding a configuration in Stage 1 cannot silently
# shift the mapping.

# MEASURED from pbs_out/ over 27 completed jobs: walltime 7:54-9:17, CPU
# 384-396% (4 cores well used). 30 min is >3x the worst case, and a short
# walltime is what gets a job scheduled when the queues are busy.
#PBS -l walltime=0:30:00
# MEASURED from pbs_out/ over 27 completed jobs: peak memory min 25.5 GB,
# median 26.5 GB, MAX 34.4 GB. 48gb is ~1.4x the worst case observed.
#
# Memory scales with the NUMBER OF GENES TESTED, frozen at 2468 for every arm,
# so this carries across all 27 jobs. Re-measure if the gene filter changes.
#PBS -l select=1:ncpus=4:mem=48gb
#PBS -N DEG_linear_distance_3d
#PBS -J 1-27
# PBS streams are KEPT, not sent to /dev/null. The Rscript output already goes
# to a per-job log, so these carry only the set -x trace and the PBS epilogue --
# and the epilogue is the only place resources_used.mem is recorded.
#PBS -o pbs_out/
#PBS -e pbs_out/

set -x

# set +u around the conda hook: with `set -u` an unbound ADDR2LINE kills
# `conda activate dgenv`. Stage 1 (which needs spatstat) is run separately,
# outside this environment.
set +u
eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"
conda activate dgenv
set -u

cd $PBS_O_WORKDIR
mkdir -p pbs_out

proj_dir=<PROJECT_ROOT>/phf1_v2
resource_dir="$proj_dir/R"
input_dir="$proj_dir/celltype_sce_neighbours"
output_dir="$proj_dir/deg"
dist_dir="$proj_dir/dist3d"
log_root="$proj_dir/R/logs"

CELLTYPE="Exc-IT-L2-3-CBLN2-HOPX"
MANIFEST="$dist_dir/MANIFEST.tsv"

if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: $MANIFEST not found. Run Stage 1 first, LOCALLY:"
  echo "  Rscript R/phf1_3d_distance.R --R 10 --R_corner 3"
  echo "(Stage 1 needs spatstat.)"
  exit 1
fi

# Expand the manifest into one line per (config, rep). Identity first, so that
# index 1 is always the correctness check regardless of manifest ordering.
JOBS=$(awk -F'\t' 'NR>1 {print $1"\t"$5"\t"$2}' "$MANIFEST" \
  | sort -k3,3 -t$'\t' \
  | awk -F'\t' '$3=="identity"{for(i=1;i<=$2;i++) printf "%s\trep_%02d\t%s\n",$1,i,$3}
                $3!="identity"{next}')
JOBS+=$'\n'
JOBS+=$(awk -F'\t' 'NR>1 && $2!="identity" {print $1"\t"$5"\t"$2}' "$MANIFEST" \
  | awk -F'\t' '{for(i=1;i<=$2;i++) printf "%s\trep_%02d\t%s\n",$1,i,$3}')
JOBS=$(echo "$JOBS" | sed '/^$/d')

N_JOBS=$(echo "$JOBS" | wc -l)
echo "Manifest expands to $N_JOBS jobs"

# Fail loudly rather than analysing the wrong table: the '#PBS -J' range above
# is static and the manifest is not.
if [ "${PBS_ARRAY_INDEX:-1}" -gt "$N_JOBS" ]; then
  echo "ERROR: PBS_ARRAY_INDEX=$PBS_ARRAY_INDEX exceeds the $N_JOBS jobs in the"
  echo "manifest. Update '#PBS -J 1-$N_JOBS' and resubmit."
  exit 1
fi

LINE=$(echo "$JOBS" | sed -n "${PBS_ARRAY_INDEX:-1}p")
CONFIG=$(echo "$LINE" | cut -f1)
REPCOL=$(echo "$LINE" | cut -f2)
TAG=$(echo "$LINE"    | cut -f3)

CONFIG_TAG="${CONFIG#aug_}_${REPCOL}"

log_dir="$log_root/de_linear_distance_3d_${CONFIG_TAG}"
mkdir -p "$log_dir"

echo "index=${PBS_ARRAY_INDEX:-1} config=$CONFIG col=$REPCOL tag=$TAG"

Rscript "$resource_dir/deg_dream_linear_distance_3d_phf1.r" \
  --sce              "$input_dir/${CELLTYPE}_sce_neighbours.qs" \
  --confounding_vars Sex,Age,PMI \
  --output_dir       "$output_dir" \
  --ncores           4 \
  --max_dist_um      1000 \
  --log_distance \
  --dist_table       "$dist_dir/${CONFIG}.tsv" \
  --dist_col         "$REPCOL" \
  --config_tag       "$CONFIG_TAG" \
  --dilution_tsv     "$proj_dir/results/tangle_3d/armB_dilution_constant.tsv" \
  > "$log_dir/${CELLTYPE}.log" 2>&1

status=$?

# The identity job self-checks against the canonical table. Anything other than
# agreement to floating-point noise means the fork is not a pure distance swap.
if [ "${PBS_ARRAY_INDEX:-1}" -eq 1 ] && [ $status -eq 0 ]; then
  Rscript -e '
    ct  <- "'"$CELLTYPE"'"
    new <- file.path("'"$output_dir"'",
             paste0("de_linear_distance_3d_", "'"$CONFIG_TAG"'"), ct,
             paste0(ct, "_dist_to_phf1_um_scaled.tsv"))
    ref <- file.path("'"$output_dir"'", "de_linear_distance", ct,
             paste0(ct, "_dist_to_phf1_um_scaled.tsv"))
    a <- read.delim(ref); b <- read.delim(new)
    m <- merge(a[, c("gene","logFC","padj")], b[, c("gene","logFC","padj")],
               by = "gene", suffixes = c("_ref","_new"))
    d <- max(abs(m$logFC_ref - m$logFC_new))
    s <- sum((m$padj_ref < 0.05) != (m$padj_new < 0.05))
    cat(sprintf("\n=== IDENTITY CHECK ===\ngenes matched: %d / %d\n", nrow(m), nrow(a)))
    cat(sprintf("max |dlogFC| = %.3g\nsignificance calls differing at padj<0.05: %d\n", d, s))
    if (d > 1e-8 || s > 0) {
      cat("FAILED. The fork is not a pure distance swap. Do NOT interpret the\n")
      cat("other array jobs until this is resolved.\n"); quit(status = 1)
    }
    cat("PASSED - the fork changes nothing but the distance.\n")
  ' 2>&1 | tee -a "$log_dir/${CELLTYPE}.log"
fi

exit $status
