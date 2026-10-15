#!/bin/bash
# ---------------------------------------------------------------------------
# run_traj_pipeline.sh
#
# Figure panels: 2B-2D, 2G; Tables S4, S5
# Runs the tangle-associated trajectory workflow for one celltype, in order.
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
# run_traj_pipeline.sh
#
# The whole tangle-associated trajectory workflow for ONE celltype, in order, with the
# conda env per step. Every script takes the celltype, and files are named by ct_tag()
# (R/cbln2_embedding_utils.R): Exc-IT-L2-3-CBLN2-HOPX -> "cbln2",
# Exc-IT-L3-5-CHGA-IL1RAPL2 -> "chga".
#
#   1  embedding + PHF1-p95 UMAP      R/plot_cbln2_umap_phf1_p95.R         $TRAJ_ENV (SCT needs glmGamPoi)
#   2  trajectory search (tangle %)   R/cbln2_trajectory_optimise.R        $TRAJ_ENV (slingshot)
#   3  fit + describe trajectory      R/cbln2_slingshot_trajectory.R       $TRAJ_ENV
#   4  % PHF1+ by bin, mean + CI      R/plot_cbln2_traj_phf1_by_bin_mean.R $TRAJ_ENV
#   5  dream DEG along pseudotime     R/deg_dream_linear_pseudotime_phf1.r $DEG_ENV
#   8  enrichR ORA                    R/pathway_enrichment_cbln2_trajectory.R  $POST_ENV, needs INTERNET (Enrichr API)
#   9  p for tangle % rising          R/cbln2_trajectory_selection_null.R  $TRAJ_ENV (whole search re-run on shuffled
#                                     labels; then re-runs step 4 so the figure carries the p)
#   Steps 6-7 of the working pipeline were exploratory and are omitted from this release;
#   the remaining steps keep their original numbers.
#   All envs default to dgenv (TRAJ_ENV / DEG_ENV / POST_ENV to override).
#   -  pathway dotplot: interactive two-phase selection, run by hand:
#        R/plot_dotplot_pathway_linear_pseudotime_phf1.r with `celltype` set in its CONFIG.
#
# RUN in an interactive session (bash, conda available; no need to activate an env first):
#   bash R/run_traj_pipeline.sh Exc-IT-L2-3-CBLN2-HOPX 8
#   bash R/run_traj_pipeline.sh Exc-IT-L2-3-CBLN2-HOPX 8 3      # resume from step 3
# Args: celltype, cores (default 8 -- set to your allocation), first step (default 1).
# Each step logs to R/logs/traj_pipeline/<celltype>/<step>.log and the script stops at the
# first failure (set -e). If step 8 fails for lack of internet, run it where there is:
#   Rscript R/pathway_enrichment_cbln2_trajectory.R <celltype>
#
# Settings: donors with >= 5 PHF1+ cells enter the trajectory objective and the AUC
# summaries (MIN_POS_DONOR below; the CBLN2 donor set is the same at 5 or 10).

set -eo pipefail   # no -u: conda's activate scripts reference unset variables

CELLTYPE="${1:-Exc-IT-L3-5-CHGA-IL1RAPL2}"
NCORES="${2:-8}"
FROM="${3:-1}"
MIN_POS_DONOR=5
# Conda env per stage. The whole workflow runs in dgenv (Seurat/SCT+glmGamPoi, harmony,
# slingshot, variancePartition, enrichR, ggrepel, limma). Override if needed, e.g.
#   TRAJ_ENV=Seuratenv bash R/run_traj_pipeline.sh <celltype> 8
TRAJ_ENV="${TRAJ_ENV:-dgenv}"   # steps 1-4 (embedding, trajectory)
DEG_ENV="${DEG_ENV:-dgenv}"     # step 5 (dream)
POST_ENV="${POST_ENV:-dgenv}"   # step 8 (enrichment)

base=<PROJECT_ROOT>/phf1_v2
cd "$base"
log_dir="$base/R/logs/traj_pipeline/$CELLTYPE"
mkdir -p "$log_dir"

eval "$(<CONDA_ROOT>/bin/conda shell.bash hook)"

if [ ! -f "celltype_sce/${CELLTYPE}_sce.qs" ]; then
  echo "ERROR: celltype_sce/${CELLTYPE}_sce.qs not found."; exit 1
fi

step () {   # step <n> <env> <label> <command...>
  local n="$1" env="$2" label="$3"; shift 3
  if [ "$n" -lt "$FROM" ]; then echo "[skip] step $n $label"; return; fi
  conda activate "$env"
  echo "[$(date +%H:%M:%S)] step $n ($env): $label  -> $log_dir/step${n}.log"
  "$@" > "$log_dir/step${n}.log" 2>&1 || { echo "FAILED at step $n -- see $log_dir/step${n}.log"; exit 1; }
  conda deactivate
}

echo "Trajectory pipeline: $CELLTYPE  (cores $NCORES, from step $FROM, min_pos_donor $MIN_POS_DONOR)"
step 1 "$TRAJ_ENV" "embedding + PHF1-p95 UMAP" \
  Rscript R/plot_cbln2_umap_phf1_p95.R --celltype "$CELLTYPE"
step 2 "$TRAJ_ENV" "trajectory search (follows tangle %)" \
  Rscript R/cbln2_trajectory_optimise.R --celltype "$CELLTYPE" --n_cores "$NCORES" --min_pos_donor "$MIN_POS_DONOR"
step 3 "$TRAJ_ENV" "fit + describe trajectory" \
  Rscript R/cbln2_slingshot_trajectory.R --celltype "$CELLTYPE"
step 4 "$TRAJ_ENV" "% PHF1+ by pseudotime bin, mean + CI" \
  Rscript R/plot_cbln2_traj_phf1_by_bin_mean.R "$CELLTYPE"
step 5 "$DEG_ENV" "dream DEG along pseudotime" \
  Rscript R/deg_dream_linear_pseudotime_phf1.r --celltype "$CELLTYPE" --confounding_vars Sex,Age,PMI \
          --output_dir deg --ncores "$NCORES"
step 8 "$POST_ENV" "enrichR ORA (needs internet)" \
  Rscript R/pathway_enrichment_cbln2_trajectory.R "$CELLTYPE"

step 9 "$TRAJ_ENV" "selection-adjusted p for tangle % along pseudotime (+ redraw step 4)" \
  bash -c "Rscript R/cbln2_trajectory_selection_null.R --celltype '$CELLTYPE' --n_cores $NCORES && Rscript R/plot_cbln2_traj_phf1_by_bin_mean.R '$CELLTYPE'"

echo "[$(date +%H:%M:%S)] done. Pathway dotplot is interactive: set celltype in R/plot_dotplot_pathway_linear_pseudotime_phf1.r."
