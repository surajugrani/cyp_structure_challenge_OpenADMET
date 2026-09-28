#!/bin/bash
#SBATCH --job-name=icm_scoreonly
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#
# Analysis 08: score the crystal poses with ICM (Score + RTCNN) via `_dockScan -n`,
# the ICM counterpart to analysis 01's Glide inplace scoring.
#
#     sbatch scripts/run_score_only.sh
set -uo pipefail
module load icm/3.9.4
module load schrodinger/2025-3      # structure reader only
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/01_score_crystal_poses.py
