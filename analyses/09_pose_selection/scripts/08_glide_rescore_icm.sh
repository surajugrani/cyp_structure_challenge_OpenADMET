#!/bin/bash
#SBATCH --job-name=gmin_icm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=analyses/09_pose_selection/slurm-gonicm-%j.out
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/09_pose_selection/scripts/08_glide_rescore_icm.py
