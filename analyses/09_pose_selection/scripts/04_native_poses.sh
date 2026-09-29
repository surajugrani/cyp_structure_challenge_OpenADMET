#!/bin/bash
#SBATCH --job-name=icm_native
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=analyses/09_pose_selection/slurm-native-%j.out
set -uo pipefail
module load icm/3.9.4
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/09_pose_selection/scripts/04_native_poses.py
