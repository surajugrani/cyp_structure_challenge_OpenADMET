#!/bin/bash
#SBATCH --job-name=glide_xtal
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=analyses/09_pose_selection/slurm-gxtal-%j.out
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/09_pose_selection/scripts/05_glide_crystal.py
