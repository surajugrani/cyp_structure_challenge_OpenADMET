#!/bin/bash
#SBATCH --job-name=a25_resc
#SBATCH --partition=shared
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --output=analyses/25_glide_on_icm_poses/slurm-%j.out
# Derived: 5 concurrent Glide runs = 5 CPUs. Analysis 23 peaked at 7.5G doing MORE per
# task (prepwizard+grid+score, 4 concurrent) = ~1.9G each; scoring alone is lighter, so
# 5 x ~1.5G = 8G. 20 ligands x 10 poses, grids already exist, ~20 min expected.
# RESOURCES: see the derivation comment above; verified against sacct MaxRSS.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/25_glide_on_icm_poses/scripts/01_rescore.py "$@"
echo "EXIT=$?"
