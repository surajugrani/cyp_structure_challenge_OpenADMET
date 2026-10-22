#!/bin/bash
#SBATCH --job-name=a25_1cpu
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 03:00:00
#SBATCH --output=analyses/25_glide_on_icm_poses/slurm-%j.out
# Derived: 5 concurrent Glide runs = 5 CPUs. Analysis 23 peaked at 7.5G doing MORE per
# task (prepwizard+grid+score, 4 concurrent) = ~1.9G each; scoring alone is lighter, so
# 5 x ~1.5G = 8G. 20 ligands x 10 poses, grids already exist, ~20 min expected.
# RESOURCES: serial variant - 1 concurrent Glide run -> 1 CPU; ~1.5G measured per run
#            -> 3G; 200 mininplace scorings at ~20 s each -> ~70 min, 3h margin.
#            Submitted to race the 5-CPU job 48196428, which is queued on (Priority).
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/25_glide_on_icm_poses/scripts/01_rescore.py --workers 1 "$@"
echo "EXIT=$?"
