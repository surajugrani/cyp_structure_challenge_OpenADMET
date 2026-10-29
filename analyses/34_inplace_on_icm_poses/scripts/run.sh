#!/bin/bash
#SBATCH --job-name=a34_inplace
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 02:00:00
#SBATCH --output=analyses/34_inplace_on_icm_poses/slurm-%j.out
# RESOURCES: serial, 1 ligand at a time -> 1 CPU; one Glide process, measured under 2G
#            in analysis 28 -> 3G, the smallest ask that backfills soonest. No PrepWizard
#            and no grid generation here: both are reused from analysis 28, so the only
#            work is one inplace pass over 10 poses per ligand, ~1-2 min each; 20 ligands
#            -> well under an hour. 2 h covers a slow node. Resumable - it skips any
#            ligand already in the CSV.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/34_inplace_on_icm_poses/scripts/01_score.py "$@"
echo "EXIT=$?"
