#!/bin/bash
#SBATCH --job-name=cyp_rmsd
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
#
# Recompute pose RMSDs from existing redocking output, without re-docking.
#     sbatch scripts/run_rmsd.sh
set -euo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/02_rmsd.py
