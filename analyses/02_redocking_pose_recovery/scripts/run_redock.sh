#!/bin/bash
#SBATCH --job-name=cyp_redock
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#
# Analysis 02: redock each crystal ligand into its own receptor and measure how close
# the docked poses are to the crystal pose.
#
# Reuses the prepped receptors and grids from analysis 01, so no preparation is redone.
#
#     sbatch scripts/run_redock.sh              # all 22, 10 poses each
#     sbatch scripts/run_redock.sh --only 8SPD  # single structure

set -euo pipefail

module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1

cd "$SLURM_SUBMIT_DIR"

echo "=== Step 1: redock (Glide SP confgen) ==="
$SCHRODINGER/run python3 scripts/01_redock.py "$@"

echo "=== Step 2: symmetry-corrected RMSD to crystal pose ==="
$SCHRODINGER/run python3 scripts/02_rmsd.py

echo "=== Done ==="
