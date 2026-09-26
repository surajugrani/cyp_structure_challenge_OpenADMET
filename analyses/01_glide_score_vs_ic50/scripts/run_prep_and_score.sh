#!/bin/bash
#SBATCH --job-name=cyp_glide_ic50
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out
#
# Steps 3 + 4 of analysis 01: PrepWizard, Glide grids, and inplace/mininplace scoring
# of the crystal poses for the 22 curated CYP3A4 structures.
#
# Submit from the analysis folder:
#     sbatch scripts/run_prep_and_score.sh
# Optionally restrict to a few structures:
#     sbatch scripts/run_prep_and_score.sh --only 2V0M,4K9T

set -euo pipefail

module load schrodinger/2025-3

# Compute nodes lack libnsl.so.1; a copy lives in $HOME (same workaround as the PXR runs)
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}

cd "$SLURM_SUBMIT_DIR"
echo "SCHRODINGER = $SCHRODINGER"
echo "Working in  = $PWD"

echo "=== Step 3: prepare structures ==="
$SCHRODINGER/run python3 scripts/03_prep_structures.py "$@"

echo "=== Step 4: grids + inplace/mininplace scoring ==="
$SCHRODINGER/run python3 scripts/04_grid_and_score.py "$@"

echo "=== Done ==="
