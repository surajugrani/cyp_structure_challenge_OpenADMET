#!/bin/bash
#SBATCH --job-name=cyp_glide_test
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-test-%j.out
#
# Smoke test: run the full step 3 + 4 pipeline on a single structure before
# committing the whole set to the queue. Never run these on the login node -
# PrepWizard minimization and Glide grid generation are CPU-heavy.
#
#     sbatch scripts/test_one.sh 8SPD
#
# Check progress with:  squeue -u $USER    and    tail -f slurm-test-<jobid>.out

set -euo pipefail

PDB=${1:-8SPD}

module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}

cd "$SLURM_SUBMIT_DIR"
echo "Testing pipeline on $PDB"

$SCHRODINGER/run python3 scripts/03_prep_structures.py --only "$PDB"
$SCHRODINGER/run python3 scripts/04_grid_and_score.py --only "$PDB"

echo "=== Scores ==="
cat results/glide_scores.csv
