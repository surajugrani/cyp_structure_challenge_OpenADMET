#!/bin/bash
#SBATCH --job-name=cyp_correlate
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
#
# Step 5: correlate Glide scores with pIC50.
#
#     sbatch scripts/run_correlate.sh                  # docking_score + evdw
#     sbatch scripts/run_correlate.sh --terms metal,ecoul
#
# Light on CPU, but submitted like everything else rather than run on the login node.

set -euo pipefail

module load schrodinger/2025-3
export PYTHONUNBUFFERED=1

cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/05_correlate.py "$@"
