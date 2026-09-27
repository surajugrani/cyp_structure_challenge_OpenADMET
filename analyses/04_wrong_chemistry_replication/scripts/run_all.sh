#!/bin/bash
#SBATCH --job-name=cyp_wrongchem
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#
# Analysis 04: repeat analyses 01 + 02 under the metal-term-activating (chemically wrong)
# heme representation, fe3_monoanion.
#
#     sbatch scripts/run_all.sh --only 8SPD   # smoke test
#     sbatch scripts/run_all.sh               # all 22
set -euo pipefail
module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/01_build_score_redock.py "$@"
