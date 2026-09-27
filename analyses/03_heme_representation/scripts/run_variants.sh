#!/bin/bash
#SBATCH --job-name=cyp_heme_var
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
#
# Analysis 03: rebuild the receptor under several heme representations, regenerate the
# grid for each, rescore the crystal pose in place, and record every Glide term.
#
#     sbatch scripts/run_variants.sh --pdbs 8SPD,4K9X,2V0M
#
# Gate test first: if no variant gives a non-zero metal term, the idea is dead and there
# is no point running all 22.
set -euo pipefail
module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/01_variants_and_score.py "$@"
