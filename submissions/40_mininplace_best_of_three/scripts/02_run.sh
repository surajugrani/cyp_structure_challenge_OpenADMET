#!/bin/bash
#SBATCH --job-name=s40_graft
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/40_mininplace_best_of_three/slurm-val-%j.out
# RESOURCES: RDKit and PoseBusters over 20 structures, single-threaded. Prior
#            validation passes on this set took ~3 min at under 2G; 4G and 30 min
#            is the smallest ask that comfortably covers it and backfills at once.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 \
    submissions/40_mininplace_best_of_three/scripts/02_graft.py
echo "EXIT=$?"
