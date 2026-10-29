#!/bin/bash
#SBATCH --job-name=s34_val
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/34_glide_minimised/slurm-val-%j.out
# RESOURCES: RDKit and PoseBusters over 20 structures, single-threaded. Prior
#            validation passes on this set took ~3 min at under 2G; 4G and 30 min
#            is the smallest ask that comfortably covers it and backfills at once.
set -uo pipefail
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
    submissions/34_glide_minimised/scripts/02_validate.py
echo "EXIT=$?"
