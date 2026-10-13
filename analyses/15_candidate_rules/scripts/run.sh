#!/bin/bash
#SBATCH --job-name=a15_prev
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --output=analyses/15_candidate_rules/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/15_candidate_rules/scripts/01_preview.py
echo "EXIT=$?"
