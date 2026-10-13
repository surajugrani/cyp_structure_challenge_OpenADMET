#!/bin/bash
#SBATCH --job-name=a14_ic50
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:30:00
#SBATCH --output=analyses/14_glide_vs_ic50/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/14_glide_vs_ic50/scripts/03_fit_model.py
echo "EXIT=$?"
