#!/bin/bash
#SBATCH --job-name=a22_er
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:40:00
#SBATCH --output=analyses/22_sub16_vs_17/slurm-a22_er-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/22_sub16_vs_17/scripts/03_energy_ratio_16_17.py
