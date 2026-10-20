#!/bin/bash
#SBATCH --job-name=a22_1920
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH --output=analyses/22_sub16_vs_17/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/22_sub16_vs_17/scripts/05_sub19_vs_20.py
echo "EXIT=$?"
