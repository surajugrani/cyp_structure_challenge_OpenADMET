#!/bin/bash
#SBATCH --job-name=a16_fe
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --output=analyses/16_fe_distance_swap/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/16_fe_distance_swap/scripts/01_fe_dist.py
echo "EXIT=$?"
