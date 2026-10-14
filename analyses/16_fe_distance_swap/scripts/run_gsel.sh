#!/bin/bash
#SBATCH --job-name=a16_gsel
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:15:00
#SBATCH --output=analyses/16_fe_distance_swap/slurm-gsel-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/16_fe_distance_swap/scripts/03_glide_selection_fe.py
echo "EXIT=$?"
