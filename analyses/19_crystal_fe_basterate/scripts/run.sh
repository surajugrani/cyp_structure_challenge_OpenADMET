#!/bin/bash
#SBATCH --job-name=a19_fe
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:15:00
#SBATCH --output=analyses/19_crystal_fe_basterate/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/19_crystal_fe_basterate/scripts/01_crystal_fe.py
echo "EXIT=$?"
