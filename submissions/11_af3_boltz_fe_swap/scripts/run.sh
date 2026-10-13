#!/bin/bash
#SBATCH --job-name=s11_swap
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/11_af3_boltz_fe_swap/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/11_af3_boltz_fe_swap/scripts/01_build.py
echo "EXIT=$?"
