#!/bin/bash
#SBATCH --job-name=s10_hliptm
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/10_af3_hliptm/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/10_af3_hliptm/scripts/01_select_and_package.py
echo "EXIT=$?"
