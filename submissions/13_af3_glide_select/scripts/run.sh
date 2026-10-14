#!/bin/bash
#SBATCH --job-name=s13_gsel
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/13_af3_glide_select/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/13_af3_glide_select/scripts/01_select_and_package.py
echo "EXIT=$?"
