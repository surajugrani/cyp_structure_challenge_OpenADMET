#!/bin/bash
#SBATCH --job-name=s14_rtcnn
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/14_af3_icm_rtcnn/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/14_af3_icm_rtcnn/scripts/01_select_and_package.py
echo "EXIT=$?"
