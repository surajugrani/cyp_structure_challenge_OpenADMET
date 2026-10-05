#!/bin/bash
#SBATCH --job-name=a13_rmsd2
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH --output=analyses/13_pose_divergence/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/13_pose_divergence/scripts/02_pool_rmsd_notemplate.py
echo "EXIT=$?"
