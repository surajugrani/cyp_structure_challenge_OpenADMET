#!/bin/bash
#SBATCH --job-name=a21_lddt
#SBATCH -p shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --output=analyses/21_pocket_lddt_vs_crystals/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/pxr-eval/bin/python \
  analyses/21_pocket_lddt_vs_crystals/scripts/01_pocket_lddt.py
echo "EXIT=$?"
