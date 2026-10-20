#!/bin/bash
#SBATCH --job-name=a24_ser
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:15:00
#SBATCH --output=analyses/24_preprint_residues/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/24_preprint_residues/scripts/02_ser119_contacts.py
