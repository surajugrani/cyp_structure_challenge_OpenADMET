#!/bin/bash
#SBATCH --job-name=pbchk
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:15:00
#SBATCH --output=analyses/22_sub16_vs_17/slurm-pbchk-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/22_sub16_vs_17/scripts/02_energy_ratio_probe.py
