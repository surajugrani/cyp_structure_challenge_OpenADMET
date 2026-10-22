#!/bin/bash
#SBATCH --job-name=a22_1924
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH --output=analyses/22_sub16_vs_17/slurm-%j.out
# RESOURCES: 1 serial process -> 1 CPU; RDKit + numpy over 20 structures, measured
#            2-4G for this kind of work -> 4G; the 19-vs-20 run took well under a
#            minute -> 20 min is already generous.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/22_sub16_vs_17/scripts/06_sub19_vs_24.py
echo "EXIT=$?"
