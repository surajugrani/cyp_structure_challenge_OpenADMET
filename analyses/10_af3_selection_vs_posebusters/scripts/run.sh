#!/bin/bash
#SBATCH --job-name=af3_sel
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH --output=analyses/10_af3_selection_vs_posebusters/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module load schrodinger/2025-3
$SCHRODINGER/run python3 analyses/10_af3_selection_vs_posebusters/scripts/01_convert_all.py
~/micromamba/envs/cyp_posebusters_env/bin/python \
  analyses/10_af3_selection_vs_posebusters/scripts/02_pb_and_analyse.py
echo "EXIT=$?"
