#!/bin/bash
#SBATCH --job-name=metal_cons
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 08:00:00
#SBATCH --output=analyses/12_metal_constraint/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module load schrodinger/2025-3
$SCHRODINGER/run python3 analyses/12_metal_constraint/scripts/01_dock_constrained.py
echo "EXIT=$?"
