#!/bin/bash
#SBATCH --job-name=pymol_pse
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --output=tools/slurm-pymol-%j.out
set -uo pipefail
module load pymol/3.1.0
cd "$SLURM_SUBMIT_DIR"
pymol -cq tools/make_pymol_session.py
echo "EXIT=$?"
