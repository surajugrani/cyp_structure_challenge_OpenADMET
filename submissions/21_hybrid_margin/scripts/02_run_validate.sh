#!/bin/bash
#SBATCH --job-name=s21_val
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/21_hybrid_margin/slurm-%j.out
# Resources derived: single-threaded validation of 20 structures; PoseBusters peaks ~2G.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/21_hybrid_margin/scripts/01_validate_package.py
echo "EXIT=$?"
# Always build the PyMOL session alongside the zip (no-op if no zip was written).
bash tools/make_pse.sh submissions/21_hybrid_margin
