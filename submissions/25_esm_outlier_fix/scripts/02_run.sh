#!/bin/bash
#SBATCH --job-name=s25_fix
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/25_esm_outlier_fix/slurm-%j.out
# RESOURCES: 1 serial process -> 1 CPU; PoseBusters + RDKit measured at 2-4G in
#            submissions 22 and 24 -> 4G; two 20-structure passes, those ran ~2.5 min
#            each -> 30 min.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/25_esm_outlier_fix/scripts/01_build.py
echo "EXIT=$?"
# Always build the PyMOL session alongside the zip (no-op if no zip was written).
bash tools/make_pse.sh submissions/25_esm_outlier_fix
