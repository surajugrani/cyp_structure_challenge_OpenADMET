#!/bin/bash
#SBATCH --job-name=s22_swap
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/22_icm_fe_swap/slurm-%j.out
# RESOURCES: 1 serial process -> 1 CPU; PoseBusters + RDKit over 20 structures measured
#            at 2-4G -> 4G; comparable validation jobs ran ~2.5 min -> 30 min.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/22_icm_fe_swap/scripts/01_build.py
echo "EXIT=$?"
# Always build the PyMOL session alongside the zip (no-op if no zip was written).
bash tools/make_pse.sh submissions/22_icm_fe_swap
