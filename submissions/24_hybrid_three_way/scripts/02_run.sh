#!/bin/bash
#SBATCH --job-name=s24_3way
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/24_hybrid_three_way/slurm-%j.out
# RESOURCES: 1 serial process -> 1 CPU; PoseBusters + RDKit over 20 structures measured
#            at 2-4G in submission 22 -> 4G; that job ran ~2.5 min -> 30 min.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/24_hybrid_three_way/scripts/01_build.py
echo "EXIT=$?"
# Always build the PyMOL session alongside the zip (no-op if no zip was written).
bash tools/make_pse.sh submissions/24_hybrid_three_way
