#!/bin/bash
#SBATCH --job-name=a35_reccheck
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --output=analyses/35_icm_receptor_check/slurm-%j.out
# RESOURCES: 20 ICM object->PDB conversions, each a few seconds, then pure-Python
#            distance arithmetic over 200 poses. 1 CPU is all it can use and the ICM
#            licence allows one process anyway. 4G/1h is well over what it needs and
#            backfills immediately.
set -uo pipefail
module load icm/3.9.4
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
python3 analyses/35_icm_receptor_check/scripts/02_check.py
echo "EXIT=$?"
