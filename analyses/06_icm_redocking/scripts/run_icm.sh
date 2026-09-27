#!/bin/bash
#SBATCH --job-name=cyp_icm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=slurm-%j.out
#
# Analysis 06: dock all CYP3A4 structures with ICM-Pro.
#
# Uses a separate licence from Schrodinger, so this can run alongside the Glide job.
# Resumable and incremental - resubmit to continue where it stopped.
#
#     sbatch scripts/run_icm.sh --limit 2   # smoke test
#     sbatch scripts/run_icm.sh             # everything with inputs prepared
set -euo pipefail
module load icm/3.9.4
module load schrodinger/2025-3   # for the structure reader and spyrmsd path only
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/02_icm_dock.py "$@"
