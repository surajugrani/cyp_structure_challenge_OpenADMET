#!/bin/bash
#SBATCH --job-name=a18_glide
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=analyses/18_glide_redock_predicted/slurm-%j.out
#
# 20 Glide redockings, ONE AT A TIME, single serial job (licence limit).
# Resumable: resubmit to continue.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/18_glide_redock_predicted/scripts/01_glide_redock.py "$@"
echo "EXIT=$?"
