#!/bin/bash
#SBATCH --job-name=a20_icm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=analyses/20_icm_score_af3/slurm-%j.out
#
# 100 ICM in-place scorings, ONE AT A TIME (licence limit). Resumable.
set -uo pipefail
module load icm/3.9.4
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/20_icm_score_af3/scripts/01_score_samples.py "$@"
echo "EXIT=$?"
