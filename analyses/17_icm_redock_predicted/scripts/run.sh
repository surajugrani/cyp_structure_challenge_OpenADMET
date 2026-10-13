#!/bin/bash
#SBATCH --job-name=a17_icm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=analyses/17_icm_redock_predicted/slurm-%j.out
#
# 40 dockings (20 ligands x 2 receptors), ONE AT A TIME, single serial job.
# Resumable: resubmit to continue where it stopped.
set -uo pipefail
module load icm/3.9.4
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/17_icm_redock_predicted/scripts/01_icm_redock.py "$@"
echo "EXIT=$?"
