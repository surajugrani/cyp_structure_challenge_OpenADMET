#!/bin/bash
#SBATCH --job-name=a29_icm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=analyses/29_icm_on_esm/slurm-%j.out
# RESOURCES: ONE ICM complex at a time (licence limit) -> 1 CPU, never more. Measured
#            baseline says 16G is ample for ICM scoring; this is scoring not docking,
#            so 4G, and a smaller ask backfills sooner. Analysis 20 did 100 structures
#            serially in ~1 h, so 60 is ~40 min; 2h is already 3x that, and the job is
#            resumable so a timeout only costs a resubmit.
set -uo pipefail
module load icm/3.9.4
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/29_icm_on_esm/scripts/01_icm_score.py "$@"
echo "EXIT=$?"
