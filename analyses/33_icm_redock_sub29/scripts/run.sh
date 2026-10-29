#!/bin/bash
#SBATCH --job-name=a33_icm
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 06:00:00
#SBATCH --output=analyses/33_icm_redock_sub29/slurm-%j.out
# RESOURCES: ONE ICM complex at a time in a single serial job - the ICM licence does not
#            allow concurrent docking, so 1 CPU is all that can be used. Analysis 17 ran
#            the same protocol at ~12-15 min per complex including project setup and
#            10-conformer docking; 20 complexes -> ~4-5 h, so 6 h with margin. Memory
#            measured under 3G there -> 4G. The job is resumable and skips anything
#            already in the CSV, so a timeout only costs a resubmission.
set -uo pipefail
module load schrodinger/2025-3
module load icm/3.9.4
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/33_icm_redock_sub29/scripts/01_redock.py "$@"
echo "EXIT=$?"
