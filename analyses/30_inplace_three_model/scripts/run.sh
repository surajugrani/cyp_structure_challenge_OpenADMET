#!/bin/bash
#SBATCH --job-name=a30_inplace
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 06:00:00
#SBATCH --output=analyses/30_inplace_three_model/slurm-%j.out
# RESOURCES: serial, 1 structure at a time -> 1 CPU; measured 2.1G peak in analysis 28
#            -> 3G; 39 structures at ~6 min each (prep + grid + ONE scoring pass)
#            -> ~4 h. Resumable: it skips anything already in the CSV, so a timeout
#            costs only a resubmit, which is why the ask is tight rather than padded.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/30_inplace_three_model/scripts/01_score.py --workers 1 "$@"
echo "EXIT=$?"
