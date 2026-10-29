#!/bin/bash
#SBATCH --job-name=a32_inplace
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 03:00:00
#SBATCH --output=analyses/32_inplace_depth2048/slurm-%j.out
# RESOURCES: serial, 1 structure at a time -> 1 CPU; one PrepWizard or Glide process at
#            a time, measured under 2G in analysis 28 -> 3G, the smallest ask that
#            backfills soonest. Per structure ~4.5 min prep + ~1 min grid + ~1 min for
#            the single inplace pass; 20 structures serial -> ~2.2 h. The job is
#            resumable: it skips anything already in the CSV, so a short ask that has to
#            be repeated costs less than a long one waiting for a slot.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/32_inplace_depth2048/scripts/01_score.py --workers 1 "$@"
echo "EXIT=$?"
