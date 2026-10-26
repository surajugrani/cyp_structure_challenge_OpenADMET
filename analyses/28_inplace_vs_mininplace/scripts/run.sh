#!/bin/bash
#SBATCH --job-name=a28_modes
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 04:00:00
#SBATCH --output=analyses/28_inplace_vs_mininplace/slurm-%j.out
# RESOURCES: serial, 1 structure at a time -> 1 CPU; one PrepWizard or Glide process
#            measured <2G -> 3G, the smallest ask that backfills soonest; per structure
#            ~4.5 min prep + ~1 min grid + 2 scoring passes ~1 min each, 60 serial
#            -> ~7.5 h. Asking 4h rather than 10h deliberately: the job is resumable,
#            it skips anything already in the CSV, and a short ask backfills far sooner
#            than a long one waits for a slot. Resubmit to finish.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/28_inplace_vs_mininplace/scripts/01_score_both.py --workers 1 "$@"
echo "EXIT=$?"
