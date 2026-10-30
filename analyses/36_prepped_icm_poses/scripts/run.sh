#!/bin/bash
#SBATCH --job-name=a36_prep
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 24:00:00
#SBATCH --output=analyses/36_prepped_icm_poses/slurm-%j.out
# RESOURCES: one PrepWizard, gridgen or Glide process at a time inside this job -> 1 CPU;
#            each measured under 2G in analysis 28 -> 3G, the smallest ask that backfills
#            soonest. Three of these run concurrently (shards 0/1/2), which is the agreed
#            Schrodinger limit and leaves headroom in the shared licence pool - analysis
#            32 lost two ligands to GLIDE_MAIN exhaustion when more were in flight.
#            Per pose: ~4.5 min prep + ~1 min grid + ~0.5 min score = ~6 min; ~67 poses
#            per shard -> ~7 h. 24 h covers licence waits. Resumable: anything already
#            recorded ok in the shard CSV is skipped.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/36_prepped_icm_poses/scripts/01_score.py "$@"
echo "EXIT=$?"
