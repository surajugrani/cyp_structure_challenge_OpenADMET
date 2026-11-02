#!/bin/bash
#SBATCH --job-name=a36_prep
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 4:00:00
#SBATCH --output=analyses/36_prepped_icm_poses/slurm-%j.out
# RESOURCES: one Glide process at a time inside this job -> 1 CPU; measured under 2G in
#            analysis 28 -> 3G, the smallest ask that backfills soonest. Three of these
#            run concurrently (shards 0/1/2), the agreed Schrodinger limit, which leaves
#            headroom in the shared licence pool - analysis 32 lost two ligands to
#            GLIDE_MAIN exhaustion when more were in flight.
#            The first run completed PrepWizard and gridgen for all 200 poses and they
#            are cached on disk, so only Glide inplace is left: ~0.5 min per pose, ~67
#            poses per shard -> ~35 min. 4 h covers licence waits. Resumable: anything
#            already recorded ok in the shard CSV is skipped.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/36_prepped_icm_poses/scripts/01_score.py "$@"
echo "EXIT=$?"
