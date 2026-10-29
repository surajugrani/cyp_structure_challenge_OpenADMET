#!/bin/bash
#SBATCH --job-name=a32_retry
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 06:00:00
#SBATCH --output=analyses/32_inplace_depth2048/slurm-retry-%j.out
# RESOURCES: the last 2 of 20 ligands cannot start because the shared Schrodinger
#            GLIDE_MAIN pool is exhausted by other users - the grids are built and the
#            work is ~2 min once a token frees. This waits rather than failing: up to 60
#            attempts 5 min apart, exiting the moment both ligands score. 1 CPU and 3G
#            because only one Glide process ever runs; 6 h bounds the wait. Idle time is
#            sleep, not compute.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
CSV=analyses/32_inplace_depth2048/results/inplace_depth2048.csv
for i in $(seq 1 60); do
  n=$(grep -cE ',ok,' "$CSV" || true)
  if [ "$n" -ge 20 ]; then echo "all 20 scored after $((i-1)) attempts"; break; fi
  # drop rows that recorded a licence failure so the resumable script retries them
  grep -vE ',error:' "$CSV" > "$CSV.t" && mv "$CSV.t" "$CSV"
  rm -rf analyses/32_inplace_depth2048/work/31_depth2048/*/inplace
  echo "--- attempt $i ($n/20 scored)"
  $SCHRODINGER/run python3 analyses/32_inplace_depth2048/scripts/01_score.py --workers 1
  n=$(grep -cE ',ok,' "$CSV" || true)
  if [ "$n" -ge 20 ]; then echo "all 20 scored on attempt $i"; break; fi
  sleep 300
done
grep -cE ',ok,' "$CSV"
echo "EXIT=$?"
