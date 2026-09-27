#!/bin/bash
#SBATCH --job-name=cyp_icm_follow
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=slurm-follow-%j.out
#
# Analysis 06: dock with ICM, following behind the Glide run.
#
# ICM reuses analysis 05's LigPrep ligands (option b - ligand prep held constant), so it
# can only process structures Glide has already prepared. Rather than requiring repeated
# manual resubmission, this job loops: prepare inputs for whatever is ready, dock them,
# then sleep and check again. It exits when every structure in the set is done.
#
# ICM uses a separate licence from Schrodinger, so this runs alongside the Glide job.
#
#     sbatch scripts/run_icm_follow.sh
set -euo pipefail
module load icm/3.9.4
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"

TOTAL=$(tail -n +2 ../05_full_pdb_redocking/data/structure_set.csv | wc -l)
IDLE=0
while true; do
  $SCHRODINGER/run python3 scripts/01_make_icm_inputs.py  >/dev/null 2>&1 || true
  $SCHRODINGER/run python3 scripts/02_icm_dock.py

  DONE=$( [ -f results/icm_summary.csv ] && tail -n +2 results/icm_summary.csv | wc -l || echo 0 )
  echo "--- progress: $DONE/$TOTAL docked with ICM"
  [ "$DONE" -ge "$TOTAL" ] && { echo "all structures done"; break; }

  # Stop if Glide has finished and we still cannot progress (nothing new to do)
  if ! squeue -u "$USER" -h -o %j | grep -q cyp_fullredock; then
    IDLE=$((IDLE+1))
    [ "$IDLE" -ge 2 ] && { echo "Glide finished and no new inputs appeared; stopping"; break; }
  fi
  sleep 600
done

echo "=== final ==="
$SCHRODINGER/run python3 - <<'PY'
import pandas as pd
d = pd.read_csv("results/icm_summary.csv")
ok = d[d.status == "ok"]
print(f"{len(ok)}/{len(d)} succeeded")
for label, col in (("top-1", "rmsd_top1"), ("best-of-N", "rmsd_best")):
    for cut in (2.0, 1.0):
        n = (ok[col] <= cut).sum()
        print(f"  {label:<10} RMSD <= {cut:.1f} A: {n:3d}/{len(ok)} ({100*n/len(ok):.0f} %)")
    print(f"  {label:<10} median: {ok[col].median():.2f} A")
PY
