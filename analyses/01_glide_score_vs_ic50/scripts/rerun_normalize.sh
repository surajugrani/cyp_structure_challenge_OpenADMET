#!/bin/bash
#SBATCH --job-name=cyp_renorm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
#
# Re-normalize the heme in all already-prepped structures and re-score everything.
#
# PrepWizard produced two different heme states across the set (see README): most
# structures came out as the porphyrin dianion with deprotonated propionates
# (charge -2), but 8SPD kept an extra pyrrole N-H (charge -1). normalize_heme()
# in step 3 now enforces the -2 state for every structure.
#
# Step 3 skips PrepWizard where a prepped file already exists, so this is fast;
# step 4 reuses existing grids and only rebuilds those that were deleted.
#
#     sbatch scripts/rerun_normalize.sh

set -euo pipefail

module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1   # so tail -f on the log is actually useful

cd "$SLURM_SUBMIT_DIR"

$SCHRODINGER/run python3 scripts/03_prep_structures.py
$SCHRODINGER/run python3 scripts/04_grid_and_score.py

echo "=== Heme charges after normalization (all should be -2) ==="
$SCHRODINGER/run python3 -c "
import pandas as pd
d = pd.read_csv('results/prep_summary.csv')
print(d[['pdb_id','heme_fixed_propionates','heme_fixed_pyrrole_N','heme_formal_charge']].to_string(index=False))
print()
bad = d[d['heme_formal_charge'] != -2]
print('INCONSISTENT:' if len(bad) else 'All heme charges consistent (-2).')
print(bad[['pdb_id','heme_formal_charge']].to_string(index=False) if len(bad) else '')
"
