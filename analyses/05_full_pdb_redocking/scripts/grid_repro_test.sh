#!/bin/bash
#SBATCH --job-name=cyp_gridrepro
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-gridrepro-%j.out
#
# Glide docking is deterministic (verified), yet analyses 02 and 05 gave different poses
# for 4D6Z from identical receptor/ligand coordinates and identical gridgen.in.
# Question: is GRID GENERATION reproducible? Rebuild the grid from the same inputs and
# re-dock; if the answer changes, grid generation is the source of the variation.
set -euo pipefail
module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"/structures/4D6Z

for i in 1 2; do
  $SCHRODINGER/utilities/generate_glide_grids -NOJOBID -WAIT \
     -rec_file 4D6Z_reclig.mae -lig_asl "res. LIG" -j 4D6Z_rg${i} >/dev/null 2>&1 || true
  sed -e "s/^GRIDFILE.*/GRIDFILE   4D6Z_rg${i}-gridgen.zip/" 4D6Z_redock.in > 4D6Z_rg${i}_dock.in
  $SCHRODINGER/glide -NOJOBID -WAIT 4D6Z_rg${i}_dock.in >/dev/null 2>&1 || true
done

cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 - <<'PY'
import os
from schrodinger import structure
d="structures/4D6Z"
print("grid zip sizes:")
for f in sorted(os.listdir(d)):
    if f.endswith("-gridgen.zip"): print(f"  {f}: {os.path.getsize(d+'/'+f)}")
print("\npose RMSDs:")
for tag in ["redock","rg1_dock","rg2_dock"]:
    pv=f"{d}/4D6Z_{tag}_pv.mae"
    if not os.path.exists(pv): print(f"  {tag}: missing"); continue
    r=[st.property.get("r_i_glide_rmsd_to_input")
       for k,st in enumerate(structure.StructureReader(pv)) if k>0]
    print(f"  {tag:10s} top1={r[0]:.3f} best={min(r):.3f} all={[round(x,2) for x in r]}")
PY
