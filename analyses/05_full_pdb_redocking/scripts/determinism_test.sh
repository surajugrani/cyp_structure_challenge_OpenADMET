#!/bin/bash
#SBATCH --job-name=cyp_determ
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-determ-%j.out
#
# Is Glide SP confgen deterministic? Analysis 05 reproduced 6MA7 exactly but not 4D6Z,
# from inputs that differ only by a title string. Re-dock the same ligand into the same
# grid N times and compare - this sets how much run-to-run noise sits under every pose
# recovery number we report.
set -euo pipefail
module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"

for PDB in 4D6Z 6MA7; do
  D=structures/$PDB
  for i in 1 2 3; do
    sed "s/^GRIDFILE.*/GRIDFILE   ${PDB}_grid-gridgen.zip/" $D/${PDB}_redock.in > $D/${PDB}_det${i}.in
    ( cd $D && $SCHRODINGER/glide -NOJOBID -WAIT ${PDB}_det${i}.in >/dev/null 2>&1 || true )
  done
done

$SCHRODINGER/run python3 - <<'PY'
import numpy as np, os
from schrodinger import structure
for pdb in ["4D6Z","6MA7"]:
    d=f"structures/{pdb}"
    ref=np.array([a.xyz for a in structure.StructureReader.read(f"{d}/{pdb}_lig.sdf").atom
                  if a.element!="H"])
    print(f"\n{pdb}:")
    for i in (1,2,3):
        pv=f"{d}/{pdb}_det{i}_pv.mae"
        if not os.path.exists(pv):
            print(f"  run {i}: no output"); continue
        r=[]
        for k,st in enumerate(structure.StructureReader(pv)):
            if k==0: continue
            r.append(st.property.get("r_i_glide_rmsd_to_input"))
        print(f"  run {i}: n={len(r)} top1={r[0]:.3f} best={min(r):.3f} "
              f"all={[round(x,2) for x in r]}")
PY
