#!/bin/bash
#SBATCH --job-name=s26_maxacc
#SBATCH -p shared
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/26_esmfold2_maxacc/slurm-%j.out
# Three passes, same pipeline as submission 23: mmCIF -> PDB, PrepWizard, then convert
# again so the minimised coordinates are grafted in, validate and zip.
# RESOURCES: 5 concurrent PrepWizard runs -> 5 CPUs (user caps at 5); PrepWizard
#            measured <2G each -> 10G; ~4.5 min each, 20 structures 5-wide -> ~20 min
#            plus two cheap conversion passes, 2h margin.
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
PB=~/micromamba/envs/cyp_posebusters_env/bin/python
S=submissions/26_esmfold2_maxacc
W=$S/work
mkdir -p "$W"

echo "=== pass 1: mmCIF -> PDB (no graft yet) ==="
$PB $S/scripts/01_convert_validate.py

echo "=== pass 2: PrepWizard, 5 concurrent ==="
prep() {
  local f=$1 P
  P=$(basename "$f" .pdb)
  if [ -s "$W/${P}_pw.pdb" ]; then echo "$P: already done"; return; fi
  cp "$f" "$W/$P.pdb"
  ( cd "$W"
    $SCHRODINGER/utilities/prepwizard -NOJOBID -rehtreat -captermini \
      -propka_pH 7.4 -epik_pH 7.4 -f OPLS_2005 -rmsd 0.3 \
      "$P.pdb" "${P}_pw.pdb" >"${P}_pw.log" 2>&1 )
  if [ -s "$W/${P}_pw.pdb" ]; then echo "$P: ok"; else echo "$P: FAILED"; fi
}
for f in "$S"/submission/*.pdb; do
  while [ "$(jobs -rp | wc -l)" -ge 5 ]; do wait -n; done
  prep "$f" &
done
wait
echo "prepped=$(ls $W/*_pw.pdb 2>/dev/null | wc -l)/20"

echo "=== pass 3: graft minimised coordinates, validate, zip ==="
$PB $S/scripts/01_convert_validate.py
echo "EXIT=$?"
# Always build the PyMOL session alongside the zip (no-op if no zip was written).
bash tools/make_pse.sh submissions/26_esmfold2_maxacc
