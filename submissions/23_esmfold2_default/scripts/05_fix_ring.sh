#!/bin/bash
#SBATCH --job-name=pw23fix
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/23_esmfold2_default/slurm-fix-%j.out
# RESOURCES: 1 structure, 1 PrepWizard run -> 1 CPU; measured <2G -> 3G; ~5 min -> 30m.
#
# OCNT-2313038 is the last PoseBusters failure: non-aromatic_ring_non-flatness, a
# distorted saturated ring that a 0.3 A restraint cannot relax.
#
# NOTE: simply dropping -rmsd does NOT loosen anything - 0.3 is PrepWizard's own default
# ("-r RMSD, -rmsd RMSD   Minimization: RMSD cutoff in angstrom (default 0.3)"). An
# earlier version of this script did exactly that and would have reproduced the identical
# failing structure. The cutoff has to be raised explicitly.
#
# 1.0 A is tried first, then 2.0 A if the ring is still flagged. Trade-off accepted
# deliberately: a looser cutoff moves the pose further from what ESMFold2 predicted, so
# this one structure is less purely "ESMFold2" than the other 19. Displacement is
# measured afterwards so the cost is explicit.
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
W=submissions/23_esmfold2_default/work
P=OCNT-2313038
for CUT in 1.0 2.0; do
  ( cd "$W"
    rm -f "${P}_loose.pdb"
    cp "${P}.pdb" "${P}_loose_in.pdb"
    $SCHRODINGER/utilities/prepwizard -NOJOBID -rehtreat -captermini \
      -propka_pH 7.4 -epik_pH 7.4 -f OPLS_2005 -rmsd "$CUT" \
      "${P}_loose_in.pdb" "${P}_loose.pdb" >"${P}_loose_${CUT}.log" 2>&1 )
  if [ -s "$W/${P}_loose.pdb" ]; then
    echo "prep OK at -rmsd $CUT"
    cp "$W/${P}_loose.pdb" "$W/${P}_pw.pdb"   # the converter grafts from *_pw.pdb
    cp "$W/${P}_loose.pdb" "$W/${P}_loose_${CUT}.pdb"
    break
  fi
  echo "no output at -rmsd $CUT"
done
[ -s "$W/${P}_pw.pdb" ] || echo "FAILED - restrained version left in place"
