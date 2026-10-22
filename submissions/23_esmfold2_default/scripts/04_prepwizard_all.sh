#!/bin/bash
#SBATCH --job-name=pw23
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 03:00:00
#SBATCH --output=submissions/23_esmfold2_default/slurm-pw-%j.out
# RESOURCES: serial - 1 PrepWizard at a time -> 1 CPU; measured <2G each -> 3G;
#            ~4.5 min each x 20 serial -> ~90 min, 3h margin. Smallest possible ask,\n#            to backfill soonest on a busy queue.
#
# Same -rmsd 0.3 restraint as every other submission, so the ESMFold2 set is prepared
# identically to the Boltz and AF3 ones. PrepWizard is the known fix for the three
# PoseBusters failures here (two minimum_distance_to_protein clashes and one
# non-aromatic ring distortion): it took submission 01 from 17/20 to 20/20 and the AF3
# samples from 57/100 to 95/100, moving backbone only ~0.15 A median.
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
IN=submissions/23_esmfold2_default/submission
W=submissions/23_esmfold2_default/work
mkdir -p "$W"
MAXJOBS=1

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

for f in "$IN"/*.pdb; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done
  prep "$f" &
done
wait
echo "DONE prepped=$(ls $W/*_pw.pdb 2>/dev/null | wc -l)/20"
