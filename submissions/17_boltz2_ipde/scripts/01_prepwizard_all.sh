#!/bin/bash
#SBATCH --job-name=pw17
#SBATCH -p shared
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/17_boltz2_ipde/slurm-pwall-%j.out
# PrepWizard on the 20 tuned-Boltz model_0 predictions, restrained to 0.3 A global
# heavy-atom RMSD - identical options to submission 03, so 15 vs 03 isolates the Boltz
# sampling settings and nothing else.
#
# TEN AT A TIME.
# allows this many concurrent PrepWizard runs; serial was ~4.5 min each (80 min total).
# single job with background subshells, not an array job.
#
# Resumable: a structure with a non-empty *_pw.pdb is skipped, so the two already
# finished by the cancelled serial run are kept.
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
IN=submissions/17_boltz2_ipde/raw_sel
W=submissions/17_boltz2_ipde/work
mkdir -p "$W"
MAXJOBS=6

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
