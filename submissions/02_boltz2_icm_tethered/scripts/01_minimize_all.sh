#!/bin/bash
#SBATCH --job-name=icm_min20
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-min20-%j.out
# ICM self-tethered Cartesian minimisation of all 20 Boltz-2 predictions.
#
# Recipe, from man/icm/icm-commands.html and man/icm/terms.html:
#   convertObject stores every atom's input coordinates as a SELFTETHER, so no
#   `set tether` is needed - the "ts" term restrains to them. ("tz" is for tethers
#   between two objects and is what made four earlier attempts fail.)
#   `minimize cartesian ... selftether=as_` chooses which atoms are tethered.
#   "mmff" is NOT a valid term in icm/3.9-4 ("unrecognized term") and silently falls
#   back to defaults, which moved the backbone 2.3 A - so ICM's own terms are used.
#   TOOLS.tsWeight is left at its default of 1.0: on the test structure that fixed the
#   clash (2.09 -> 2.78 A) for 0.89 A median backbone movement, which was accepted.
#
# Serial, one structure at a time. Resumable: existing outputs are skipped.
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"
IN=submissions/01_boltz2_default/submission
W=submissions/02_boltz2_icm_tethered/work
mkdir -p "$W"
n=0
for f in "$IN"/*.pdb; do
  P=$(basename "$f" .pdb); n=$((n+1))
  if [ -s "$W/${P}_min.pdb" ]; then echo "[$n/20] $P: already done"; continue; fi
  echo "[$n/20] $P: minimizing"
  cp "$f" "$W/$P.pdb"
  ( cd "$W"
    icm64 <<ICM 2>&1 | grep -E "Error>|Minimizing|atoms written"
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_wrk.
set terms "vw,14,hb,el,to,ss,ts"
minimize cartesian 2000 "vw,14,hb,el,to,ss,ts" selftether=a_wrk_1.//!h*
write pdb a_wrk_1. "${P}_min.pdb"
quit
ICM
  )
done
echo "DONE minimized=$(ls $W/*_min.pdb 2>/dev/null | wc -l)/20"
