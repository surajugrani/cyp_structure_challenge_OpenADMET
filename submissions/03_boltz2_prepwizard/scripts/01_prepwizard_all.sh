#!/bin/bash
#SBATCH --job-name=pw_all20
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH -t 04:00:00
#SBATCH --output=submissions/03_boltz2_prepwizard/slurm-pwall-%j.out
# PrepWizard on all 20 Boltz-2 predictions, restrained to 0.3 A global heavy-atom RMSD.
#
# Validated on OCNT-2395328 first: closest contact 2.09 -> 2.81 A (PoseBusters PASS),
# backbone moved only 0.15 A median, ligand-Fe distance unchanged (4.40 -> 4.42 A),
# and chain/residue naming preserved (A/B/C, HEM 43, LIG 27, 4103/4103 atoms matched).
# Compare ICM self-tethered minimisation, which moved heavy atoms 1.0-7.6 A median and
# corrupted 6 ligands - not used.
#
# Serial, one at a time (~3 min each). Resumable, but note the work dir must be clean:
# a stale output from an earlier broken test silently caused a structure to be skipped
# during the ICM attempt.
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
IN=submissions/01_boltz2_default/submission
W=submissions/03_boltz2_prepwizard/work
mkdir -p "$W"
n=0
for f in "$IN"/*.pdb; do
  P=$(basename "$f" .pdb); n=$((n+1))
  if [ -s "$W/${P}_pw.pdb" ]; then echo "[$n/20] $P: already done"; continue; fi
  echo "[$n/20] $P: prepwizard"
  cp "$f" "$W/$P.pdb"
  ( cd "$W"
    $SCHRODINGER/utilities/prepwizard -NOJOBID -rehtreat -captermini \
      -propka_pH 7.4 -epik_pH 7.4 -f OPLS_2005 -rmsd 0.3 \
      "$P.pdb" "${P}_pw.pdb" 2>&1 | grep -E "Final rmsd|successful preparations|DONE|Error|error" )
done
echo "DONE prepped=$(ls $W/*_pw.pdb 2>/dev/null | wc -l)/20"
