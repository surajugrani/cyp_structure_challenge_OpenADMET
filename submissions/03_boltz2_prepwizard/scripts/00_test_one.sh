#!/bin/bash
#SBATCH --job-name=pw_test
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/03_boltz2_prepwizard/slurm-pwtest-%j.out
# PrepWizard on one PoseBusters failure. Options match analysis 01's proven call
# (-f OPLS_2005, not "-f 2005", which schrodinger/2025-3 rejects).
# -rmsd 0.3 is a GLOBAL heavy-atom RMSD convergence threshold, so individual clashing
# side-chain atoms can move further than 0.3 A while the average stays restrained.
# -noepik / -noprotassign are NOT used: default H-bond optimisation is part of what
# relieves the clash. -nometaltreat left default so the haem Fe is not re-typed oddly.
set -uo pipefail
module load schrodinger/2025-3
cd "$SLURM_SUBMIT_DIR"
W=submissions/03_boltz2_prepwizard/work
P=OCNT-2395328
cp submissions/01_boltz2_default/submission/$P.pdb $W/$P.pdb
cd $W
$SCHRODINGER/utilities/prepwizard -NOJOBID -rehtreat -captermini \
  -propka_pH 7.4 -epik_pH 7.4 -f OPLS_2005 -rmsd 0.3 \
  $P.pdb ${P}_pw.pdb
echo "EXIT=$?"
ls -la ${P}_pw.pdb 2>/dev/null || echo "NO OUTPUT"
