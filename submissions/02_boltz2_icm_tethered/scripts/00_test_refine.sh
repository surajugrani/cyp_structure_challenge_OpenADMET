#!/bin/bash
#SBATCH --job-name=icm_refine
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-refine-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
P=OCNT-2395328
# sc  = side chains only (i_nRegIter=0) - leaves the Boltz backbone untouched
# full= 5 annealing iterations + side-chain refinement
for V in sc full; do
  [ "$V" = sc ] && ARGS="0 yes" || ARGS="5 yes"
  echo "############ refineModel $ARGS ############"
  icm64 <<ICM
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_wrk.
refineModel $ARGS
write pdb a_wrk_1. "${P}_ref_${V}.pdb"
quit
ICM
done
