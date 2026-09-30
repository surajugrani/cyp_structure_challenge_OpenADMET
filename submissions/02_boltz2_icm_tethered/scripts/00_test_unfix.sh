#!/bin/bash
#SBATCH --job-name=icm_unfix
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-unfix-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
P=OCNT-2395328
echo "############ unfix chi + refineModel 0 yes ############"
icm64 <<ICM
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_wrk.
print "free vars before unfix: ", Nof(v_//!?vt* & a_wrk_1.)
unfix v_//x*
print "free vars after unfix : ", Nof(v_//x* & a_wrk_1.)
refineModel 0 yes
write pdb a_wrk_1. "${P}_unfix_sc.pdb"
quit
ICM
