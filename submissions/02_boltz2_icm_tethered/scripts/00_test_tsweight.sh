#!/bin/bash
#SBATCH --job-name=icm_tsw
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 01:30:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-tsw-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
P=OCNT-2395328
for W in 1. 10. 100.; do
echo "########## TOOLS.tsWeight=$W ##########"
icm64 <<ICM
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_wrk.
TOOLS.tsWeight = $W
set terms "vw,14,hb,el,to,ss,ts"
minimize cartesian 2000 "vw,14,hb,el,to,ss,ts" selftether=a_wrk_1.//!h*
write pdb a_wrk_1. "${P}_w${W}pdb"
quit
ICM
done
