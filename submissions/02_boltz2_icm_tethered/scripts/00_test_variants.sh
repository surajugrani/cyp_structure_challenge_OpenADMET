#!/bin/bash
#SBATCH --job-name=icm_tvar
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:40:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-tvar-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
P=OCNT-2395328
for V in A B C; do
echo "############ VARIANT $V ############"
icm64 <<ICM
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_1.
set tether a_1.
if ("$V" == "A") minimize cartesian
if ("$V" == "B") then
  set terms "vw,14,hb,el,to,sf,tz"
  minimize cartesian "vw,14,hb,el,to,sf,tz"
endif
if ("$V" == "C") then
  set vsel v_//!?vt*
  set terms "vw,14,hb,el,to,sf,tz"
  minimize 2000 tzWeight=1.0
endif
write pdb a_1. "${P}_v${V}.pdb"
quit
ICM
done
