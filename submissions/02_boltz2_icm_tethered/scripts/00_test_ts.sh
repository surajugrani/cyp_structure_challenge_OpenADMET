#!/bin/bash
#SBATCH --job-name=icm_ts
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-ts-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
P=OCNT-2395328
# convertObject stores each atom's input coords as a selftether; the "ts" term
# restrains to them. mmff is tried separately because convertObject warned that
# a HEM atom has no MMFF type.
run () { # $1=tag $2=terms $3=extra
icm64 <<ICM
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_wrk.
set terms "$2"
minimize cartesian $3 2000 "$2" selftether=a_wrk_1.//!h*
write pdb a_wrk_1. "${P}_ts_$1.pdb"
quit
ICM
}
echo "########## A: icm terms ##########"; run icm "vw,14,hb,el,to,ss,ts" ""
echo "########## B: mmff      ##########"; run mmff "mmff,ts" "type charge"
