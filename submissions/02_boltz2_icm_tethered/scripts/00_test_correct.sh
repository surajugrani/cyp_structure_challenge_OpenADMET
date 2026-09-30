#!/bin/bash
#SBATCH --job-name=icm_tfix
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:40:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-tfix-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
P=OCNT-2395328
# w = tether weight; higher = stays closer to the Boltz pose
for W in 0.1 1.0 10.0; do
echo "############ tzWeight=$W ############"
icm64 <<ICM
call _startup
read pdb "$P.pdb" name="wrk"
convertObject a_wrk.
set tether a_wrk_1.
set terms "vw,14,hb,el,to,sf,tz"
minimize cartesian "vw,14,hb,el,to,sf,tz" tzWeight=$W
write pdb a_wrk_1. "${P}_tz${W}.pdb"
quit
ICM
done
