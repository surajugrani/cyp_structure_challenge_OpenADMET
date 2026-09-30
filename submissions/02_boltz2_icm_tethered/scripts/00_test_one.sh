#!/bin/bash
#SBATCH --job-name=icm_tether_test
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-test-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"
W=submissions/02_boltz2_icm_tethered/work
P=OCNT-2395328          # one of the three PoseBusters failures
cp submissions/01_boltz2_default/submission/$P.pdb $W/$P.pdb
cd $W
icm64 <<ICM
call _startup
read pdb "$P.pdb" name="ref"
read pdb "$P.pdb" name="wrk"
convertObject a_2.
set terms "vw,14,hb,to,el,sf,tz"
set tether a_2. a_1.//!h*
minimize 1000 "vw,14,hb,to,el,sf,tz" tzWeight=0.5
write pdb a_2. "${P}_min.pdb"
quit
ICM
