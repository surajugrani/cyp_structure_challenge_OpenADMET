#!/bin/bash
#SBATCH --job-name=icm_rtcnn
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-rtcnn-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"/../06_icm_redocking/structures/8SPD
icm64 <<'ICM' 2>&1 | grep -viE "^Startup|^ Hint|^\.\.\.ICM|Warning" | tail -20
call _startup
read object "D_8SPD_rec.ob"
read mol "8SPD_crystal_v2000.sdf"
convertObject a_2.
calcRTCNN a_1. a_2. ""
print "RTCNN_INPLACE= " + r_out
quit
ICM
