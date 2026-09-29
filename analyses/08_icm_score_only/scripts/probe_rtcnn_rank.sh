#!/bin/bash
#SBATCH --job-name=rtcnn_rank
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-rank-%j.out
set -uo pipefail
module load icm/3.9.4
OB=~/micromamba/envs/obabel_env/bin/obabel
D="$SLURM_SUBMIT_DIR/../06_icm_redocking/structures/8SPD"
cd "$D"
# split the 10-pose hitlist into single-pose SDFs
$OB 8SPD_hits.sdf -osdf -O pose_.sdf -m 2>/dev/null
ls pose_*.sdf | head -20
for f in pose_*.sdf; do
  icm64 <<ICM 2>&1 | grep -E "^RT_"
call _startup
read object "D_8SPD_rec.ob"
read mol "$f"
convertObject a_2.
calcRTCNN a_1. a_2. ""
print "RT_$f= " + r_out
quit
ICM
done
