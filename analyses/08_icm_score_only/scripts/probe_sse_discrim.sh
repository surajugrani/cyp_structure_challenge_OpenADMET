#!/bin/bash
#SBATCH --job-name=sse_discrim
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=slurm-ssed-%j.out
set -uo pipefail
module load icm/3.9.4
D="$SLURM_SUBMIT_DIR/../06_icm_redocking/structures/8SPD"
cd "$D"
for i in 1 5 10; do
  icm64 <<ICM 2>&1 | grep -E "^SCORES>" | sed "s/^/pose_$i /"
call _startup
read pdb "8SPD_rec.pdb" name=rec
read mol "pose_$i.sdf" name=lig
move a_2. a_1.
convertObject a_1.
receptor = a_2. - a_2.m
scanScoreExternal receptor a_2.m
quit
ICM
done
