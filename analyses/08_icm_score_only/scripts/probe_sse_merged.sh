#!/bin/bash
#SBATCH --job-name=sse_merged
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-sse-%j.out
set -uo pipefail
module load icm/3.9.4
D="$SLURM_SUBMIT_DIR/../06_icm_redocking/structures/8SPD"
cd "$D"
icm64 <<'ICM' 2>&1 | grep -viE "^Startup|^ Hint|^\.\.\.ICM|Warning" | tail -25
call _startup
read pdb "8SPD_rec.pdb" name=rec
read mol "8SPD_crystal_v2000.sdf" name=lig
move a_2. a_1.
convertObject a_1.
list object
receptor = a_2. - a_2.m
print "SCORING_NOW"
scanScoreExternal receptor a_2.m
quit
ICM
