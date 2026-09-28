#!/bin/bash
#SBATCH --job-name=icm_conv
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-conv-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"/structures/8SPD
icm64 <<'ICM' 2>&1 | grep -viE "^Startup|^ Hint|^\.\.\.ICM|Warning" | tail -12
call _startup
read object "D_8SPD_rec.ob"
read mol "8SPD_crystal_v2000.sdf"
convertObject a_2.
scanScoreExternal a_1. a_2.
quit
ICM
