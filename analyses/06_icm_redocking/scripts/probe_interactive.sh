#!/bin/bash
#SBATCH --job-name=icm_inter
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-inter-%j.out
set -uo pipefail
module load icm/3.9.4
OB=~/micromamba/envs/obabel_env/bin/obabel
cd "$SLURM_SUBMIT_DIR"/structures/8SPD
[ -f 8SPD_crystal_v2000.sdf ] || $OB 8SPD_crystal.sdf -osdf -O 8SPD_crystal_v2000.sdf 2>/dev/null
echo "##### piping to icm64 stdin (interactive-like)"
icm64 <<'ICM' 2>&1 | grep -viE "^Startup|^ Hint|^ Info|^\.\.\.ICM|Warning" | tail -14
call _startup
read object "D_8SPD_rec.ob"
read mol "8SPD_crystal_v2000.sdf"
scanScoreExternal a_1. a_2.
quit
ICM
