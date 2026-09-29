#!/bin/bash
#SBATCH --job-name=sse_7
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-sse7-%j.out
set -uo pipefail
module load icm/3.9.4
OB=~/micromamba/envs/obabel_env/bin/obabel
cd "$SLURM_SUBMIT_DIR"
OUT=results/sse_rigid7_raw.txt
: > $OUT
for PDB in 2V0M 8SPD 7LXL 3UA1 6MA6 6MA7 4D6Z; do
  d=../06_icm_redocking/structures/$PDB
  ( cd "$d"
    [ -f ${PDB}_crystal_v2000.sdf ] || $OB ${PDB}_crystal.sdf -osdf -O ${PDB}_crystal_v2000.sdf 2>/dev/null
    icm64 <<ICM 2>&1 | grep -E "^SCORES>" | sed "s/^/${PDB} /"
call _startup
read pdb "${PDB}_rec.pdb" name=rec
read mol "${PDB}_crystal_v2000.sdf" name=lig
move a_2. a_1.
convertObject a_1.
receptor = a_2. - a_2.m
scanScoreExternal receptor a_2.m
quit
ICM
  ) >> "$SLURM_SUBMIT_DIR/$OUT"
done
echo "scored: $(grep -c SCORES $OUT)/7"
