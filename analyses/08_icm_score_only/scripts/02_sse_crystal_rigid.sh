#!/bin/bash
#SBATCH --job-name=sse_rigid
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-sserigid-%j.out
#
# Score the crystal poses with scanScoreExternal (merged-object recipe), which unlike
# `_dockScan -n` needs no docking project and returns the full component breakdown.
# Run over the whole affinity set so the rigid subset is not a special case.
set -uo pipefail
module load icm/3.9.4
OB=~/micromamba/envs/obabel_env/bin/obabel
cd "$SLURM_SUBMIT_DIR"
OUT=results/sse_crystal_raw.txt
: > $OUT
for d in ../06_icm_redocking/structures/*/; do
  PDB=$(basename "$d")
  [ -f "$d/${PDB}_rec.pdb" ] || continue
  [ -f "$d/${PDB}_crystal.sdf" ] || continue
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
echo "done: $(grep -c SCORES $OUT) structures scored"
