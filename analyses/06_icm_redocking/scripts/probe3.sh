#!/bin/bash
#SBATCH --job-name=icm_probe3
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-probe3-%j.out
set -euo pipefail
module load icm/3.9.4
OB=~/micromamba/envs/obabel_env/bin/obabel
cd "$SLURM_SUBMIT_DIR"
for PDB in 8SPD 6MA6; do
  cd structures/$PDB
  [ -f ${PDB}_crystal_v2000.sdf ] || $OB ${PDB}_crystal.sdf -osdf -O ${PDB}_crystal_v2000.sdf 2>/dev/null
  cat > /tmp/ss_$PDB.icm <<ICM
call _startup
read object "D_${PDB}_rec.ob"
read mol "${PDB}_crystal_v2000.sdf"
print "PDBID= ${PDB}"
scanScoreExternal a_1. a_2.
quit
ICM
  echo "########## $PDB"
  icm64 -s /tmp/ss_$PDB.icm 2>&1 | grep -viE "^Startup|^ Hint|^ Info|^\.\.\.ICM|Warning" | tail -8
  cd ../..
done
