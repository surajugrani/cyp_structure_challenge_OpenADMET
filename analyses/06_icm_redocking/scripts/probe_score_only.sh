#!/bin/bash
#SBATCH --job-name=icm_scoreonly
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-scoreonly-%j.out
set -euo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"/structures/8SPD
cat > /tmp/so_$$.icm <<'ICM'
call _startup
read object "D_8SPD_rec.ob"
read mol "8SPD_crystal_v2000.sdf"
print "OBJECTS:"
list object
print "---RUN scanScoreExternal---"
scanScoreExternal a_1. a_2.
print "SCORE_VAR= " + Score
print "---RUN calcRTCNN---"
calcRTCNN a_1. a_2. ""
quit
ICM
# crystal ligand must be V2000 for ICM
if [ ! -f 8SPD_crystal_v2000.sdf ]; then
  ~/micromamba/envs/obabel_env/bin/obabel 8SPD_crystal.sdf -osdf -O 8SPD_crystal_v2000.sdf 2>/dev/null
fi
icm64 -s /tmp/so_$$.icm 2>&1 | grep -vE "^Startup|^ Hint|^\.\.\.ICM" | tail -45
