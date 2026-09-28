#!/bin/bash
#SBATCH --job-name=icm_probe4
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-probe4-%j.out
set -euo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"/structures/8SPD
# variant A: convert the ligand object first
cat > /tmp/vA.icm <<'ICM'
call _startup
read object "D_8SPD_rec.ob"
read mol "8SPD_crystal_v2000.sdf"
convertObject a_2. no no no no no no no ""
print "VARIANT_A"
scanScoreExternal a_1. a_2.
quit
ICM
# variant B: merge ligand into the receptor object, then score mol selections
cat > /tmp/vB.icm <<'ICM'
call _startup
read object "D_8SPD_rec.ob"
read mol "8SPD_crystal_v2000.sdf"
move a_2.1 a_1.
print "VARIANT_B nmol=" + Nof(a_1.*)
scanScoreExternal a_1.1 a_1.2
quit
ICM
for v in A B; do echo "##### VARIANT $v"; icm64 -s /tmp/v$v.icm 2>&1 | grep -viE "^Startup|^ Hint|^ Info|^\.\.\.ICM|Warning" | tail -8; done
