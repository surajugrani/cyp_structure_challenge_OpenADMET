#!/bin/bash
#SBATCH --job-name=icm_probe5
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-probe5-%j.out
set -uo pipefail
module load icm/3.9.4
OB=~/micromamba/envs/obabel_env/bin/obabel
cd "$SLURM_SUBMIT_DIR"/structures/8SPD
[ -f 8SPD_crystal_v2000.sdf ] || $OB 8SPD_crystal.sdf -osdf -O 8SPD_crystal_v2000.sdf 2>/dev/null
rm -f 8SPD_native_hits.sdf
icm64 /opt/applications/icm/3.9.4/_dockScan D_8SPD input=8SPD_crystal_v2000.sdf \
   confs=1 effort=1. seed=42 output=8SPD_native_hits.sdf -n -s -a -S > /tmp/n.log 2>&1 || true
echo "=== exists:"; ls -l 8SPD_native_hits.sdf 2>/dev/null || echo MISSING
echo "=== tags:"; grep -oE '^> *<(Score|RTCNN|mfScore|NAME)>' 8SPD_native_hits.sdf 2>/dev/null | sort -u
echo "=== values:"; awk '/^> *<(Score|RTCNN|mfScore)>/{t=$0; getline v; print t" "v}' 8SPD_native_hits.sdf 2>/dev/null
echo "=== log tail:"; tail -6 /tmp/n.log
