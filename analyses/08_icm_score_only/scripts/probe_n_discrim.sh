#!/bin/bash
#SBATCH --job-name=n_discrim
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-discrim-%j.out
set -uo pipefail
module load icm/3.9.4
D="$SLURM_SUBMIT_DIR/../06_icm_redocking/structures/8SPD"
cd "$D"
for i in 1 5 10; do
  rm -f nchk_$i.sdf
  icm64 /opt/applications/icm/3.9.4/_dockScan D_8SPD input=pose_$i.sdf \
     confs=1 effort=1. seed=42 output=nchk_$i.sdf -n -s -a -S >/dev/null 2>&1 || true
  v=$(awk '/^> *<RTCNN>/{getline; print}' nchk_$i.sdf 2>/dev/null)
  s=$(awk '/^> *<Score>/{getline; print}' nchk_$i.sdf 2>/dev/null)
  echo "pose_$i -> RTCNN=${v:-NONE} Score=${s:-NONE}"
done
