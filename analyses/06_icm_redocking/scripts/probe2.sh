#!/bin/bash
#SBATCH --job-name=icm_probe2
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-probe2-%j.out
set -euo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"
cat > /tmp/p2_$$.icm <<'ICM'
call _startup
show macro scanScoreExternal
print "===SEP==="
show macro calcRTCNN
quit
ICM
icm64 -s /tmp/p2_$$.icm 2>&1 | grep -vE "^Startup|^ Hint|^ Info|^\.\.\.ICM"
