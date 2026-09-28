#!/bin/bash
#SBATCH --job-name=icm_probe
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-probe-%j.out
set -euo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"
cat > /tmp/probe_$$.icm <<'ICM'
call _startup
help scanScoreExternal
print "---SEP---"
help RTCNN
print "---SEP2---"
show macro scanScoreExternal
quit
ICM
icm64 -s /tmp/probe_$$.icm 2>&1 | grep -vE "^Startup|^ Hint|^ Info|^\.\.\.ICM"
