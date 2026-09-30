#!/bin/bash
#SBATCH --job-name=icm_objs
#SBATCH -p shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH -t 00:20:00
#SBATCH --output=submissions/02_boltz2_icm_tethered/slurm-objs-%j.out
set -uo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR/submissions/02_boltz2_icm_tethered/work"
icm64 <<'ICM'
call _startup
read pdb "OCNT-2395328.pdb" name="wrk"
print "=== objects after read ==="
list objects
convertObject a_wrk.
print "=== objects after convertObject ==="
list objects
ICM
