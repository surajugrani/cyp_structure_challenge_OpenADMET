#!/bin/bash
#SBATCH --job-name=s37_pse
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/37_icm_redock_glide_select/slurm-pse-%j.out
# RESOURCES: PyMOL loading 20 structures and writing one session, single-threaded.
#            Prior sessions on this set took ~1 min at under 2G; 4G and 30 min is the
#            smallest ask that covers it and backfills immediately.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
bash tools/make_pse.sh submissions/37_icm_redock_glide_select
echo "EXIT=$?"
