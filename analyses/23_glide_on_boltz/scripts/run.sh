#!/bin/bash
#SBATCH --job-name=a23_glide
#SBATCH --partition=shared
#SBATCH --cpus-per-task=5
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH --output=analyses/23_glide_on_boltz/slurm-%j.out
# 100 Boltz samples: PrepWizard -> grid -> mininplace score.
# Resources derived: 5 concurrent Glide pipelines (user-approved licence limit) = 5 CPUs;
# PrepWizard/Glide peak ~2G each, so 5 x 2G rounded up = 12G. Submission 04's serial run
# of the same 100-sample pipeline took ~9 h, so 5 concurrent needs ~2 h; 8 h is the
# safety margin on a resumable job.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/23_glide_on_boltz/scripts/01_glide_score.py "$@"
echo "EXIT=$?"
