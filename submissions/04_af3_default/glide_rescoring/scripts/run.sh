#!/bin/bash
#SBATCH --job-name=af3_glide
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH --output=submissions/04_af3_default/glide_rescoring/slurm-%j.out
# PrepWizard + grid + Glide mininplace for all 100 AF3 samples. ~6 min each, serial,
# resumable (samples already in results/glide_scores.csv are skipped).
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module load schrodinger/2025-3
$SCHRODINGER/run python3 \
  submissions/04_af3_default/glide_rescoring/scripts/01_prep_grid_score.py
echo "EXIT=$?"
