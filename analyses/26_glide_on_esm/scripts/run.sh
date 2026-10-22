#!/bin/bash
#SBATCH --job-name=a26_esm
#SBATCH --partition=shared
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH -t 03:00:00
#SBATCH --output=analyses/26_glide_on_esm/slurm-%j.out
# 20 ESMFold2 structures: PrepWizard (-> .mae) -> grid -> mininplace score.
# PrepWizard is back in the loop because the receptor must reach Glide as .mae; a .pdb
# loses bond orders and gridgen dies in Lewis-structure perception.
# RESOURCES: 5 concurrent pipelines -> 5 CPUs (user caps at 5); PrepWizard and Glide
#            both measured <2G each -> 10G; PrepWizard ~4.5 min + grid/score ~3 min,
#            20 structures 5-wide -> ~35 min, 3h margin on a resumable job.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/26_glide_on_esm/scripts/01_glide_score.py "$@"
echo "EXIT=$?"
