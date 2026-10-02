#!/bin/bash
#SBATCH --job-name=s08_prep
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 12:00:00
#SBATCH --output=submissions/08_af3_1template/slurm-prep-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module load schrodinger/2025-3
$SCHRODINGER/run python3 submissions/08_af3_1template/scripts/04_convert_and_prep.py
echo "EXIT=$?"
