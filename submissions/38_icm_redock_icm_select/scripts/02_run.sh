#!/bin/bash
#SBATCH --job-name=s38_exp
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 00:30:00
#SBATCH --output=submissions/38_icm_redock_icm_select/slurm-exp-%j.out
# RESOURCES: reads 20 cached complex.mae files and rewrites 20 PDBs, single-threaded,
#            no Glide or PrepWizard. Schrodinger's structure reader on a 4k-atom
#            complex sits well under 2G, so 3G and 30 min is the smallest ask that
#            covers it and backfills immediately.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 \
    submissions/38_icm_redock_icm_select/scripts/02_export_prepped.py
echo "EXIT=$?"
