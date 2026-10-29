#!/bin/bash
#SBATCH --job-name=s34_build
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH --output=submissions/34_glide_minimised/slurm-%j.out
# RESOURCES: reads 20 existing Glide pose-viewer files and rewrites 20 PDBs. No docking,
#            no scoring, no minimisation - the poses were produced by analysis 28 and are
#            already on disk. 1 CPU, 4G and 20 min is well over what a file rewrite needs;
#            the small ask backfills immediately.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 submissions/34_glide_minimised/scripts/01_build.py
echo "EXIT=$?"
