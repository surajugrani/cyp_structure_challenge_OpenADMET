#!/bin/bash
#SBATCH --job-name=cyp_ligprep
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=slurm-ligprep-%j.out
#
# Step 0: prepare every ligand from SMILES, up front. Shared input to the Glide run
# (analysis 05) and the ICM run (analysis 06), so both can then proceed independently.
set -euo pipefail
module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/00_ligprep_all.py
