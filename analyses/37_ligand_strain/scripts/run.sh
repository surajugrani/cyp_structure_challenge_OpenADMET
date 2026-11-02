#!/bin/bash
#SBATCH --job-name=a37_strain
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --output=analyses/37_ligand_strain/slurm-%j.out
# RESOURCES: RDKit and PoseBusters over 60 structures, single-threaded. Validation
#            passes on 20 structures took ~3 min at under 2G, so 60 is ~9 min; 4G and
#            1 h is the smallest ask that covers conformer generation and backfills now.
set -uo pipefail
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
    analyses/37_ligand_strain/scripts/01_strain.py
echo "EXIT=$?"
