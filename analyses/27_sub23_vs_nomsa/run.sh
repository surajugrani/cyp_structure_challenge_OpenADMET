#!/bin/bash
#SBATCH --job-name=a27_rmsd
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH --output=analyses/27_sub23_vs_nomsa/slurm-%j.out
# RESOURCES: 1 serial process -> 1 CPU; RDKit + numpy over 20 structures, measured at
#            2-4G for this work -> 4G; the same comparison for 19 vs 24 ran in well
#            under a minute -> 20 min is generous.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python tools/ligand_rmsd.py \
  submissions/23_esmfold2_default/submission \
  submissions/27_esmfold2_nomsa/submission \
  --label-a 23 --label-b 27 \
  --out analyses/27_sub23_vs_nomsa/results/sub23_vs_27.csv
echo "EXIT=$?"
