#!/bin/bash
#SBATCH --job-name=a31_rep
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=analyses/31_msa_depth/slurm-%j.out
# RESOURCES: same RDKit-only comparison as run.sh, which finished in 10 seconds.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
PY=~/micromamba/envs/cyp_posebusters_env/bin/python
A=analyses/31_msa_depth
$PY $A/scripts/01_convert.py $A/replicate/raw $A/structures/replicate
# The control: same settings as submission 23, so any difference is pure run-to-run
# noise. Also against depth2048, to see whether the knob moves poses any further than
# a plain rerun does.
$PY tools/ligand_rmsd.py $A/structures/replicate submissions/23_esmfold2_default/submission \
    --label-a replicate --label-b sub23 --out $A/results/replicate_vs_sub23.csv
$PY tools/ligand_rmsd.py $A/structures/replicate $A/structures/depth2048 \
    --label-a replicate --label-b depth2048 --out $A/results/replicate_vs_depth2048.csv
echo "EXIT=$?"
