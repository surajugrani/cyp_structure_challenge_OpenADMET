#!/bin/bash
#SBATCH --job-name=cyp_fullredock
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=slurm-%j.out
#
# Analysis 05: prepare all 102 CYP3A4 structures and dock each ligand FROM SMILES
# (LigPrep), scoring pose recovery with spyrmsd against the crystal pose.
#
# SERIAL by design - only one Schrodinger licence token is available, so parallel array
# tasks would fail on licence checkout rather than finish sooner. ~100 structures at
# roughly 5 min each is a long run, hence the 36 h wall clock.
#
# The script is resumable and writes results after every structure, so if this job is
# interrupted just resubmit it and the run continues where it stopped.
#
#     sbatch scripts/run_all.sh                 # full set
#     sbatch scripts/run_all.sh --limit 2       # smoke test
#
# Progress:  tail -f slurm-<jobid>.out    |    wc -l results/redock_summary.csv
set -euo pipefail
module load schrodinger/2025-3
export LD_LIBRARY_PATH=$HOME:${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
# spyrmsd installed with --no-deps so it uses Schrodinger's own numpy/networkx
export PYTHONPATH=/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 scripts/02_prep_dock.py "$@"
