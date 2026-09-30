#!/bin/bash
#SBATCH --job-name=pw_validate
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/03_boltz2_prepwizard/slurm-validate-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/03_boltz2_prepwizard/scripts/02_graft_validate_pb.py
echo "EXIT=$?"
