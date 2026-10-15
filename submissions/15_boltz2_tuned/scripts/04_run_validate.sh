#!/bin/bash
#SBATCH --job-name=pw_validate
#SBATCH -p shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/15_boltz2_tuned/slurm-validate-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/15_boltz2_tuned/scripts/03_graft_validate_pb.py
echo "EXIT=$?"
