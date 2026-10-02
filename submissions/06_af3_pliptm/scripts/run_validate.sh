#!/bin/bash
#SBATCH --job-name=s06_val
#SBATCH -p shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/06_af3_pliptm/slurm-val-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/06_af3_pliptm/scripts/02_validate_pb.py
echo "EXIT=$?"
