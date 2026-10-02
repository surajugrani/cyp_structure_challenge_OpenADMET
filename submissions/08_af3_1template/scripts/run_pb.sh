#!/bin/bash
#SBATCH --job-name=s08_pb
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 04:00:00
#SBATCH --output=submissions/08_af3_1template/slurm-pb-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/08_af3_1template/scripts/05_pb_select_package.py
echo "EXIT=$?"
