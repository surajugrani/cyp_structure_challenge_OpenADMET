#!/bin/bash
#SBATCH --job-name=s12_icm
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH --output=submissions/12_icm_redock_af3/slurm-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/12_icm_redock_af3/scripts/01_build.py "$@"
echo "EXIT=$?"
