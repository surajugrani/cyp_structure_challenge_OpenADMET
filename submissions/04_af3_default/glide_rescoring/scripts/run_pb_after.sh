#!/bin/bash
#SBATCH --job-name=pb_after
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH --output=submissions/04_af3_default/glide_rescoring/slurm-pbafter-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/04_af3_default/glide_rescoring/scripts/03_pb_after_prepwizard.py
echo "EXIT=$?"
