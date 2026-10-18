#!/bin/bash
#SBATCH --job-name=pw_validate
#SBATCH -p shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/18_boltz2_glide/slurm-validate-%j.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/18_boltz2_glide/scripts/01_graft_validate_pb.py
echo "EXIT=$?"

# Always build the PyMOL session alongside the zip (no-op if no zip was written).
bash tools/make_pse.sh submissions/18_boltz2_glide
