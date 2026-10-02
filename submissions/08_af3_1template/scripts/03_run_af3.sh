#!/bin/bash
#SBATCH --job-name=af3_1t
#SBATCH --partition=alphafold
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --exclude=nodeb201
#SBATCH --output=submissions/08_af3_1template/slurm-af3-%j.out
# AF3 with the 38FY + 38FU cryo-EM templates. Same recipe as submission 04 otherwise:
# PXR's run_af3.sh settings, shared ColabFold MSA, --run_data_pipeline=false.
# nodeb201 excluded - it cannot see /gpfs/home and dies with RaisedSignal:53.
set -uo pipefail
module load alphafold/3.0
export TF_FORCE_UNIFIED_MEMORY=1
cd "$SLURM_SUBMIT_DIR"
D=submissions/08_af3_1template
python "$(which run_alphafold.py)" \
  --input_dir="$D/inputs" \
  --output_dir="$D/results" \
  --run_data_pipeline=false \
  --run_inference=true
echo "EXIT=$?"
