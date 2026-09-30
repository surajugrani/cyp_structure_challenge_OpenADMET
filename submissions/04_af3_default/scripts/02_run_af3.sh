#!/bin/bash
#SBATCH --job-name=af3_cyp
#SBATCH --partition=alphafold
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=20:00:00
#SBATCH --exclude=nodeb201
#SBATCH --output=submissions/04_af3_default/slurm-af3-%j.out
# AlphaFold3, default settings, all 20 blinded ligands.
#
# Copied from the PXR recipe (slurm/run_af3.sh) - same partition, gres, cpus, walltime,
# module, TF_FORCE_UNIFIED_MEMORY, and `python $(which run_alphafold.py)` with
# --input_dir/--output_dir. AF3 iterates the directory itself, one complex at a time.
#
# Three deliberate deviations, each with a reason:
#   (mem stays at PXR's 16G: complexes run one at a time so peak memory is set by a
#    single 503 aa complex, not by the 20, and a larger request only makes the job harder
#    to backfill and slower to start.)
#   --exclude=nodeb201: that node cannot see /gpfs/home (autofs stub), so a job landing
#     there dies instantly with RaisedSignal:53 before writing a log. It is often the only
#     idle node in this partition so the scheduler picks it first. nodeb217/311/601 work.
#   --run_data_pipeline=false: every input already carries unpairedMsaPath (one shared
#     ColabFold MSA, reused from submission 01) and templates: [], so there is nothing for
#     the pipeline to do. This is what keeps the MSA from being searched 20 times.
set -uo pipefail
module load alphafold/3.0
export TF_FORCE_UNIFIED_MEMORY=1
cd "$SLURM_SUBMIT_DIR"
D=submissions/04_af3_default
python "$(which run_alphafold.py)" \
  --input_dir="$D/inputs" \
  --output_dir="$D/results" \
  --run_data_pipeline=false \
  --run_inference=true
echo "EXIT=$?"
