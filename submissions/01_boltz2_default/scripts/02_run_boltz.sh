#!/bin/bash
#SBATCH --job-name=boltz_cyp01
#SBATCH -p gpu
#SBATCH --gres=gpu:gtx1080ti:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --output=submissions/01_boltz2_default/slurm-boltz-%j.out
# Boltz-2, default settings, all 20 blinded ligands, each with the haem cofactor (CCD HEM).
#
# The env's binary is called directly rather than via `micromamba activate`: the
# micromamba shell function is defined in ~/.bashrc, which a batch script does not source,
# so activation fails with "command not found" under sbatch.
#
# gtx1080ti (11 GB) is requested rather than the generic gpu:1, which can land on an 8 GB
# GTX 1080. The rtxpro6000 nodes would be much faster but both are currently drained or
# reserved, so a job requesting them never starts.
#
# --use_msa_server generates MSAs remotely, so no local a3m is needed.
# --output_format pdb, NOT mmcif. The .cif files in the tutorial repo's
# inputs/boltz-2_cyp3a4_predictions/ are example *inputs*, not example submissions:
# validation/structure_validation.py accepts only files ending .pdb and rejects a zip with
# "Zip file contains no PDB files." Those .cif files did settle the chain layout though -
# protein, HEM, ligand - which is what these YAMLs produce.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
# torch in this env is built against CUDA 11.8; without the matching module the runtime
# libraries are not on the path and torch fails with "CUDA driver initialization failed"
# even though nvidia-smi sees the GPU. The PXR run script loaded this too.
module load cuda/11.8

BOLTZ="$HOME/micromamba/envs/boltz/bin/boltz"
[ -x "$BOLTZ" ] || { echo "FATAL: boltz not found at $BOLTZ"; exit 1; }
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
# `srun` is required, not optional. Without it the script body runs outside the job-step
# cgroup, so nvidia-smi reports both GPUs on the node while only one is allocated, and
# torch dies with "CUDA driver initialization failed" even though the same torch call
# succeeds under srun on the same node. PyTorch Lightning warns about this explicitly.
srun --gres=gpu:gtx1080ti:1 "$BOLTZ" predict submissions/01_boltz2_default/inputs \
  --out_dir submissions/01_boltz2_default/results \
  --output_format pdb \
  --use_msa_server
echo "EXIT=$?"
