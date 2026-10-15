#!/bin/bash
#SBATCH --job-name=boltz_cyp15
#SBATCH -p gpu
#SBATCH --gres=gpu:gtx1080ti:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH --output=submissions/15_boltz2_tuned/slurm-boltz-%j.out
# Boltz-2 with PXR's tuned sampling, all 20 blinded ligands, each with the haem cofactor.
#
# Submission 03 (and 01) ran pure defaults: recycling 3, sampling 200, 1 diffusion sample.
# On the PXR challenge those knobs were worth a real amount:
#
#   PXR 8   ColabFold MSA, defaults          0.4439
#   PXR 9   recycling 5, sampling 300        0.4621   <- best, +0.0182
#   PXR 10  recycling 6, sampling 400        0.4548   more is NOT better
#   PXR 11  recycling 5, sampling 300 + 3 templates   0.4516
#
# So 5/300 is the sweet spot; 6/400 gave back a third of the gain. +0.018 is worth
# chasing here - it is more than double the 0.0070 that separated every selection-metric
# variant we tried (submissions 05, 06, 10).
#
# --diffusion_samples 5 is added beyond what PXR varied. Analysis 21 found Boltz's
# advantage over AF3 is that it ADAPTS the pocket per ligand (AF3 matched crystal 6DAA
# for 15 of 20 ligands; Boltz spread over five), so a diffusion ensemble should give
# genuine conformational spread. We cannot yet pick the best member - analyses 09, 10, 20
# and submissions 13/14 all failed at that - so this first run simply takes Boltz's OWN
# top-ranked sample and compares it against default Boltz (submission 03). The ensemble
# is kept on disk for when a working selector exists.
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
# Job 48107149 died here with "CUDA driver initialization failed" on nodeb0419, the same
# node where submission 01 ran fine with this identical recipe - so it is a node/driver
# state issue, not the script. These two lines make the next failure diagnosable.
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}  node=$(hostname)"
srun --gres=gpu:gtx1080ti:1 python -c "import torch; print('torch sees', torch.cuda.device_count(), 'GPU(s); cuda ok:', torch.cuda.is_available())" 2>&1 | tail -2
# `srun` is required, not optional. Without it the script body runs outside the job-step
# cgroup, so nvidia-smi reports both GPUs on the node while only one is allocated, and
# torch dies with "CUDA driver initialization failed" even though the same torch call
# succeeds under srun on the same node. PyTorch Lightning warns about this explicitly.
srun --gres=gpu:gtx1080ti:1 "$BOLTZ" predict submissions/15_boltz2_tuned/inputs \
  --out_dir submissions/15_boltz2_tuned/results \
  --output_format pdb \
  --use_msa_server \
  --recycling_steps 5 \
  --sampling_steps 300 \
  --diffusion_samples 5
echo "EXIT=$?"
