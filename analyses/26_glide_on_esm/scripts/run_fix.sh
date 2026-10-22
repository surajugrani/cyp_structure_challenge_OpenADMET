#!/bin/bash
#SBATCH --job-name=a26_fix
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=analyses/26_glide_on_esm/slurm-fix-%j.out
# RESOURCES: 1 structure, serial -> 1 CPU; grid+mininplace+PoseBusters ~2-3G -> 4G;
#            ~5 min -> 30m.
set -uo pipefail
module load schrodinger/2025-3
export PYTHONPATH=/gpfs/home/sugrani/micromamba/envs/cyp_posebusters_env/lib/python3.11/site-packages:${PYTHONPATH:-}
cd "$SLURM_SUBMIT_DIR"
$SCHRODINGER/run python3 analyses/26_glide_on_esm/scripts/02_fix_2313038.py
echo "EXIT=$?"
