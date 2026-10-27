#!/bin/bash
#SBATCH --job-name=a31_depth
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH --output=analyses/31_msa_depth/slurm-%j.out
# RESOURCES: pure RDKit geometry on 20 small structures, single-threaded -> 1 CPU.
#            Analysis 27 ran the same ligand_rmsd comparison in under a minute at
#            well below 4G; 30 min is already generous. No licence, no Glide.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
PY=~/micromamba/envs/cyp_posebusters_env/bin/python
A=analyses/31_msa_depth

echo "=== mmCIF -> PDB for the depth-2048 run ==="
$PY $A/scripts/01_convert.py submissions/31_esmfold2_depth2048/raw $A/structures/depth2048

echo "=== ligand RMSD vs the three existing ESMFold2 sets ==="
for s in 23_esmfold2_default:d1024 26_esmfold2_maxacc:maxacc 27_esmfold2_nomsa:nomsa; do
  d=${s%%:*}; lab=${s##*:}
  $PY tools/ligand_rmsd.py $A/structures/depth2048 submissions/$d/submission \
      --label-a depth2048 --label-b $lab --out $A/results/vs_${lab}.csv
done
echo "EXIT=$?"
