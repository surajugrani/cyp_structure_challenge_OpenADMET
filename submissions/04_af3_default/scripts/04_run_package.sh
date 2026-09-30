#!/bin/bash
#SBATCH --job-name=af3_pack
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --output=submissions/04_af3_default/slurm-pack-%j.out
# Select best AF3 sample by protein-ligand iPTM, convert mmCIF -> PDB with the ligand
# renamed to LIG, then run the organizers' validation checks and PoseBusters.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module load schrodinger/2025-3
echo "########## selection + conversion ##########"
$SCHRODINGER/run python3 submissions/04_af3_default/scripts/03_select_and_package.py
echo "########## validation + PoseBusters ##########"
~/micromamba/envs/cyp_posebusters_env/bin/python \
  submissions/04_af3_default/scripts/05_validate_pb.py
echo "EXIT=$?"
