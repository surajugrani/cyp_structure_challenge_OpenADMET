#!/bin/bash
#SBATCH --job-name=cyp_icmout
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-icmout-%j.out
# Does _dockScan output=<.sdf> give scores correctly bound to each pose?
# Needed because parsing SCORES> from stdout mis-assigns scores to poses.
set -euo pipefail
module load icm/3.9.4
cd "$SLURM_SUBMIT_DIR"/structures/4D6Z
rm -f hits.sdf D_4D6Z_4D6Z_ligand_v20001.ob
icm64 /opt/applications/icm/3.9.4/_dockScan D_4D6Z input=4D6Z_ligand_v2000.sdf \
      confs=10 effort=3. seed=42 output=hits.sdf -s -a -S > dockscan.log 2>&1 || true
echo "=== hits.sdf present:"; ls -l hits.sdf 2>/dev/null || echo MISSING
echo "=== molecules and SD tags:"
grep -c '^\$\$\$\$' hits.sdf 2>/dev/null || echo 0
grep -oE '^> *<[^>]+>' hits.sdf 2>/dev/null | sort -u | head -20
