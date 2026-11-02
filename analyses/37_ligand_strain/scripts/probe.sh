#!/bin/bash
#SBATCH --job-name=pb_probe
#SBATCH -p shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -t 00:10:00
#SBATCH --output=analyses/37_ligand_strain/slurm-probe-%j.out
# RESOURCES: imports posebusters and prints its config file contents, no structures
#            loaded. Trivial; 1 CPU, 2G and 10 min is the smallest ask available.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
~/micromamba/envs/cyp_posebusters_env/bin/python -c "
import posebusters, os, yaml
d = os.path.join(os.path.dirname(posebusters.__file__), 'config')
print('configs:', sorted(os.listdir(d)))
for f in sorted(os.listdir(d)):
    if f.endswith(('.yml','.yaml')):
        c = yaml.safe_load(open(os.path.join(d,f)))
        mods = c.get('modules', [])
        names = [m.get('name', m.get('function','?')) if isinstance(m,dict) else m for m in mods]
        print(f, '->', names)
"
echo "EXIT=$?"
