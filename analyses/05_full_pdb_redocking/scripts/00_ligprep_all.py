"""
Analysis 05 step 0: run LigPrep for every ligand in the set, up front.

Ligand preparation is a shared input to BOTH the Glide pipeline (analysis 05) and the ICM
pipeline (analysis 06, which reuses these exact states so ligand prep is held constant).
Doing it as an explicit first pass decouples the two docking runs: they can then proceed
independently instead of ICM having to trail Glide structure by structure.

Epik ionization/tautomers at pH 7.4. `-g` respects the stereochemistry given in the SMILES
rather than enumerating it - every ligand in this set has stereo-defined SMILES from the
RCSB chemical component dictionary.

Serial - one Schrodinger licence token.

Usage (via scripts/run_ligprep.sh):
    $SCHRODINGER/run python3 scripts/00_ligprep_all.py
"""

import json
import os
import subprocess

import pandas as pd
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")
LIGPREP = os.path.join(os.environ["SCHRODINGER"], "ligprep")


def require_compute_node(what):
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(f"REFUSING to run {what} outside SLURM. Use: sbatch scripts/run_ligprep.sh")


def main():
    s = pd.read_csv(os.path.join(DATA, "structure_set.csv"))
    smiles = json.load(open(os.path.join(DATA, "ligand_smiles.json")))
    os.makedirs(RESULTS, exist_ok=True)

    rows = []
    for i, row in enumerate(s.itertuples(), 1):
        pdb, comp = row.pdb_id, row.comp_id
        wdir = os.path.join(STRUCT_DIR, pdb)
        os.makedirs(wdir, exist_ok=True)
        out = os.path.join(wdir, f"{pdb}_ligprep.maegz")
        rec = {"pdb_id": pdb, "comp_id": comp}
        try:
            if os.path.exists(out):
                rec["n_states"] = sum(1 for _ in structure.StructureReader(out))
                rec["status"] = "cached"
            else:
                smi = (smiles.get(comp) or {}).get("smiles_stereo")
                if not smi:
                    raise RuntimeError(f"no SMILES for {comp}")
                smi_file = os.path.join(wdir, f"{pdb}_lig.smi")
                with open(smi_file, "w") as f:
                    f.write(f"{smi} {pdb}_{comp}\n")
                require_compute_node("LigPrep")
                subprocess.run([LIGPREP, "-NOJOBID", "-WAIT",
                                "-ismi", os.path.basename(smi_file),
                                "-omae", os.path.basename(out),
                                "-epikx", "-ph", "7.4", "-pht", "2.0", "-g", "-bff", "14"],
                               cwd=wdir, check=False)
                if not os.path.exists(out):
                    raise RuntimeError("LigPrep produced no output")
                rec["n_states"] = sum(1 for _ in structure.StructureReader(out))
                rec["status"] = "ok"
            print(f"[{i}/{len(s)}] {pdb} ({comp}): {rec['n_states']} state(s) "
                  f"[{rec['status']}]", flush=True)
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"[{i}/{len(s)}] {pdb} ({comp}): ERROR {e}", flush=True)
        rows.append(rec)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, "ligprep_status.csv"), index=False)
    good = df[df["status"].isin(["ok", "cached"])]
    print(f"\n{len(good)}/{len(df)} ligands prepared "
          f"({(df.status == 'cached').sum()} cached). "
          f"States per ligand: min {good.n_states.min():.0f}, "
          f"median {good.n_states.median():.0f}, max {good.n_states.max():.0f}")


if __name__ == "__main__":
    main()
