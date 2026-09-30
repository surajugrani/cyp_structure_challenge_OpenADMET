"""
Submission 01, step 4: run PoseBusters on the packaged structures.

PoseBusters is a hard gate - the organizers score a failing structure LDDT-PLI 0 - and it
is run separately from the submission validator replicated in 03_package.py, so passing
that does not imply passing this.

This mirrors the organizers' own invocation in evaluation/evaluate_predictions.py:

    pred_mol   = Chem.MolFromPDBFile(model_path, removeHs=False, sanitize=False)
    fragments  = Chem.SplitMolByPDBResidues(pred_mol); fragments.pop("LIG")
    pb_protein = everything else combined          -> mol_cond
    template   = Chem.MolFromSmiles(smiles)        -> mol_true
    PoseBusters(config="dock").bust(mol_pred=ligand, mol_true=template, mol_cond=pb_protein)

Note `fragments.pop("LIG")` removes only the ligand, so the haem stays in `mol_cond` and
ligand-haem clashes are checked - which matters for a CYP.

Their `evaluation/posebusters_utils.py` is vendored unmodified under scripts/organizer_code/
and used directly rather than reimplemented, because its
`assign_bond_orders_from_template` deliberately tries every topological match instead of
RDKit's first-match-wins, and its docstring documents three ways the naive call
misassigns bond orders on PDB-derived ligands.

The organizers allow up to POSEBUSTERS_MAX_FAILURES failed checks, a constant defined in
their private evaluation config, not in the public tutorial repo - so the count is
reported here rather than compared against a threshold. Treat any failure as a risk.

Usage:
    ~/micromamba/envs/cyp_posebusters_env/bin/python scripts/04_posebusters.py

Outputs:
    results/posebusters_report.csv   one row per structure, with any failed check names
"""

import csv
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SUBDIR = os.path.join(HERE, "submission")
RESULTS = os.path.join(HERE, "results")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def read_smiles():
    path = os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")
    with open(path) as f:
        return {r["Molecule_Name"].strip(): r["SMILES"].strip()
                for r in csv.DictReader(f) if r["Molecule_Name"].strip()}


def main():
    smiles = read_smiles()
    pb = PoseBusters(config="dock")
    rows = []

    for name in sorted(smiles):
        path = os.path.join(SUBDIR, f"{name}.pdb")
        rec = {"name": name}
        if not os.path.exists(path):
            rec.update(status="MISSING", n_failed=None, failed_checks="")
            rows.append(rec)
            continue

        pred = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
        frags = Chem.SplitMolByPDBResidues(pred)
        ligand = frags.pop("LIG", None)
        protein = None
        for frag in frags.values():
            protein = frag if protein is None else Chem.CombineMols(protein, frag)

        if ligand is None or protein is None:
            rec.update(status="PARSE_FAIL", n_failed=None, failed_checks="")
            rows.append(rec)
            print(f"  {name}: PARSE_FAIL")
            continue

        template = Chem.MolFromSmiles(smiles[name])
        try:
            ligand = assign_bond_orders_from_template(template, ligand)
        except Exception as e:
            rec["bond_order_note"] = f"{type(e).__name__}: {e}"

        res = pb.bust(mol_pred=ligand, mol_true=template, mol_cond=protein)
        bool_cols = [c for c in res.columns if res[c].dtype == bool]
        row = res[bool_cols].iloc[0]
        failed = sorted(row.index[~row].tolist())
        rec.update(status="pass" if not failed else "FAIL",
                   n_failed=len(failed), failed_checks="; ".join(failed))
        rows.append(rec)
        print(f"  {name}: {'pass' if not failed else f'{len(failed)} failed - ' + ', '.join(failed)}")

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "posebusters_report.csv")
    df.to_csv(out, index=False)
    n_pass = int((df.status == "pass").sum())
    print(f"\n{n_pass}/{len(df)} pass all PoseBusters checks")
    if n_pass != len(df):
        bad = df[df.status != "pass"]
        print("\nfailed checks by frequency:")
        cnt = {}
        for s in bad.failed_checks.dropna():
            for c in [x for x in s.split("; ") if x]:
                cnt[c] = cnt.get(c, 0) + 1
        for c, n in sorted(cnt.items(), key=lambda kv: -kv[1]):
            print(f"  {n:2d}x {c}")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
