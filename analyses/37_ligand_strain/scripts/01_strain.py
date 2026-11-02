"""Ligand internal strain for the three ESMFold2 runs, from PoseBusters.

The stored validation CSVs record only which PoseBusters checks failed, not the numbers
behind them. The `dock` configuration does run an Energy ratio module - the pose's
conformer energy against an ensemble of generated conformers - but `bust` returns only
the boolean verdict unless `full_report=True` is passed. A first attempt without it
produced 60 rows with no numeric columns at all. This re-runs with the full report and
keeps every non-boolean column.

The question it serves: does the best submission so far contain ligands that are
internally strained yet score well with Glide? Glide gscore rewards fit to the receptor
and does not price the ligand's own conformational cost, so a strained pose can win a
selection it should lose.

pb_fail's molecule preparation is imported from submission 23's script rather than
copied, so the two cannot diverge in how the ligand is split out and bond orders are
assigned.
"""
import argparse
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "submissions", "23_esmfold2_default", "scripts",
                   "02_convert_validate.py")
OUT_CSV = os.path.join(HERE, "results", "strain.csv")
SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}

spec = importlib.util.spec_from_file_location("conv23", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def numbers(path, smi, pb):
    """Everything pb.bust returns that is not a boolean check."""
    from rdkit import Chem
    pred = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    if pred is None:
        return {"status": "unparseable"}
    frags = Chem.SplitMolByPDBResidues(pred)
    lig = frags.pop("LIG", None)
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    if lig is None or prot is None:
        return {"status": "no LIG or protein"}
    tpl = Chem.MolFromSmiles(smi)
    try:
        lig = m.assign_bond_orders_from_template(tpl, lig)
    except Exception:
        pass
    r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot, full_report=True)
    row = r.iloc[0]
    out = {"status": "ok"}
    for c in r.columns:
        if r[c].dtype != bool:
            try:
                out[c] = float(row[c])
            except (TypeError, ValueError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.parse_args()
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = m.PoseBusters(config="dock")

    rows, cols = [], ["ligand", "set", "status"]
    for setname, subdir in SETS.items():
        for lig in sorted(smiles):
            p = os.path.join(ROOT, "submissions", subdir, "submission", f"{lig}.pdb")
            if not os.path.exists(p):
                print(f"{setname} {lig}: MISSING", flush=True)
                continue
            d = numbers(p, smiles[lig], pb)
            d.update(ligand=lig, set=setname)
            for k in d:
                if k not in cols:
                    cols.append(k)
            rows.append(d)
            if not rows[:-1]:
                print("numeric columns: " + ", ".join(
                    k for k in d if k not in ("status", "ligand", "set")), flush=True)
            print(f"{setname} {lig}: ok", flush=True)

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {os.path.relpath(OUT_CSV, ROOT)}  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
