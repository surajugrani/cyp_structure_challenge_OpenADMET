"""Did PrepWizard fix the PoseBusters failures in the AF3 samples?

The Glide rescoring pipeline PrepWizard'd all 100 AF3 samples as a side effect. Those
prepped complexes were never PoseBusters-checked - the pass/fail numbers in analysis 10
are for the RAW AF3 poses. PrepWizard took Boltz from 17/20 to 20/20 in submission 03, so
the question is whether it does the same here.

Compares pass/fail per sample, before (raw AF3) and after (PrepWizard), and reports the
ligand-Fe distance change, since several AF3 failures were Fe contacts that were too short
and PrepWizard could push the ligand off the iron.
"""
import csv, glob, os, sys
import numpy as np
from rdkit import Chem, RDLogger

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
WORK = os.path.join(HERE, "work")
RES = os.path.join(HERE, "results")
RAW = os.path.join(ROOT, "analyses", "10_af3_selection_vs_posebusters",
                   "results", "sample_level.csv")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")


def fe_dist(mol):
    c = mol.GetConformer(); L = []; FE = None
    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        ri = a.GetPDBResidueInfo(); x = np.array(c.GetAtomPosition(a.GetIdx()))
        if ri.GetResidueName().strip() == "LIG":
            L.append(x)
        elif ri.GetName().strip() == "FE":
            FE = x
    if FE is None or not L:
        return None
    return float(np.linalg.norm(np.array(L) - FE, axis=1).min())


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    raw = {(r["ligand"], r["sample"]): r["pb_pass"] == "True"
           for r in csv.DictReader(open(RAW))}
    pb = PoseBusters(config="dock")

    rows = []
    files = sorted(glob.glob(os.path.join(WORK, "*", "*_pw.pdb")))
    print(f"{len(files)} PrepWizard'd samples found")
    for i, f in enumerate(files, 1):
        base = os.path.basename(f)[:-7]          # strip _pw.pdb
        lig_name, tag = base.split("__")
        pred = Chem.MolFromPDBFile(f, removeHs=False, sanitize=False)
        if pred is None:
            rows.append({"ligand": lig_name, "sample": tag, "after_pass": None,
                         "note": "unparseable"})
            continue
        frags = Chem.SplitMolByPDBResidues(pred)
        lig = frags.pop("LIG", None)
        prot = None
        for fr in frags.values():
            prot = fr if prot is None else Chem.CombineMols(prot, fr)
        if lig is None or prot is None:
            rows.append({"ligand": lig_name, "sample": tag, "after_pass": None,
                         "note": "no LIG after prepwizard"})
            continue
        tpl = Chem.MolFromSmiles(smiles[lig_name])
        try:
            lig = assign_bond_orders_from_template(tpl, lig)
        except Exception:
            pass
        r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot)
        bc = [c for c in r.columns if r[c].dtype == bool]
        failed = sorted(r[bc].iloc[0].index[~r[bc].iloc[0]].tolist())
        rows.append({"ligand": lig_name, "sample": tag,
                     "before_pass": raw.get((lig_name, tag)),
                     "after_pass": not failed,
                     "after_failed": "; ".join(failed),
                     "fe_dist_after": fe_dist(pred), "note": ""})
        if i % 20 == 0:
            print(f"  {i}/{len(files)}", flush=True)

    with open(os.path.join(RES, "pb_after_prepwizard.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader(); w.writerows(rows)

    done = [r for r in rows if r.get("after_pass") is not None
            and r.get("before_pass") is not None]
    b = sum(r["before_pass"] for r in done)
    a = sum(r["after_pass"] for r in done)
    fixed = sum(1 for r in done if not r["before_pass"] and r["after_pass"])
    broke = sum(1 for r in done if r["before_pass"] and not r["after_pass"])
    print(f"\n=== sample level (n={len(done)}) ===")
    print(f"  before PrepWizard: {b}/{len(done)} pass")
    print(f"  after  PrepWizard: {a}/{len(done)} pass   (fixed {fixed}, broke {broke})")

    byl = {}
    for r in done:
        byl.setdefault(r["ligand"], []).append(r)
    print(f"\n=== ligand level (n={len(byl)}) ===")
    print(f"  ligands with >=1 passing sample BEFORE: "
          f"{sum(any(x['before_pass'] for x in v) for v in byl.values())}/{len(byl)}")
    print(f"  ligands with >=1 passing sample AFTER : "
          f"{sum(any(x['after_pass'] for x in v) for v in byl.values())}/{len(byl)}")
    fe = [r["fe_dist_after"] for r in done if r["fe_dist_after"]]
    print(f"\nligand-Fe distance after PrepWizard: min {min(fe):.2f} "
          f"median {np.median(fe):.2f} max {max(fe):.2f}")


if __name__ == "__main__":
    main()
