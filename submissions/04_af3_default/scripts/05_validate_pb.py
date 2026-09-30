"""Submission 04, step 4: run the organizers' checks + PoseBusters on the AF3 PDBs.

Same checks as submissions 01/03: .pdb files, exactly one residue named LIG, at most 3
chains, HEM present with 43 atoms, no duplicate residue numbers / atom names / atom
serials, ligand connectivity matching the given SMILES, then PoseBusters(config="dock")
with the haem left in mol_cond so ligand-haem clashes are checked.
"""
import csv, os, sys, zipfile
from collections import defaultdict
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SUB = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_af3_default.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")


def errs_for(lines, name, smi, tmp):
    e = []
    at = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if len({l[21] for l in at}) > 3:
        e.append("more than 3 chains")
    lig = [l for l in at if l[17:20].strip() == "LIG"]
    if not lig:
        e.append("no LIG residue")
    elif len({(l[21], l[22:27]) for l in lig}) > 1:
        e.append("more than one LIG residue")
    nhem = len([l for l in at if l[17:20].strip() == "HEM"])
    if nhem != 43:
        e.append(f"HEM has {nhem} atoms")
    seen = {}
    for l in at:
        k = (l[21], l[22:27])
        if k in seen and seen[k] != l[17:20].strip():
            e.append("duplicate residue number")
        seen[k] = l[17:20].strip()
    per = defaultdict(list)
    for l in at:
        per[(l[21], l[22:27])].append(l[12:16].strip())
    if any(len(v) != len(set(v)) for v in per.values()):
        e.append("duplicate atom names in a residue")
    s = [l[6:11] for l in at]
    if len(s) != len(set(s)):
        e.append("duplicate atom serials")
    p = os.path.join(tmp, f"{name}_lig.pdb")
    open(p, "w").writelines(lig + ["END\n"])
    ref = Chem.MolFromSmiles(smi)
    probe = Chem.MolFromPDBFile(p, removeHs=True, sanitize=False)
    if ref is None or probe is None:
        e.append("ligand parse failure")
    else:
        try:
            AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except ValueError:
            e.append("ligand connectivity != SMILES")
    return e


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    tmp = os.path.join(RES, "_tmp"); os.makedirs(tmp, exist_ok=True)
    pb = PoseBusters(config="dock")
    rows, nstruct, npb = [], 0, 0
    for n in sorted(smiles):
        f = os.path.join(SUB, f"{n}.pdb")
        if not os.path.exists(f):
            print(f"  {n}: MISSING"); rows.append({"name": n, "status": "MISSING"}); continue
        lines = open(f).readlines()
        se = errs_for(lines, n, smiles[n], tmp)
        pred = Chem.MolFromPDBFile(f, removeHs=False, sanitize=False)
        frags = Chem.SplitMolByPDBResidues(pred)
        lig = frags.pop("LIG", None)
        prot = None
        for fr in frags.values():
            prot = fr if prot is None else Chem.CombineMols(prot, fr)
        tpl = Chem.MolFromSmiles(smiles[n])
        try:
            lig = assign_bond_orders_from_template(tpl, lig)
        except Exception:
            pass
        r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot)
        bc = [c for c in r.columns if r[c].dtype == bool]
        failed = sorted(r[bc].iloc[0].index[~r[bc].iloc[0]].tolist())
        nstruct += not se; npb += not failed
        rows.append({"name": n, "struct_ok": not se, "struct_errors": "; ".join(se),
                     "pb_pass": not failed, "pb_failed": "; ".join(failed)})
        print(f"  {n}: struct {'ok' if not se else 'FAIL '+'; '.join(se)} | "
              f"PB {'pass' if not failed else 'FAIL '+', '.join(failed)}")
    with open(os.path.join(RES, "validation_report.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader(); w.writerows(rows)
    print(f"\nstructural {nstruct}/{len(smiles)}   PoseBusters {npb}/{len(smiles)}")
    if nstruct == len(smiles) and npb == len(smiles):
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for n in sorted(smiles):
                z.write(os.path.join(SUB, f"{n}.pdb"), f"{n}.pdb")
        print(f"wrote {ZIP}")
    else:
        print("zip not written - failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
