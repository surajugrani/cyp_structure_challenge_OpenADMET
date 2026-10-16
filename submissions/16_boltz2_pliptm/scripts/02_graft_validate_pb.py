"""
Submission 03, step 2: graft PrepWizard coordinates onto the Boltz PDBs, re-run the
organizers' submission checks, and re-run PoseBusters.

Grafting instead of submitting PrepWizard's own PDB: PrepWizard adds hydrogens and
rewrites the file, and the validator is sensitive to atom naming, duplicate atom names
within a residue and duplicate atom serials. The Boltz files already pass all of that, so
only the x/y/z columns of matching heavy atoms are replaced, keyed on
(chain, residue number, atom name).

Usage:
    ~/micromamba/envs/cyp_posebusters_env/bin/python scripts/02_graft_validate_pb.py

Outputs:
    submission/<name>.pdb
    submission_boltz2_pliptm.zip     (only if all 20 pass both check sets)
    results/prepwizard_report.csv
"""

import csv
import os
import sys
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "submissions", "16_boltz2_pliptm", "raw_sel")
WORK = os.path.join(HERE, "work")
OUT = os.path.join(HERE, "submission")
RESULTS = os.path.join(HERE, "results")
ZIPPATH = os.path.join(HERE, "submission_boltz2_pliptm.zip")

sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def key(l):
    return (l[21], l[22:27].strip(), l[12:16].strip())


def structural_errors(lines, name, smiles, tmpdir):
    errs = []
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if len({l[21] for l in atoms}) > 3:
        errs.append("more than 3 chains")
    lig = [l for l in atoms if l[17:20].strip() == "LIG"]
    if not lig:
        errs.append("no LIG residue")
    elif len({(l[21], l[22:27]) for l in lig}) > 1:
        errs.append("more than one LIG residue")
    if len([l for l in atoms if l[17:20].strip() == "HEM"]) != 43:
        errs.append("HEM atom count != 43")
    seen = {}
    for l in atoms:
        k = (l[21], l[22:27])
        if k in seen and seen[k] != l[17:20].strip():
            errs.append("duplicate residue number")
        seen[k] = l[17:20].strip()
    per = defaultdict(list)
    for l in atoms:
        per[(l[21], l[22:27])].append(l[12:16].strip())
    if any(len(v) != len(set(v)) for v in per.values()):
        errs.append("duplicate atom names in a residue")
    ser = [l[6:11] for l in atoms]
    if len(ser) != len(set(ser)):
        errs.append("duplicate atom serials")
    p = os.path.join(tmpdir, f"{name}_lig.pdb")
    with open(p, "w") as f:
        f.writelines(lig + ["END\n"])
    ref = Chem.MolFromSmiles(smiles)
    probe = Chem.MolFromPDBFile(p, removeHs=True, sanitize=False)
    if ref is None or probe is None:
        errs.append("ligand parse failure")
    else:
        try:
            AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except ValueError:
            errs.append("ligand connectivity != SMILES")
    return errs


def geometry(path):
    """min ligand-protein distance and min ligand-Fe distance, heavy atoms."""
    m = Chem.MolFromPDBFile(path, removeHs=True, sanitize=False)
    c = m.GetConformer()
    L, P, FE = [], [], None
    for a in m.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        ri = a.GetPDBResidueInfo()
        xyz = np.array(c.GetAtomPosition(a.GetIdx()))
        if ri.GetResidueName().strip() == "LIG":
            L.append(xyz)
        else:
            P.append(xyz)
            if ri.GetName().strip() == "FE":
                FE = xyz
    L, P = np.array(L), np.array(P)
    return (float(np.linalg.norm(L[:, None] - P[None], axis=2).min()),
            float(np.linalg.norm(L - FE, axis=1).min()) if FE is not None else np.nan)


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    tmpdir = os.path.join(RESULTS, "_tmp")
    os.makedirs(tmpdir, exist_ok=True)
    pb = PoseBusters(config="dock")

    rows = []
    for name in sorted(smiles):
        src = os.path.join(SRC, f"{name}.pdb")
        pw = os.path.join(WORK, f"{name}_pw.pdb")
        rec = {"name": name}
        if not (os.path.exists(src) and os.path.exists(pw)):
            rec.update(status="MISSING")
            rows.append(rec)
            continue

        new = {}
        for l in open(pw):
            if l.startswith(("ATOM", "HETATM")):
                el = l[76:78].strip() or l[12:16].strip()[0]
                if el != "H":
                    new[key(l)] = (l[30:38], l[38:46], l[46:54])

        out_lines, disp, unmatched = [], [], 0
        for l in open(src):
            if l.startswith(("ATOM", "HETATM")):
                k = key(l)
                if k in new:
                    old = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                    x, y, z = new[k]
                    disp.append(np.linalg.norm(
                        old - np.array([float(x), float(y), float(z)])))
                    l = l[:30] + x + y + z + l[54:]
                else:
                    unmatched += 1
            out_lines.append(l)
        dst = os.path.join(OUT, f"{name}.pdb")
        with open(dst, "w") as f:
            f.writelines(out_lines)

        errs = structural_errors(out_lines, name, smiles[name], tmpdir)
        before = geometry(src)
        after = geometry(dst)

        pred = Chem.MolFromPDBFile(dst, removeHs=False, sanitize=False)
        frags = Chem.SplitMolByPDBResidues(pred)
        lig = frags.pop("LIG", None)
        prot = None
        for fr in frags.values():
            prot = fr if prot is None else Chem.CombineMols(prot, fr)
        tpl = Chem.MolFromSmiles(smiles[name])
        try:
            lig = assign_bond_orders_from_template(tpl, lig)
        except Exception:
            pass
        res = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot)
        bc = [c for c in res.columns if res[c].dtype == bool]
        row = res[bc].iloc[0]
        pb_failed = sorted(row.index[~row].tolist())

        rec.update(unmatched=unmatched,
                   med_disp=float(np.median(disp)), max_disp=float(np.max(disp)),
                   min_lig_prot_before=before[0], min_lig_prot_after=after[0],
                   lig_fe_before=before[1], lig_fe_after=after[1],
                   struct_ok=not errs, struct_errors="; ".join(errs),
                   pb_pass=not pb_failed, pb_failed="; ".join(pb_failed))
        rows.append(rec)
        flag = "ok" if (not errs and not pb_failed) else \
               ("PB:" + ",".join(pb_failed) if pb_failed else "STRUCT:" + "; ".join(errs))
        print(f"  {name}: disp {rec['med_disp']:.2f}/{rec['max_disp']:.2f} A | "
              f"clash {before[0]:.2f}->{after[0]:.2f} | Fe {before[1]:.2f}->{after[1]:.2f} | {flag}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, "prepwizard_report.csv"), index=False,
              float_format="%.3f")
    ns = int(df.struct_ok.sum())
    npb = int(df.pb_pass.sum())
    print(f"\nstructural checks: {ns}/{len(df)}   PoseBusters: {npb}/{len(df)}")
    if ns == len(df) and npb == len(df):
        with zipfile.ZipFile(ZIPPATH, "w", zipfile.ZIP_DEFLATED) as z:
            for name in sorted(smiles):
                z.write(os.path.join(OUT, f"{name}.pdb"), f"{name}.pdb")
        print(f"wrote {ZIPPATH}")
    else:
        print("zip not written - see failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
