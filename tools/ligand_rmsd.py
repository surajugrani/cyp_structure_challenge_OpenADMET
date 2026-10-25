"""Pocket-superposed, symmetry-corrected ligand RMSD between two submissions.

This is the same quantity the organizers score as BiSyRMSD: superpose on the binding-site
CA atoms, then compare ligand poses with symmetry taken into account, so a flipped
phenyl does not read as a 2 A error.

It had been copied into analysis 22 three times (19 vs 20, 06 vs 18, 19 vs 24) with only
the two directory names changed. One tool instead:

    python3 tools/ligand_rmsd.py A_DIR B_DIR [--out results.csv] [--label-a X]

Each DIR holds one <ligand>.pdb per ligand, chain-A protein / HEM / LIG, as every
submission/ folder here does. Needs RDKit, so run it on a compute node.

Identical files give exactly 0.000, which is the control: a submission assembled by
copying should return 0 for every structure it did not change.
"""
import argparse
import csv
import os
import sys

import numpy as np
from rdkit import Chem, RDLogger

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from rdkit.Chem import rdMolAlign  # noqa: E402
RDLogger.DisableLog("rdApp.*")
POCKET_R = 10.0


def parse(p):
    """CA atoms by residue, ligand lines, and the iron position."""
    ca, lig, fe = {}, [], None
    for l in open(p):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if el == "FE":
            fe = xyz
        if l[17:20].strip() == "LIG":
            lig.append(l)
        elif l[12:16].strip() == "CA":
            ca[(l[21], l[22:27])] = xyz
    return ca, lig, fe


def xyz(ls):
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in ls])


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def ligmol(ls, smi):
    m = Chem.MolFromPDBBlock("".join(ls) + "END\n", removeHs=True, sanitize=False)
    if m is None:
        return None
    try:
        return assign_bond_orders_from_template(Chem.MolFromSmiles(smi), m)
    except Exception:
        return m


def pair_rmsd(pa, pb, smi):
    """Superpose B's pocket onto A's, then symmetry-corrected ligand RMSD."""
    ca_a, lig_a, fe_a = parse(pa)
    ca_b, lig_b, fe_b = parse(pb)
    la = xyz(lig_a)
    keys = [k for k, v in ca_a.items()
            if np.min(np.linalg.norm(la - np.array(v), axis=1)) <= POCKET_R
            and k in ca_b]
    if len(keys) < 20:
        return None, None, None
    R, t = kabsch(np.array([ca_b[k] for k in keys]),
                  np.array([ca_a[k] for k in keys]))
    lb = (R @ xyz(lig_b).T).T + t
    ma, mb = ligmol(lig_a, smi), ligmol(lig_b, smi)
    d = None
    if ma is not None and mb is not None and \
            ma.GetNumAtoms() == len(la) and mb.GetNumAtoms() == len(lb):
        try:
            pr = Chem.Mol(mb)
            cf = pr.GetConformer()
            for i, p in enumerate(lb):
                cf.SetAtomPosition(i, p.tolist())
            d = float(rdMolAlign.CalcRMS(pr, ma))
        except Exception:
            pass
    if d is None and len(la) == len(lb):
        d = float(np.sqrt(((la - lb) ** 2).sum(1).mean()))
    fa = None if fe_a is None else float(
        np.min(np.linalg.norm(la - np.array(fe_a), axis=1)))
    fb = None if fe_b is None else float(
        np.min(np.linalg.norm(xyz(lig_b) - np.array(fe_b), axis=1)))
    return d, fa, fb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a_dir")
    ap.add_argument("b_dir")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--out")
    a = ap.parse_args()
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    rows = []
    for n in sorted(smiles):
        pa, pb = os.path.join(a.a_dir, f"{n}.pdb"), os.path.join(a.b_dir, f"{n}.pdb")
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        d, fa, fb = pair_rmsd(pa, pb, smiles[n])
        rows.append({"ligand": n,
                     "rmsd": None if d is None else round(d, 3),
                     f"fe_{a.label_a}": None if fa is None else round(fa, 2),
                     f"fe_{a.label_b}": None if fb is None else round(fb, 2)})

    print(f"{'ligand':16} {'RMSD (A)':>9} {'Fe ' + a.label_a:>8} {'Fe ' + a.label_b:>8}")
    for r in sorted(rows, key=lambda x: -(x["rmsd"] or 0)):
        print(f"{r['ligand']:16} {str(r['rmsd']):>9} "
              f"{str(r[f'fe_{a.label_a}']):>8} {str(r[f'fe_{a.label_b}']):>8}")
    v = [r["rmsd"] for r in rows if r["rmsd"] is not None]
    moved = [x for x in v if x > 0.001]
    if v:
        import statistics as st
        print(f"\n{len(v)} compared, {len(moved)} differ at all")
        if moved:
            print(f"among those that differ: median {st.median(moved):.2f} A, "
                  f"mean {st.mean(moved):.2f} A, max {max(moved):.2f} A")
            print(f"  >1 A: {sum(1 for x in moved if x > 1)}   "
                  f">2 A: {sum(1 for x in moved if x > 2)}   "
                  f">4 A: {sum(1 for x in moved if x > 4)}")
    if a.out and rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {os.path.relpath(a.out, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
