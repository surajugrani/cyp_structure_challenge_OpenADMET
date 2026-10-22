"""Ligand RMSD between submission 19 (AF3+Boltz hybrid) and submission 24 (three-way).

Submission 24 is submission 19 with 5 structures replaced by ESMFold2, chosen because
Glide scored the ESMFold2 pose best of the three models. The Fe-distance figure shows
almost nothing: all 5 land within 0.2 A of the same distance from the haem iron, so the
coordination count is unchanged at 11/20.

That figure is the wrong instrument here. A crude order-independent check put the 5
swapped ligands 5.5-7.1 A apart by centroid, which says they sit at a similar radius
from the iron but elsewhere in the pocket. This measures it properly: pocket-superposed,
symmetry-corrected RMSD, the same machinery as analysis 13 and the way BiSyRMSD is
scored.

The 15 unchanged structures are copied byte-for-byte from submission 19 and must come
back at exactly 0.000 - they are the control for the measurement itself.

19  0.4331 / 3.6754 / 0.7823 - best LDDT-PLI before ESMFold2
24  pending
"""
import csv, math, os, sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A = os.path.join(ROOT, "submissions", "19_hybrid_glide", "submission")
B = os.path.join(ROOT, "submissions", "24_hybrid_three_way", "submission")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
RDLogger.DisableLog("rdApp.*")
POCKET_R = 10.0


def parse(p):
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


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    # which model submission 24 took, so the table says why a row moved
    picks = {r["ligand"]: r for r in csv.DictReader(open(os.path.join(
        ROOT, "submissions", "24_hybrid_three_way", "results", "selection.csv")))}

    rows = []
    for n in sorted(smiles):
        ca_a, lig_a, fe_a = parse(os.path.join(A, f"{n}.pdb"))
        ca_b, lig_b, fe_b = parse(os.path.join(B, f"{n}.pdb"))
        la = xyz(lig_a)
        keys = [k for k, v in ca_a.items()
                if np.min(np.linalg.norm(la - np.array(v), axis=1)) <= POCKET_R
                and k in ca_b]
        d = None
        if len(keys) >= 20:
            R, t = kabsch(np.array([ca_b[k] for k in keys]),
                          np.array([ca_a[k] for k in keys]))
            lb = (R @ xyz(lig_b).T).T + t
            ma, mb = ligmol(lig_a, smiles[n]), ligmol(lig_b, smiles[n])
            if ma is not None and mb is not None and \
                    ma.GetNumAtoms() == len(la) and mb.GetNumAtoms() == len(lb):
                try:
                    pr = Chem.Mol(mb); cf = pr.GetConformer()
                    for i, p in enumerate(lb):
                        cf.SetAtomPosition(i, p.tolist())
                    d = float(rdMolAlign.CalcRMS(pr, ma))
                except Exception:
                    pass
            if d is None and len(la) == len(lb):
                d = float(np.sqrt(((la - lb) ** 2).sum(1).mean()))
        pk = picks.get(n, {})
        rows.append({"ligand": n,
                     "differs": pk.get("changed_vs_19") == "yes",
                     "sub24_model": pk.get("source", ""),
                     "sub19_model": pk.get("sub19_took", ""),
                     "rmsd_19_to_24": None if d is None else round(d, 3),
                     "fe_19": None if fe_a is None else round(
                         float(np.min(np.linalg.norm(la - np.array(fe_a), axis=1))), 2),
                     "fe_24": None if fe_b is None else round(
                         float(np.min(np.linalg.norm(xyz(lig_b) - np.array(fe_b),
                                                     axis=1))), 2)})

    with open(os.path.join(HERE, "results", "sub19_vs_24.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    diff = rows
    same = []
    print(f"{'ligand':16} {'RMSD 19->24':>12} {'Fe 19':>6} {'Fe 24':>6}  source")
    for r in sorted(diff, key=lambda x: -(x["rmsd_19_to_24"] or 0)):
        print(f"{r['ligand']:16} {str(r['rmsd_19_to_24']):>12} "
              f"{str(r['fe_19']):>6} {str(r['fe_24']):>6}  {r['sub24_model']}"
              f"{'   <-- swapped' if r['differs'] else ''}")
    nz = [r for r in diff if (r["rmsd_19_to_24"] or 0) > 0.001]
    print(f"\n{len(diff)} structures compared, {len(nz)} actually moved")
    v = [r["rmsd_19_to_24"] for r in nz if r["rmsd_19_to_24"] is not None]
    if v:
        import statistics as st
        print(f"movement among the {len(v)} that differ: median {st.median(v):.2f} A, "
              f"mean {st.mean(v):.2f} A, max {max(v):.2f} A")
        print(f"  >2 A: {sum(1 for x in v if x > 2)}   >4 A: {sum(1 for x in v if x > 4)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
