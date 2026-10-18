"""How far did each ligand move between submissions 16 and 17, and who owns the gain?

16 (PL-iPTM selection)  0.4229 / 3.9283 / 0.8060
17 (complex_iPDE)       0.4157 / 3.8021 / 0.7976

BiSyRMSD improved 0.126 A - the best of any submission - while LDDT-PLI slipped 0.0072.
Only 8 of 20 structures differ. With ~10 complexes scored, a 0.126 A mean improvement is
~1.26 A summed, which one or two large moves could account for entirely.

This measures, for every ligand, the pocket-superposed symmetry-corrected RMSD between
16's pose and 17's, so the eight changes can be ranked by how much they actually moved.
Same machinery as analysis 13.
"""
import csv, math, os, sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A = os.path.join(ROOT, "submissions", "16_boltz2_pliptm", "submission")
B = os.path.join(ROOT, "submissions", "17_boltz2_ipde", "submission")
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
    picks = {r["ligand"]: r for r in csv.DictReader(open(os.path.join(
        ROOT, "submissions", "17_boltz2_ipde", "results", "ipde_picks.csv")))}

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
                     "differs": pk.get("differs_from_16") == "True",
                     "sub16_model": pk.get("pl_iptm_pick_sub16", ""),
                     "sub17_model": pk.get("ipde_pick", ""),
                     "rmsd_16_to_17": None if d is None else round(d, 3),
                     "fe_16": None if fe_a is None else round(
                         float(np.min(np.linalg.norm(la - np.array(fe_a), axis=1))), 2),
                     "fe_17": None if fe_b is None else round(
                         float(np.min(np.linalg.norm(xyz(lig_b) - np.array(fe_b),
                                                     axis=1))), 2)})

    with open(os.path.join(HERE, "results", "sub16_vs_17.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    diff = [r for r in rows if r["differs"]]
    same = [r for r in rows if not r["differs"]]
    print(f"{'ligand':16} {'16':>4} {'17':>4} {'RMSD 16->17':>12} {'Fe16':>6} {'Fe17':>6}")
    for r in sorted(diff, key=lambda x: -(x["rmsd_16_to_17"] or 0)):
        print(f"{r['ligand']:16} {r['sub16_model']:>4} {r['sub17_model']:>4} "
              f"{str(r['rmsd_16_to_17']):>12} {str(r['fe_16']):>6} {str(r['fe_17']):>6}")
    print(f"\n{len(diff)} structures differ, {len(same)} identical")
    v = [r["rmsd_16_to_17"] for r in diff if r["rmsd_16_to_17"] is not None]
    if v:
        import statistics as st
        print(f"movement among the {len(v)} that differ: median {st.median(v):.2f} A, "
              f"mean {st.mean(v):.2f} A, max {max(v):.2f} A")
        print(f"  >2 A: {sum(1 for x in v if x > 2)}   >4 A: {sum(1 for x in v if x > 4)}")
    z = [r["rmsd_16_to_17"] for r in same if r["rmsd_16_to_17"] is not None]
    if z:
        print(f"sanity - identical picks should be ~0: max {max(z):.3f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
