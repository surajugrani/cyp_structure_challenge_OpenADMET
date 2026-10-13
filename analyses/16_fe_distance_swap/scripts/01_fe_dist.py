"""Fe-ligand distance for submission 06 (AF3) vs submission 03 (Boltz-2), per ligand.

Hypothesis under test: AF3 misplaces the ligand entirely for a few structures - wrong
cavity or a peripheral site rather than the catalytic one - and Boltz-2 gets the pocket
right for those. If so, the Fe-ligand distance exposes it without needing any ground
truth, because CYP3A4 turns substrates over at the haem iron.

Threshold care: 12 of the 20 blinded ligands are type II binders that coordinate Fe
through an imidazole / pyridine / triazole N at 2.06-2.32 A. The other 8 are type I and
sit in the catalytic pocket WITHOUT coordinating, typically 3.5-5 A out, which is correct
for them. So the diagnostic is not "fails to coordinate" - it is "too far to be in the
catalytic site at all". Reported here as the raw distribution so the cut can be chosen
from the data rather than assumed.

Measures min distance from the haem Fe to any ligand heavy atom, in the exact files each
submission shipped.
"""
import csv, os, sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S06 = os.path.join(ROOT, "submissions", "06_af3_pliptm", "submission")
S03 = os.path.join(ROOT, "submissions", "03_boltz2_prepwizard", "submission")


def fe_dist(path):
    fe, lig = None, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if el == "FE":
            fe = xyz
        elif l[17:20].strip() == "LIG":
            lig.append(xyz)
    if fe is None or not lig:
        return None, len(lig)
    return float(np.min(np.linalg.norm(np.array(lig) - np.array(fe), axis=1))), len(lig)


def main():
    names = sorted(r["Molecule_Name"].strip() for r in csv.DictReader(
        open(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
        if r["Molecule_Name"].strip())
    rows = []
    for n in names:
        a, na = fe_dist(os.path.join(S06, f"{n}.pdb")) if os.path.exists(
            os.path.join(S06, f"{n}.pdb")) else (None, 0)
        b, nb = fe_dist(os.path.join(S03, f"{n}.pdb")) if os.path.exists(
            os.path.join(S03, f"{n}.pdb")) else (None, 0)
        rows.append({"ligand": n,
                     "af3_06_fe": None if a is None else round(a, 2),
                     "boltz_03_fe": None if b is None else round(b, 2),
                     "boltz_closer_by": None if (a is None or b is None)
                     else round(a - b, 2),
                     "n_lig_atoms": na})
    with open(os.path.join(HERE, "results", "fe_distance.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print(f"{'ligand':16} {'AF3 (06)':>9} {'Boltz (03)':>11} {'diff':>7}")
    for r in sorted(rows, key=lambda x: -(x["af3_06_fe"] or 0)):
        print(f"{r['ligand']:16} {str(r['af3_06_fe']):>9} {str(r['boltz_03_fe']):>11} "
              f"{str(r['boltz_closer_by']):>7}")
    a = np.array([r["af3_06_fe"] for r in rows if r["af3_06_fe"] is not None])
    b = np.array([r["boltz_03_fe"] for r in rows if r["boltz_03_fe"] is not None])
    print(f"\nAF3  (06): median {np.median(a):.2f}  n<=2.5A {int((a<=2.5).sum())}/{len(a)}  "
          f"n>6A {int((a>6).sum())}  max {a.max():.2f}")
    print(f"Boltz(03): median {np.median(b):.2f}  n<=2.5A {int((b<=2.5).sum())}/{len(b)}  "
          f"n>6A {int((b>6).sum())}  max {b.max():.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
