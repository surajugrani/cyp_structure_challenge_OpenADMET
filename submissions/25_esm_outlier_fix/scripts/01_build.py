"""Submission 25: submission 23 with its three worst-justified structures replaced.

Submission 23 (pure ESMFold2) is the best result of the project at 0.4527 LDDT-PLI, and
everything tried since has made it worse - submission 24 let Glide choose across three
co-folding models and landed 0.0179 BELOW shipping 23 unaltered. So this changes as
little as possible: 3 of 20 structures, each taken from submission 19.

    OCNT-2315379   Fe 6.89 -> 4.05 A    has an N donor yet sits 2.4 A further from the
                                        haem than anything else in the set
    OCNT-2395460   Fe 1.79 -> 3.71 A    NO N donor yet 1.79 A from the iron - too close
                                        to be a real contact, and analysis 19 found 0 of
                                        7 such crystal ligands coordinate at all
    OCNT-2395328   Fe 4.03 -> 3.01 A    Glide -3.520 -> -6.234, the largest per-ligand
                                        gap between any two models in the whole set

The first two are corrections of impossible values and move in OPPOSITE directions, which
is what separates them from the Fe-distance heuristic that lost in submissions 11 and 22:
that rule moved ligands toward the iron on the assumption that nearer is better, whereas
these are corrected toward chemical plausibility.

OCNT-2395328 was initially held back as a variant, on the argument that it has no
nitrogen donor so ESMFold2's 4.03 A is the more defensible distance and moving it to
3.01 A repeats the failed heuristic. That argument looked only at Fe distance and was
incomplete: Glide scores submission 19's structure at -6.234 against ESMFold2's -3.520,
a 2.7 kcal/mol gap, by far the largest disagreement between models on any ligand here.
Two independent signals pointing opposite ways is a reason to prefer the one with more
behind it, and the Glide gap is the stronger of the two - Fe distance has never once
predicted a leaderboard gain, while the 3-4 A regime is where cross-model Glide
differences are at least real. Shipped at the user's direction.

COPY-ONLY: every structure is copied unchanged from a submission that already passed the
full gate 20/20 (chain layout, ligand connectivity against SMILES, PoseBusters `dock`).
PoseBusters scores each complex in isolation, so a byte-identical copy cannot return a
different verdict; re-running it would prove nothing. Copies are verified by sha256
instead, which is the only failure mode a copy actually has. No SLURM job needed.
"""
import csv
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from copy_check import copy_verified, report  # noqa: E402

S23 = os.path.join(ROOT, "submissions", "23_esmfold2_default", "submission")
S19 = os.path.join(ROOT, "submissions", "19_hybrid_glide", "submission")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_esm_outlier_fix.zip")

SWAP = {
    "OCNT-2315379": "Fe 6.89 -> 4.05 A; N donor present but furthest from haem in the set",
    "OCNT-2395460": "Fe 1.79 -> 3.71 A; no N donor, 1.79 A is too close to be real",
    "OCNT-2395328": "Glide -3.520 -> -6.234, largest cross-model gap in the set",
}


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    rows = []
    for n in sorted(smiles):
        swapped = n in SWAP
        src = os.path.join(S19 if swapped else S23, f"{n}.pdb")
        r = copy_verified(src, os.path.join(OUT, f"{n}.pdb"), ligand=n,
                          source="sub_19" if swapped else "sub_23_esmfold2")
        r["reason"] = SWAP.get(n, "")
        rows.append(r)
        if swapped:
            print(f"  {n}: from submission 19 - {SWAP[n]}")

    nsw = sum(1 for r in rows if r["reason"])
    print(f"\n{nsw}/20 from submission 19, {20 - nsw}/20 unchanged from submission 23")
    if not report(rows):
        raise SystemExit("copy verification failed - zip not written")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for r in rows:
            z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")
    print(f"wrote {os.path.relpath(ZIP, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
