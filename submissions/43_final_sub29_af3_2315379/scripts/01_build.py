"""Submission 43: the final entry - submission 29 with OCNT-2315379 taken from AF3.

Submission 29 is the project best at 0.4638 LDDT-PLI / 3.7277 BiSyRMSD / 0.8001 LDDT-LP:
per ligand, the best Glide `inplace` gscore over three ESMFold2 runs. Nineteen of its
twenty structures ship here unchanged.

The exception is OCNT-2315379, where submission 29's pick leaves the ligand 6.84 A from
the haem iron. All five candidates we have scored, `inplace`:

    source                Fe      gscore   evdw     ecoul    lipo    hbond   rewards
    27 nomsa  (sub29)     6.84   -7.805   -42.10   -1.625   -3.924    0.00   -2.135
    23 default + MSA      6.89   -7.770   -42.35   -0.375   -3.953    0.00   -2.246
    06 AF3, PL-iPTM       4.05   -7.469   -40.62   -3.956   -3.979   -0.32   -1.149
    18 Boltz-2, Glide     3.81   -6.610   -39.23   -3.374   -3.495   -0.32   -0.930
    26 max_accuracy       4.00   -6.490   -40.85   -3.573   -3.311   -0.32   -0.883

The AF3 structure is the compromise: it is the best-scoring pose that actually puts the
ligand near the metal, 0.336 gscore behind submission 29's pick and 0.979 ahead of the
ESMFold2 option at the same Fe distance (set 26, which submission 39 ships). The three
near-Fe poses all gain a hydrogen bond and roughly 2 kcal of electrostatics but give it
back on vdW and on Glide's reward terms - that, not a clash, is why they score worse.

What this is NOT supported by. Finding 4 in the README: short Fe distance is not a proxy
for correctness, and every Fe-distance swap this project has shipped lost ground
(submission 11; 22 at -0.0171; 25 at -0.0062). Finding 3: mixing models is the regime
where Glide selection fails, and AF3 06 has by far the worst pocket of the three models,
0.6985 LDDT-LP against 0.8001 here. One structure of twenty will not move LDDT-LP much,
but the direction is known and it is the wrong one.

What this IS supported by. 6.84 A puts the ligand outside contact with the metal
entirely, and OCNT-2315379 carries a nitrogen donor, so by finding 4's own reasoning
(analysis 19: among crystal structures only nitrogen-bearing ligands coordinate) it is
one of the ligands that plausibly should coordinate. Submission 34 showed LDDT-PLI
barely registers sub-angstrom moves, so a change this large is one of the few swaps that
could register at all.

The honest summary: the measurable risk is one structure's worth of a rule that has lost
three times out of three, against a pose that is physically implausible as it stands.
Submission 39 makes the same swap from ESMFold2 set 26 instead, and OCNT-2315379 sits
outside the scored subset (finding 5), so neither swap can be told apart from submission
29 on the public leaderboard - only on the final full-set evaluation.

COPY-ONLY: every structure is copied unchanged from a submission that already passed the
full gate 20/20, verified by sha256. PoseBusters is not re-run - a byte-identical copy
cannot return a different verdict.
"""
import collections
import csv
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from copy_check import copy_verified, report  # noqa: E402

SWAP = "OCNT-2315379"
DONOR = "06_af3_pliptm"
BASE = "29_glide_inplace_best_of_three"
SEL = os.path.join(ROOT, "submissions", BASE, "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_final_sub29_af3_2315379.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    picks = {r["ligand"]: r["picked_set"] for r in csv.DictReader(open(SEL))}
    if len(picks) != 20:
        raise SystemExit(f"expected 20 ligands, found {len(picks)}")
    if SWAP not in picks:
        raise SystemExit(f"{SWAP} is not in submission 29's selection")

    rows, copies = [], []
    for lig in sorted(picks):
        swapped = lig == SWAP
        src_dir = DONOR if swapped else BASE
        src = os.path.join(ROOT, "submissions", src_dir, "submission", f"{lig}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{lig}.pdb"),
                                    ligand=lig, source=src_dir))
        rows.append({"ligand": lig, "source": src_dir,
                     "sub29_pick": picks[lig],
                     "swapped": "yes" if swapped else "no"})

    if not report(copies):
        raise SystemExit("copy verification failed - zip not written")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for r in rows:
            z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")

    tally = collections.Counter(r["source"] for r in rows)
    print("sources: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"swapped {sum(1 for r in rows if r['swapped'] == 'yes')}/20")
    print(f"wrote {os.path.relpath(ZIP, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
