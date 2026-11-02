"""Submission 40: best Glide `mininplace` gscore over the three ESMFold2 runs.

The third and last of the Glide selectors analysis 28 made available on this pool:

    inplace gscore     shipped, 0.4638 - the project best
    emodel             shipped, 0.4493 - only available under mininplace
    mininplace gscore  this one, never built

`mininplace` restrained-minimises the pose before scoring it, so its gscore rates the
pose the ligand would relax into rather than the pose as co-folded. That is the one
difference from the 0.4638 submission; the structures shipped are the unminimised
originals either way, which submission 34 established is immaterial - minimising the
ligands moved LDDT-PLI by 0.0001.

Changes 7 of 20 picks against the `inplace` selection. All 60 scores already exist in
analysis 28, so this costs no compute.

Expectations are low. Both alternatives to `inplace` gscore have now lost: `emodel` by
0.0145, and every attempt to widen or narrow the candidate pool by more. The pool is a
local optimum and `inplace` gscore is the selector that found it.

This step only records the selection and stages the winning sets' structures. The poses
shipped are then replaced by 02_graft.py with the MINIMISED ones the selection was
actually computed from, so the object scored and the object submitted are the same thing.
Submission 22 shipped raw ICM poses while ranking them on minimised scores and its
-0.0171 therefore indicted nothing cleanly; this does not repeat that.
"""
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from copy_check import copy_verified, report  # noqa: E402

SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
SCORES = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "results",
                      "scores_both_modes.csv")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_mininplace_best_of_three.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(SCORES)):
        if r["status"] == "ok" and r["mode"] == "mininplace" and r["gscore"]:
            by[r["ligand"]][r["set"]] = float(r["gscore"])
    if len(by) != 20:
        raise SystemExit(f"expected 20 ligands scored, found {len(by)}")

    prev = {r["ligand"]: r["picked_set"] for r in csv.DictReader(open(SUB29))}

    rows, copies = [], []
    for lig in sorted(by):
        vals = by[lig]
        if len(vals) != len(SETS):
            raise SystemExit(f"{lig}: scored in {len(vals)} of the 3 sets")
        win = min(vals, key=vals.get)          # gscore: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{lig}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{lig}.pdb"),
                                    ligand=lig, source=win))
        rows.append({"ligand": lig, "picked_set": win,
                     "mininplace_23": round(vals["23_default_msa"], 3),
                     "mininplace_26": round(vals["26_max_accuracy"], 3),
                     "mininplace_27": round(vals["27_default_nomsa"], 3),
                     "mininplace_best": round(vals[win], 3),
                     "inplace_pick": prev[lig],
                     "changed_vs_inplace": "yes" if win != prev[lig] else "no"})

    if not report(copies):
        raise SystemExit("copy verification failed")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    tally = collections.Counter(r["picked_set"] for r in rows)
    chg = sum(1 for r in rows if r["changed_vs_inplace"] == "yes")
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"changed vs the inplace selection: {chg}/20")
    print("staged 20/20 - run 02_graft.py to replace the poses with the minimised ones")
    return 0


if __name__ == "__main__":
    sys.exit(main())
