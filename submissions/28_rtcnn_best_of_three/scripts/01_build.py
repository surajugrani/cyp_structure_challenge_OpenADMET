"""Submission 28: per-ligand best ICM RTCNN score across the three ESMFold2 runs.

RTCNN is the only scoring function besides Glide that is comparable across programs: a
co-folding model's confidence is calibrated to itself, RTCNN is not. It had only ever
been applied within AF3's ensemble (submission 14, 0.3949 against PL-iPTM's 0.4291).
Analysis 29 scored all 60 ESMFold2 structures so it could be used across runs instead:

    23  default endpoint + our ColabFold MSA     0.4527   best of the project
    26  max_accuracy endpoint, its own MSA       0.4365
    27  default endpoint, no MSA                 0.4384

THE SELECTION LOOKS LIKE CHANCE. RTCNN picks 23 six times, 26 six times and 27 eight
times, against 6.67 expected from a coin. Submission 23 outscores the other two by about
0.015 LDDT-PLI, so a selector with any cross-run skill should prefer it well above
chance. It does not. This is the same diagnostic that condemned submission 24, where
Glide preferred ESMFold2 for 5 of 20 despite ESMFold2 beating both alternatives
outright, and that submission duly landed below its best parent.

The structures are also close together - analysis 27 measured a median pocket-superposed
ligand RMSD of 1.53 A between sets 23 and 27, the within-ensemble regime where every
selector tried here has failed. Built because it was asked for and the scoring was
already done; the expectation is a result near the mean of the three parents.

TWO LIGANDS WHERE THE PICK CANNOT BE HONOURED. RTCNN scored each set's own PrepWizard
output, but two of the winning sets do not SHIP their own structure for that ligand,
because it failed PoseBusters there:

    OCNT-2313038   RTCNN picks 23, which ships submission 19's structure - the ESMFold2
                   ring fails flatness in all three runs, so no set ships its own
    OCNT-2395789   RTCNN picks 26, which ships submission 23's structure after
                   max-accuracy failed bond_angles

In both cases the shipped, already-validated file is used. The alternative would be to
ship a structure that fails PoseBusters, which is not an option. So the selection is
effectively overridden on 2 of 20, and both resolve toward submission 23.

COPY-ONLY: every structure is copied unchanged from a submission that already passed the
full gate 20/20, so PoseBusters is not re-run - a byte-identical copy cannot return a
different verdict. Copies are verified by sha256.
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

SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
SCORES = os.path.join(ROOT, "analyses", "29_icm_on_esm", "results",
                      "icm_esm_scores.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_rtcnn_best_of_three.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(SCORES)):
        if r["status"] == "ok" and r["RTCNN"]:
            by[r["ligand"]][r["set"]] = float(r["RTCNN"])
    if len(by) != 20:
        raise SystemExit(f"expected 20 ligands scored, found {len(by)}")

    rows, copies = [], []
    for n in sorted(by):
        vals = by[n]
        if len(vals) != len(SETS):
            raise SystemExit(f"{n}: scored in {len(vals)} sets, expected {len(SETS)}")
        win = min(vals, key=vals.get)          # RTCNN: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        c = copy_verified(src, os.path.join(OUT, f"{n}.pdb"), ligand=n, source=win)
        copies.append(c)
        rows.append({"ligand": n, "picked_set": win,
                     "rtcnn_23": round(vals["23_default_msa"], 3),
                     "rtcnn_26": round(vals["26_max_accuracy"], 3),
                     "rtcnn_27": round(vals["27_default_nomsa"], 3),
                     "rtcnn_best": round(vals[win], 3),
                     "spread": round(max(vals.values()) - min(vals.values()), 3)})

    tally = collections.Counter(r["picked_set"] for r in rows)
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"differs from submission 23 on "
          f"{sum(1 for r in rows if r['picked_set'] != '23_default_msa')}/20")
    if not report(copies):
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
