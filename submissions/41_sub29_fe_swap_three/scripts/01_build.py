"""Submission 41: submission 29 (best Glide `inplace` gscore over three ESMFold2 runs,
the project best at 0.4638 LDDT-PLI) with three ligands swapped to a different run,
chosen because the Glide winner sits far from the haem iron.

The rule being tested: when the top-gscore pose leaves the ligand far from Fe, take the
next-best-scoring run instead, paying gscore for a shorter Fe-ligand distance.

    ligand          sub29 pick  Fe   gscore      new pick   Fe   gscore    rank of new
    OCNT-2315379    27 nomsa   6.84  -7.805  ->  26 maxacc  4.00  -6.490    3rd of 3
    OCNT-2395328    26 maxacc  4.49  -4.341  ->  27 nomsa   3.85  -4.225    2nd of 3
    OCNT-2312966    26 maxacc  4.39  -5.047  ->  23 default 4.18  -4.488    2nd of 3

Two caveats recorded before the result comes back.

First, OCNT-2315379 does not follow the stated rule. max_accuracy is the WORST of the
three runs on gscore for that ligand, not the second best; 23_default_msa is second at
-7.770 and sits at 6.89 A from iron, no closer than the pose being replaced. There is no
second-best option that also moves the ligand toward the metal, so this swap is a
straight Fe-distance swap, which is what submission 39 already ships on its own.

Second, OCNT-2312966 buys 0.21 A of Fe distance for 0.56 gscore units. That is inside
the run-to-run noise of these three models; the other two swaps are the ones with real
movement behind them.

Prior Fe-distance swaps in this project have all cost accuracy: submission 11, submission
22 at -0.0171 LDDT-PLI, submission 25 at -0.0062. Submission 39 (the OCNT-2315379 swap
alone) has not reported yet. If 39 and 41 both come back, the pair separates that one
swap from the two smaller ones.

COPY-ONLY: every structure is copied unchanged from a submission that already passed the
full gate 20/20, verified by sha256. PoseBusters is not re-run - a byte-identical copy
cannot return a different verdict.
"""
import csv
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from copy_check import copy_verified, report  # noqa: E402

BASE = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three", "submission")
SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
SWAPS = {"OCNT-2315379": "26_max_accuracy",
         "OCNT-2395328": "27_default_nomsa",
         "OCNT-2312966": "23_default_msa"}
SEL = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                   "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_sub29_fe_swap_three.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    sub29 = {r["ligand"]: r for r in csv.DictReader(open(SEL))}
    if len(sub29) != 20:
        raise SystemExit(f"expected 20 ligands in submission 29, found {len(sub29)}")
    missing = set(SWAPS) - set(sub29)
    if missing:
        raise SystemExit(f"swap targets absent from submission 29: {sorted(missing)}")

    rows, copies = [], []
    for lig in sorted(sub29):
        old = sub29[lig]["picked_set"]
        new = SWAPS.get(lig, old)
        if new == old:
            src = os.path.join(BASE, f"{lig}.pdb")
        else:
            src = os.path.join(ROOT, "submissions", SETS[new], "submission", f"{lig}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{lig}.pdb"),
                                    ligand=lig, source=new))
        g = {k[-2:]: float(sub29[lig][k])
             for k in ("gscore_23", "gscore_26", "gscore_27")}
        rows.append({"ligand": lig, "sub29_pick": old, "picked_set": new,
                     "changed": "yes" if new != old else "no",
                     "gscore_sub29": g[old.split("_")[0][-2:]],
                     "gscore_new": g[new.split("_")[0][-2:]],
                     "gscore_given_up": round(g[new.split("_")[0][-2:]]
                                              - g[old.split("_")[0][-2:]], 3)})

    chg = [r for r in rows if r["changed"] == "yes"]
    print(f"changed vs submission 29: {len(chg)}/20")
    for r in chg:
        print(f"  {r['ligand']}  {r['sub29_pick']} -> {r['picked_set']}  "
              f"gscore {r['gscore_sub29']:+.3f} -> {r['gscore_new']:+.3f} "
              f"({r['gscore_given_up']:+.3f})")
    if not report(copies):
        raise SystemExit("copy verification failed - zip not written")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(RES, "copy_manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(copies[0]))
        w.writeheader()
        w.writerows(copies)
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for r in rows:
            z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")
    print(f"wrote {os.path.relpath(ZIP, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
