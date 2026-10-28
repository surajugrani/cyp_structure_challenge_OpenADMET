"""Submission 33: best Glide `inplace` gscore over TWO ESMFold2 runs, 23 and 27.

Submission 29 ran the same rule over three runs and scored 0.4638, the project best.
This drops one of the three and asks whether the selector was better off without it.

Which one, and why. Submission 29's 20 winning picks broke down as:

    23  default endpoint + ColabFold MSA      6 picks   0.4527 PLI / 0.7985 LP
    26  max_accuracy endpoint, own MSA        9 picks   0.4365 PLI / 0.8063 LP
    27  default endpoint, no MSA              5 picks   0.4384 PLI / 0.7908 LP

Submission 26 supplies 9 of 20 winners against 6.7 expected by chance, despite having
the WORST protein-ligand interface score of the three. It also has the BEST pocket
geometry of any submission in this project, LDDT-LP 0.8063.

That pairing is the suspicious part, and there is a mechanism for it. Glide gscore
measures how well a ligand fits the receptor it is handed; it has no way to judge the
receptor. A better-built pocket therefore earns a better gscore whether or not the
ligand is correctly placed. Submission 30 paid 0.0113 for exactly this blind spot by
selecting across models of unequal pocket quality. Within this three-run pool,
submission 26 is the set with the conspicuously better pocket, so it is the one whose
picks are most likely to reflect pocket quality rather than pose accuracy.

Dropping it leaves 23 and 27, which are the cleanest pair available: same
`fold_all_atom` endpoint, same esmfold2-2026-05 weights, same sequence and chain
layout, differing only in whether our ColabFold MSA was supplied. Comparable candidate
quality is the precondition that submissions 30 and 31 both established.

Measured before building: this changes 9 of 20 structures against submission 29, giving
up a median 0.57 and at most 1.08 gscore units on the ones it changes.

The honest risk: submission 29 WON with those nine picks from set 26, so whatever bias
exists was not fatal. Removing a third of the candidate pool from a rule that is
already working can simply cost accuracy. A result below 0.4638 means the pool mattered
more than the bias; at or above it means Glide was being fooled by pocket quality.

COPY-ONLY: every structure is copied unchanged from a submission that already passed
the full gate 20/20, verified by sha256. PoseBusters is not re-run - a byte-identical
copy cannot return a different verdict.
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
        "27_default_nomsa": "27_esmfold2_nomsa"}
DROPPED = "26_max_accuracy"
SCORES = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "results",
                      "scores_both_modes.csv")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_glide_best_of_two.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(SCORES)):
        if r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]:
            by[r["ligand"]][r["set"]] = float(r["gscore"])
    if len(by) != 20:
        raise SystemExit(f"expected 20 ligands scored, found {len(by)}")

    prev = {r["ligand"]: r["picked_set"] for r in csv.DictReader(open(SUB29))}

    rows, copies = [], []
    for n in sorted(by):
        vals = {s: v for s, v in by[n].items() if s in SETS}
        if len(vals) != len(SETS):
            raise SystemExit(f"{n}: scored in {len(vals)} of the kept sets")
        win = min(vals, key=vals.get)          # gscore: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{n}.pdb"),
                                    ligand=n, source=win))
        rows.append({"ligand": n, "picked_set": win,
                     "gscore_23": round(by[n]["23_default_msa"], 3),
                     "gscore_27": round(by[n]["27_default_nomsa"], 3),
                     "gscore_26_dropped": round(by[n][DROPPED], 3),
                     "gscore_best_of_two": round(vals[win], 3),
                     "sub29_pick": prev.get(n, ""),
                     "changed_vs_sub29": "yes" if win != prev.get(n) else "no",
                     "gscore_given_up": round(vals[win] - by[n][prev[n]], 3)
                     if prev.get(n) else ""})

    tally = collections.Counter(r["picked_set"] for r in rows)
    chg = [r for r in rows if r["changed_vs_sub29"] == "yes"]
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"changes vs submission 29: {len(chg)}/20")
    print(f"differs from submission 23 alone: "
          f"{sum(1 for r in rows if r['picked_set'] != '23_default_msa')}/20")
    if chg:
        gu = sorted(r["gscore_given_up"] for r in chg)
        print(f"gscore given up on the changed ones: median {gu[len(gu)//2]:.2f}, "
              f"max {gu[-1]:.2f}")
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
