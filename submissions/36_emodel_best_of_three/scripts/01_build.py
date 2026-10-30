"""Submission 36: best `emodel` over the same three ESMFold2 runs submission 29 used.

Same pool, same structures, different ranking term.

Submission 29 selected on `gscore` under `DOCKING_METHOD inplace` and scored 0.4638,
the best of the project. Two attempts to improve it by changing the POOL have failed in
opposite directions - submission 33 removed a candidate (-0.0174), submission 35 added
one (-0.0126) - which says the three-run pool is a local optimum. What has never been
changed is the quantity being ranked.

Why `emodel` is the obvious thing to try. `gscore` is Glide's estimate of binding
affinity: it answers "how well does this ligand bind this receptor". `emodel` is Glide's
pose-ranking term, built to choose among alternative poses OF THE SAME LIGAND, weighting
the force-field interaction energy far more heavily. Choosing among alternative
structures of the same ligand is exactly our task, and `emodel` is the term designed for
it. This project has used it for nothing.

Why it was never available until now. `inplace` returns `emodel` as 0.00; only
`mininplace` computes it. The project's original reason for preferring `mininplace` was
precisely to keep `emodel` - then every selection built on it ranked by `gscore` anyway,
so the term was paid for and never used. Analysis 28 scored all 60 structures in both
modes, so the numbers have been on disk since.

The objection this used to carry, and why it is now answered. `emodel` describes the
MINIMISED pose, while the structure shipped is the unminimised one, so the quantity
ranked is not a property of the file in the zip. That is the defect that made submission
22 uninterpretable. Submission 34 settled it: every ligand in submission 29 was moved to
its minimised position, a median 0.497 A displacement, and LDDT-PLI changed by 0.0001
with LDDT-LP identical to four decimals. The minimised and unminimised poses are the
same structure as far as the metric is concerned, so ranking one and shipping the other
is no longer a meaningful mismatch.

Measured before building: 8 of 20 structures change against submission 29, and the picks
come out 23x7 / 26x8 / 27x5 - more even than submission 29's 6/9/5 and far more even
than submission 35's skew toward its weakest member.

NO QUALITY CLAIM is made by any emodel value here. It is a scoring term; the leaderboard
is the only measurement of quality this project has.

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

SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
SCORES = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "results",
                      "scores_both_modes.csv")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_emodel_best_of_three.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    em = collections.defaultdict(dict)
    gs = collections.defaultdict(dict)
    for r in csv.DictReader(open(SCORES)):
        if r["status"] != "ok":
            continue
        if r["mode"] == "mininplace" and r.get("emodel"):
            em[r["ligand"]][r["set"]] = float(r["emodel"])
        if r["mode"] == "inplace" and r.get("gscore"):
            gs[r["ligand"]][r["set"]] = float(r["gscore"])

    short = {n: sorted(set(SETS) - set(v)) for n, v in em.items()
             if len(set(SETS) & set(v)) != len(SETS)}
    if len(em) != 20 or short:
        raise SystemExit(f"incomplete emodel scores: {len(em)} ligands, missing {short}")

    prev = {r["ligand"]: r["picked_set"] for r in csv.DictReader(open(SUB29))}

    rows, copies = [], []
    for n in sorted(em):
        vals = {s: em[n][s] for s in SETS}
        win = min(vals, key=vals.get)          # emodel: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{n}.pdb"),
                                    ligand=n, source=win))
        rows.append({"ligand": n, "picked_set": win,
                     **{f"emodel_{s}": round(vals[s], 2) for s in sorted(SETS)},
                     "emodel_best": round(vals[win], 2),
                     "gscore_of_pick": round(gs[n].get(win, float("nan")), 3),
                     "sub29_pick": prev.get(n, ""),
                     "changed_vs_sub29": "yes" if win != prev.get(n) else "no"})

    tally = collections.Counter(r["picked_set"] for r in rows)
    chg = [r for r in rows if r["changed_vs_sub29"] == "yes"]
    agree = 20 - len(chg)
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"changes vs submission 29: {len(chg)}/20  "
          f"(emodel and gscore agree on {agree})")
    for r in chg:
        print(f"  {r['ligand']}  {r['sub29_pick']} -> {r['picked_set']}")

    if not chg:
        raise SystemExit("picks identical to submission 29 - nothing to submit.")
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
