"""Submission 35: best Glide `inplace` gscore over FOUR ESMFold2 runs.

Submission 29 ran this rule over three runs and scored 0.4638, the project best and the
only rule that has beaten plain ESMFold2. Submission 33 then removed a candidate from
that pool and scored 0.4464 - a loss of 0.0174, and below plain ESMFold2 as well. The
candidate it removed, set 26, was over-picked relative to chance and held the worst
interface score of the three, so on the face of it was the one to cut; cutting it cost.

That points at pool size rather than candidate quality as what the selector is living
on. This runs the converse to check it, and makes the test as sharp as the available
material allows by adding the WORST ESMFold2 run of the project as the fourth
candidate:

    23  default endpoint + ColabFold MSA      0.4527
    26  max_accuracy endpoint, own MSA        0.4365
    27  default endpoint, no MSA              0.4384
    31  default endpoint, msa_max_depth 2048  0.4116   <- added here

If pool size is the lever, a candidate 0.0411 below the best of the pool should still
help. If candidate quality is what matters, adding this one should hurt.

Set 31's structures are already a shipped, validated submission, so every copy here
comes from a set that passed the full 20/20 gate.

NOTE on OCNT-2313038: sets 23 and 31 both ship submission 19's structure for this
ligand, because ESMFold2 builds its ring non-flat and PoseBusters rejects it. For that
one ligand the fourth candidate is a duplicate of one already in the pool and cannot
change the pick. 19 of 20 ligands get a genuine fourth option.

WHAT THIS CANNOT SETTLE. The leaderboard scores 10 of the 20 structures we submit, and
we do not know which 10. If the picks that change fall in the unscored half, the result
comes back indistinguishable from submission 29 and means nothing either way. A null
result here is uninformative, not evidence for the hypothesis.

NO QUALITY CLAIM is made by any gscore below. Glide ranks candidates within one
receptor; that is the single use of it with a leaderboard result behind it. Whether a
better gscore means a better structure is what the leaderboard decides.

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
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa",
        "31_depth2048": "31_esmfold2_depth2048"}
SCORES3 = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "results",
                       "scores_both_modes.csv")
SCORES4 = os.path.join(ROOT, "analyses", "32_inplace_depth2048", "results",
                       "inplace_depth2048.csv")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_glide_best_of_four.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(SCORES3)):
        if r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]:
            by[r["ligand"]][r["set"]] = float(r["gscore"])
    for r in csv.DictReader(open(SCORES4)):
        if r["status"] == "ok" and r["gscore"]:
            by[r["ligand"]]["31_depth2048"] = float(r["gscore"])

    short = {n: sorted(set(SETS) - set(v)) for n, v in by.items()
             if len(set(SETS) & set(v)) != len(SETS)}
    if len(by) != 20 or short:
        raise SystemExit(f"incomplete scores: {len(by)} ligands, missing {short}")

    prev = {r["ligand"]: r["picked_set"] for r in csv.DictReader(open(SUB29))}

    rows, copies = [], []
    for n in sorted(by):
        vals = {s: by[n][s] for s in SETS}
        win = min(vals, key=vals.get)          # gscore: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{n}.pdb"),
                                    ligand=n, source=win))
        rows.append({"ligand": n, "picked_set": win,
                     **{f"gscore_{s}": round(vals[s], 3) for s in sorted(SETS)},
                     "gscore_best_of_four": round(vals[win], 3),
                     "sub29_pick": prev.get(n, ""),
                     "changed_vs_sub29": "yes" if win != prev.get(n) else "no",
                     "gscore_gained": round(by[n][prev[n]] - vals[win], 3)
                     if prev.get(n) else ""})

    tally = collections.Counter(r["picked_set"] for r in rows)
    chg = [r for r in rows if r["changed_vs_sub29"] == "yes"]
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"changes vs submission 29: {len(chg)}/20")
    print(f"of those, won by the new set 31: "
          f"{sum(1 for r in chg if r['picked_set'] == '31_depth2048')}")
    if chg:
        g = sorted(r["gscore_gained"] for r in chg)
        print(f"gscore gained on the changed ones: median {g[len(g)//2]:.2f}, "
              f"max {g[-1]:.2f}")
        for r in chg:
            print(f"  {r['ligand']}  {r['sub29_pick']} -> {r['picked_set']}"
                  f"  ({r['gscore_gained']:+.3f})")

    if not chg:
        raise SystemExit("picks identical to submission 29 - byte-identical zip, "
                         "nothing to submit. Not packaging.")
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
