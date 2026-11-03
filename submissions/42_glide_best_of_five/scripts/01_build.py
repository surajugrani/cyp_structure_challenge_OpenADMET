"""Submission 42: best Glide `inplace` gscore over a FIVE-candidate pool.

The pool is the three ESMFold2 runs that submission 29 selected over (23 default +
ColabFold MSA, 26 max_accuracy, 27 no MSA) plus one structure each from the two other
co-folding models: AF3 picked by PL-iPTM (submission 06) and Boltz-2 picked by Glide
(submission 18). The AF3 and Boltz entries are each already the best of that model's
five samples, so this is a best-of-five-sources rule, not a best-of-all-samples one.

Submission 30 ran the cross-model rule with a three-candidate pool - AF3 06, Boltz 18
and ESM 23 only - and scored 0.4414, which is 0.0224 under submission 29. This widens
that pool to include the two ESMFold2 runs 30 left out, so an ESM structure now has
three chances to win a ligand rather than one.

Measured before building: 9 of 20 picks change against submission 29, 6 moving to
Boltz-2 and 3 to AF3.

The case against it is finding 3 in the README, restated by submissions 24 and 30:
Glide gscore judges how well a ligand fits the receptor it is handed and cannot judge
the receptor, so across models of unequal pocket quality it selects on the wrong thing.
Submissions 33 and 35 showed the three-run ESM pool is a local optimum - both removing
a run and adding a fourth cost accuracy. This adds two candidates from models whose
pockets are measurably worse (AF3 06 at 0.6985 LDDT-LP, Boltz 18 at 0.8047, against
0.8001 for submission 29), which is exactly the regime where that bias bites.

The case for it: submission 29 proved the `inplace` mode is what made Glide selection
work at all, and every cross-model failure (13, 18, 20, 24, 30 in part) predates or
partly predates that fix. Submission 30 is the one `inplace` cross-model test and it
used the narrowest possible pool.

COPY-ONLY: every structure is copied unchanged from a submission that already passed the
full gate 20/20, verified by sha256. PoseBusters is not re-run.
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

FOLDERS = {"23_default_msa": "23_esmfold2_default",
           "26_max_accuracy": "26_esmfold2_maxacc",
           "27_default_nomsa": "27_esmfold2_nomsa",
           "06_af3_pliptm": "06_af3_pliptm",
           "18_boltz2_glide": "18_boltz2_glide"}
EXTRA = os.path.join(ROOT, "analyses", "30_inplace_three_model", "results",
                     "inplace_scores.csv")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_glide_best_of_five.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    extra = collections.defaultdict(dict)
    for r in csv.DictReader(open(EXTRA)):
        if r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]:
            if r["set"] in ("06_af3_pliptm", "18_boltz2_glide"):
                extra[r["ligand"]][r["set"]] = float(r["gscore"])

    rows, copies = [], []
    for r in csv.DictReader(open(SUB29)):
        lig, old = r["ligand"], r["picked_set"]
        pool = {"23_default_msa": float(r["gscore_23"]),
                "26_max_accuracy": float(r["gscore_26"]),
                "27_default_nomsa": float(r["gscore_27"])}
        pool.update(extra.get(lig, {}))
        if len(pool) != 5:
            raise SystemExit(f"{lig}: pool has {len(pool)} candidates, expected 5")
        win = min(pool, key=pool.get)          # gscore: more negative is better
        src = os.path.join(ROOT, "submissions", FOLDERS[win], "submission", f"{lig}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{lig}.pdb"),
                                    ligand=lig, source=win))
        rows.append({"ligand": lig, "sub29_pick": old, "picked_set": win,
                     "changed": "yes" if win != old else "no",
                     "gscore_sub29": round(pool[old], 3),
                     "gscore_best_of_five": round(pool[win], 3),
                     "gscore_gained": round(pool[win] - pool[old], 3),
                     **{f"gscore_{k}": round(v, 3) for k, v in sorted(pool.items())}})
    if len(rows) != 20:
        raise SystemExit(f"expected 20 ligands, built {len(rows)}")

    tally = collections.Counter(x["picked_set"] for x in rows)
    chg = [x for x in rows if x["changed"] == "yes"]
    print("picked per source: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"changed vs submission 29: {len(chg)}/20")
    for x in chg:
        print(f"  {x['ligand']}  {x['sub29_pick']} -> {x['picked_set']}  "
              f"{x['gscore_sub29']:+.3f} -> {x['gscore_best_of_five']:+.3f}")
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
        for x in rows:
            z.write(os.path.join(OUT, f"{x['ligand']}.pdb"), f"{x['ligand']}.pdb")
    print(f"wrote {os.path.relpath(ZIP, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
