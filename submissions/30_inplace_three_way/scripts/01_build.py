"""Three-model hybrid, selected by Glide `inplace` gscore: AF3 + Boltz-2 + ESMFold2.

This is the mininplace three-model hybrid (submission 24) redone under `inplace`.

Submission 24 took, per ligand, whichever of AF3 (submission 06), Boltz-2 (submission 18)
or ESMFold2 (submission 23) had the best `mininplace` gscore. It scored 0.4348 - which is
0.0179 BELOW simply shipping ESMFold2 alone - and was recorded at the time as killing the
cross-model selection rule.

That conclusion did not survive. The same selection over three ESMFold2 runs, scored
`inplace` instead, became the best submission of the project:

    ESMFold2 default (our MSA)          0.4527 / 3.8551 / 0.7985
    best inplace gscore of 3 ESM runs   0.4638 / 3.7277 / 0.8001   <- beat its best
                                                                      parent on all three
    best RTCNN of the same 3 runs       0.4411 / 3.7722 / 0.7965

Every Glide selection that ever failed here used `mininplace`; the one that used
`inplace` won. Analysis 28 supplies the mechanism: `mininplace` gives the WORSE score for
38 of 60 structures (median +0.177), because its restrained minimisation moves the pose
off the co-folded geometry before POSTDOCK rescores it. For structures we never docked,
scoring the pose as it stands is the right thing to do. The usual argument for
`mininplace` - that `inplace` zeroes emodel and einternal - does not apply, because this
ranks on gscore alone.

So this asks the obvious question: does the three-model hybrid also flip? The ESMFold2
result is one prospective test, and one result is not a rule. These three parents are
further apart than three runs of one model, which is the regime where submission 19
(AF3 + Boltz by mininplace) appeared to work before submission 24 undercut it.

Scores: AF3 and Boltz-2 from analysis 30, ESMFold2 from analysis 28, all `inplace`, each
computed on that structure's own receptor and grid.

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

# set key -> (submission folder, scores csv)
A30 = os.path.join(ROOT, "analyses", "30_inplace_three_model", "results",
                   "inplace_scores.csv")
A28 = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "results",
                   "scores_both_modes.csv")
SETS = {"06_af3_pliptm": "06_af3_pliptm",
        "18_boltz2_glide": "18_boltz2_glide",
        "23_default_msa": "23_esmfold2_default"}
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_inplace_three_way.zip")
SUB24 = os.path.join(ROOT, "submissions", "24_hybrid_three_way", "results",
                     "selection.csv")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(A30)):
        if r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]:
            by[r["ligand"]][r["set"]] = float(r["gscore"])
    for r in csv.DictReader(open(A28)):
        if (r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]
                and r["set"] == "23_default_msa"):
            by[r["ligand"]][r["set"]] = float(r["gscore"])

    if len(by) != 20:
        raise SystemExit(f"expected 20 ligands, found {len(by)}")

    # what submission 24 chose, to report how many picks the mode changes
    prev = {}
    if os.path.exists(SUB24):
        m = {"af3_06": "06_af3_pliptm", "boltz_18": "18_boltz2_glide",
             "esmfold2_23": "23_default_msa"}
        prev = {r["ligand"]: m.get(r["source"], r["source"])
                for r in csv.DictReader(open(SUB24))}

    rows, copies, changed = [], [], 0
    for n in sorted(by):
        vals = by[n]
        win = min(vals, key=vals.get)          # gscore: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        c = copy_verified(src, os.path.join(OUT, f"{n}.pdb"), ligand=n, source=win)
        copies.append(c)
        same_as_24 = prev.get(n) == win
        changed += (prev.get(n) is not None and not same_as_24)
        rows.append({"ligand": n, "picked_set": win,
                     "n_scored": len(vals),
                     "gscore_af3_06": vals.get("06_af3_pliptm", ""),
                     "gscore_boltz_18": vals.get("18_boltz2_glide", ""),
                     "gscore_esm_23": vals.get("23_default_msa", ""),
                     "gscore_best": round(vals[win], 3),
                     "sub24_mininplace_pick": prev.get(n, ""),
                     "changed_vs_sub24": "" if prev.get(n) is None
                                         else ("no" if same_as_24 else "yes")})

    tally = collections.Counter(r["picked_set"] for r in rows)
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    nfull = sum(1 for r in rows if r["n_scored"] == 3)
    print(f"{nfull}/20 ligands had all three sets scored")
    if prev:
        print(f"differs from the mininplace hybrid (submission 24) on {changed}/20")
    print(f"differs from ESMFold2 alone on "
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
