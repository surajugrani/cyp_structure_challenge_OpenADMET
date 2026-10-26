"""Submission 29: per-ligand best Glide `inplace` gscore across the three ESMFold2 runs.

The companion to submission 28, which used ICM RTCNN over the same three sets and scored
0.4411 / 3.7722 / 0.7965 - beating two of its three parents on LDDT-PLI and all three on
BiSyRMSD, though not submission 23.

    23  default endpoint + our ColabFold MSA     0.4527 / 3.8551 / 0.7985
    26  max_accuracy endpoint, its own MSA       0.4365 / 3.9024 / 0.8063
    27  default endpoint, no MSA                 0.4384 / 3.8017 / 0.7908

Why `inplace` rather than the `mininplace` used by every previous Glide selection here:

  * The two are NOT equivalent. Analysis 28 scored all 60 structures both ways: gscore
    correlates at r=0.791 between modes, and they pick a different structure for 7 of 20
    ligands. So the mode is a real choice, not a formality.
  * `mininplace` gives the WORSE score for 38 of 60 structures, median +0.177. The
    restrained minimisation appears to move the pose off the co-folded geometry, and
    POSTDOCK then rescores from there, so the displacement can cost more than the
    relaxation gains. On co-folded input - which is all we ship - that argues for
    scoring the pose as it stands.
  * The one argument for `mininplace` is that `inplace` returns emodel and einternal as
    0.00, confirmed here on all 60. That matters only if emodel is used. This selects on
    gscore alone, so nothing is lost.
  * Analysis 01, the only other head-to-head anywhere in this project, scored 20 crystal
    poses both ways against IC50 and found `inplace` the stronger (Spearman -0.265
    against -0.093), though neither was significant.

So this is the Glide half of the same experiment, under the mode that has never been
shipped. Against submission 28 it also asks a second question: the two scorers agree on
the winning set for only 10 of 20 ligands (Pearson between them +0.373), so if both
land near their parents' mean the agreement bound is uninformative, whereas a clear
split between 28 and 29 would say something about which scorer is tracking pose quality.

Expect it to differ from submission 23 on 14 of 20 ligands.

COPY-ONLY: every structure is copied unchanged from a submission that already passed the
full gate 20/20, so PoseBusters is not re-run - a byte-identical copy cannot return a
different verdict. Copies are verified by sha256.

Where a winning set does not SHIP its own ESMFold2 structure for a ligand - because it
failed PoseBusters there and was substituted - the shipped, already-validated file is
used, as in submission 28. The alternative would be shipping a structure known to fail.
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
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_glide_inplace_best_of_three.zip")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(SCORES)):
        if r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]:
            by[r["ligand"]][r["set"]] = float(r["gscore"])
    if len(by) != 20:
        raise SystemExit(f"expected 20 ligands scored, found {len(by)}")

    rows, copies = [], []
    for n in sorted(by):
        vals = by[n]
        if len(vals) != len(SETS):
            raise SystemExit(f"{n}: scored in {len(vals)} sets, expected {len(SETS)}")
        win = min(vals, key=vals.get)          # gscore: more negative is better
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        c = copy_verified(src, os.path.join(OUT, f"{n}.pdb"), ligand=n, source=win)
        copies.append(c)
        rows.append({"ligand": n, "picked_set": win,
                     "gscore_23": round(vals["23_default_msa"], 3),
                     "gscore_26": round(vals["26_max_accuracy"], 3),
                     "gscore_27": round(vals["27_default_nomsa"], 3),
                     "gscore_best": round(vals[win], 3),
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
