"""Submission 39: the best-of-three co-folded submission with OCNT-2315379 swapped in
from the max_accuracy ESMFold2 run.

Submission 29 - best Glide `inplace` gscore over three ESMFold2 runs - scored 0.4638,
the project best. Its pick for OCNT-2315379 came from the no-MSA run and leaves the
ligand 6.84 A from the haem iron, the largest Fe-ligand separation in the whole
submission. Across the three runs that ligand sits at 6.89 A (default + our ColabFold
MSA), 4.00 A (max_accuracy) and 6.84 A (no MSA, the pick). The max_accuracy structure is
the only candidate anywhere near the iron, and Glide ranked it LAST by 1.315 gscore units.

This swaps that one structure. Everything else is submission 29 unchanged.

THE TRACK RECORD ON THIS MOVE IS BAD, and it is on record before the result comes back:

    11  AF3 pliptm with 3 Fe-outliers from Boltz-2      0.4241
    22  hybrid with 5 Glide-ranked ICM poses            0.4160   -0.0171
    25  default ESM run with 3 structures from 19       0.4465   -0.0062

Three attempts, three losses. Fe-ligand distance says a structure is unusual, not that it
is wrong, and the leaderboard has never once rewarded acting on it. There is no evidence
here that 4.00 A is closer to the truth than 6.84 A - only that it is closer to the iron.

What makes this one marginally different from those three: it swaps a single structure
rather than three or five, and the replacement comes from INSIDE submission 29's own
candidate pool rather than from another model, so it does not reintroduce the
cross-model mixing that submissions 30 and 35 both paid for. If it loses, it loses about
a twentieth of the blend, and the Fe-distance heuristic is finished for good.

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
DONOR = "26_esmfold2_maxacc"
BASE = "29_glide_inplace_best_of_three"
SEL = os.path.join(ROOT, "submissions", BASE, "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_sub29_fe_swap_2315379.zip")


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
