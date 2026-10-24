"""Finalise submission 26: 16 max-accuracy structures + 4 from submission 23.

Pass 3 of the pipeline validated the max-accuracy set at 20/20 structural but 16/20
PoseBusters. The four failures:

    OCNT-2313038   non-aromatic_ring_non-flatness
    OCNT-2395436   bond_angles
    OCNT-2395464   bond_angles
    OCNT-2395789   bond_angles

Three of these passed cleanly in submission 23 as pure ESMFold2 output, so the defect is
specific to the max-accuracy endpoint rather than to the ligand: it fixes two chemically
impossible Fe distances (OCNT-2315379 6.89 -> 4.02 A, OCNT-2395460 1.79 -> 3.82 A) while
degrading local ligand geometry on three others. OCNT-2313038 fails the same way it did
in submission 23, where PrepWizard at -rmsd 1.0 moved the ligand only 0.19 A mean and
left the ring flagged, so re-minimising is not worth another job.

All four are therefore taken from submission 23. That keeps the submission entirely
ESMFold2 - two endpoints of the same model rather than a cross-model hybrid - and every
replacement is a structure that already passed the full gate 20/20. (OCNT-2313038 is
transitively submission 19's, as it was in 23.)

COPY-ONLY for the four: byte-identical copies of already-validated structures cannot
return a different PoseBusters verdict, so they are checked by sha256 instead. The 16
max-accuracy structures were validated in pass 3 and are taken as they stand.
"""
import csv
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from copy_check import copy_verified, report  # noqa: E402

OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
S23 = os.path.join(ROOT, "submissions", "23_esmfold2_default", "submission")
ZIP = os.path.join(HERE, "submission_esmfold2_maxacc.zip")

FROM_23 = {"OCNT-2313038": "non-aromatic_ring_non-flatness",
           "OCNT-2395436": "bond_angles",
           "OCNT-2395464": "bond_angles",
           "OCNT-2395789": "bond_angles"}


def main():
    val = {r["ligand"]: r for r in
           csv.DictReader(open(os.path.join(RES, "validation.csv")))}
    if len(val) != 20:
        raise SystemExit(f"expected 20 validated rows, found {len(val)}")

    # Confirm the four we are replacing are exactly the ones that failed, and that
    # nothing else did - otherwise the hardcoded list is stale.
    failed = {n for n, r in val.items() if r["pb_failed"] or r["struct_errors"]}
    if failed != set(FROM_23):
        raise SystemExit(f"failure set changed: {sorted(failed)} "
                         f"vs expected {sorted(FROM_23)}")

    rows, copies = [], []
    for n in sorted(val):
        if n in FROM_23:
            c = copy_verified(os.path.join(S23, f"{n}.pdb"),
                              os.path.join(OUT, f"{n}.pdb"),
                              ligand=n, source="sub_23_esmfold2_default")
            copies.append(c)
            rows.append({"ligand": n, "source": "sub_23_esmfold2_default",
                         "reason": f"max-accuracy failed {FROM_23[n]}"})
            print(f"  {n}: from submission 23 - max-accuracy failed {FROM_23[n]}")
        else:
            rows.append({"ligand": n, "source": "esmfold2_max_accuracy",
                         "reason": ""})

    print(f"\n{len(FROM_23)}/20 from submission 23, "
          f"{20 - len(FROM_23)}/20 ESMFold2 max-accuracy")
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
