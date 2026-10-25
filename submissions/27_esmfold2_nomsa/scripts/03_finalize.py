"""Finalise submission 27: 19 no-MSA structures + OCNT-2313038 from submission 23.

Validation came back 20/20 structural, 19/20 PoseBusters. The single failure is
OCNT-2313038 on non-aromatic_ring_non-flatness - the SAME defect, on the same ligand,
that submission 23 hit. Submission 23 could not fix it either: PrepWizard re-minimised
at -rmsd 1.0 moved the ligand 0.19 A mean and left the ring flagged, so another job is
not worth it, and 23 took that structure from submission 19.

Taking it from submission 23 here is what keeps this a one-variable experiment. Both
submissions then carry the identical non-ESMFold2 structure at that position, so the
19 remaining structures - folded by the same endpoint, same esmfold2-2026-05 weights,
same sequence, differing only in whether our ColabFold MSA was supplied - are the whole
of the difference between 27 and 23. (That structure is transitively submission 19's,
in both.)

That the no-MSA run reproduces 23's one PoseBusters failure, on the same ligand and the
same check, is itself a result: the defect comes from the model, not from the MSA.

COPY-ONLY for the one substitution: a byte-identical copy of an already-validated
structure cannot return a different verdict, so it is checked by sha256. The 19 no-MSA
structures were validated in the job.
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
ZIP = os.path.join(HERE, "submission_esmfold2_nomsa.zip")

FROM_23 = {"OCNT-2313038": "non-aromatic_ring_non-flatness"}


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
                         "reason": f"no-MSA run failed {FROM_23[n]}"})
            print(f"  {n}: from submission 23 - no-MSA run failed {FROM_23[n]}")
        else:
            rows.append({"ligand": n, "source": "esmfold2_max_accuracy",
                         "reason": ""})

    print(f"\n{len(FROM_23)}/20 from submission 23, "
          f"{20 - len(FROM_23)}/20 ESMFold2 no-MSA")
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
