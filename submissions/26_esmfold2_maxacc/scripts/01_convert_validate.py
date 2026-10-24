"""Submission 26: ESMFold2's fold_max_accuracy endpoint, all 20 complexes.

Submission 23 used `fold_all_atom` with our 14,301-sequence ColabFold MSA and scored
0.4527 LDDT-PLI, the best of the project. This uses a different endpoint on the same
model family, and the differences are forced rather than chosen:

  * `fold_max_accuracy` exposes no quality knobs - FoldMaxAccuracyConfig carries only
    the three include_* output flags; the server picks the settings.
  * It REFUSES a supplied MSA and builds its own, so this run is MSA-free from our side.
  * It pins its own weights, ESMFOLD2_CUTOFF_2025, not the esmfold2-2026-05 of 23.
  * It is rate-limited to 10 complexes per day, so the set took two days to fold.

Why it is worth a submission: it independently corrects BOTH structures that submission
25 had to borrow from submission 19 to fix.

    OCNT-2315379   6.89 -> 4.02 A   the far outlier in 23, with an N donor
    OCNT-2395460   1.79 -> 3.82 A   the impossible 1.79 A contact on a ligand with no
                                    nitrogen donor

It reached those from its own MSA with no knowledge of that reasoning, which is a
genuine independent check on it. Everything else barely moves: 14 of 20 agree with
submission 23 within 0.1 A, median Fe distance 2.20 against 2.24 A, coordination
unchanged at 12 of the 15 ligands that can coordinate. The one structure that moves the
other way is OCNT-2395328, 4.03 -> 4.48 A.

No structure is taken from another submission. Submission 23 had to borrow
OCNT-2313038 because its ESMFold2 ring failed PoseBusters and PrepWizard could not flatten
it; whether the max-accuracy pose has the same defect is unknown until this runs, so
SWAP_FROM_19 is empty and the validation decides.

All the conversion, grafting and validation logic is submission 23's, imported rather
than copied so the two runs cannot silently diverge; only the paths differ.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "submissions", "23_esmfold2_default", "scripts",
                   "02_convert_validate.py")


def load():
    """Import submission 23's converter and repoint it at this submission's folders."""
    spec = importlib.util.spec_from_file_location("conv23", SRC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.RAW = os.path.join(ROOT, "submissions", "23_esmfold2_default",
                         "max_accuracy", "raw")
    m.OUT = os.path.join(HERE, "submission")
    m.RES = os.path.join(HERE, "results")
    m.WORK = os.path.join(HERE, "work")
    m.ZIP = os.path.join(HERE, "submission_esmfold2_maxacc.zip")
    # Nothing is borrowed from another submission: this set is pure max-accuracy unless
    # validation says otherwise.
    m.SWAP_FROM_19 = set()
    return m


if __name__ == "__main__":
    sys.exit(load().main())
