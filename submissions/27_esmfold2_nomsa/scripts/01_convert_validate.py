"""Submission 27: the default ESMFold2 endpoint with NO MSA - the clean MSA control.

Submission 23 is `fold_all_atom` + our 14,301-sequence ColabFold MSA, 0.4527 LDDT-PLI,
the best result of the project. Submission 26 was the max-accuracy endpoint, which
builds its own MSA, and scored 0.4365 - but it also changed the endpoint and the model
weights, so it could not say what the MSA was worth.

This changes ONE thing. Same `fold_all_atom` endpoint, same esmfold2-2026-05 weights,
same sequence and chain layout as submission 23, `msa=None`. The difference against 23
is the ColabFold MSA and nothing else.

The Fe geometry says the MSA is doing almost nothing: 18 of 20 structures land within
0.3 A of submission 23, median Fe distance 2.20 against 2.24 A, coordination unchanged
at 12 of the 15 ligands that can form it. It also reproduces both of the values that
looked impossible in 23 - OCNT-2315379 at 6.91 A against 6.89, OCNT-2395460 at 1.75 A
against 1.79 - so those are properties of the endpoint and its weights, not artefacts
of the MSA we supplied.

That makes this worth scoring rather than filing away. If it lands near 0.4527 the MSA
is costing us an MSA build for nothing; if it lands well below, the MSA is buying
something the Fe distances cannot see.

Conversion, grafting and validation are submission 23's, imported and repointed rather
than copied so the runs cannot silently diverge. SWAP_FROM_19 is empty: nothing is taken
from another submission unless validation forces it.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "submissions", "23_esmfold2_default", "scripts",
                   "02_convert_validate.py")


def load():
    spec = importlib.util.spec_from_file_location("conv23", SRC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.RAW = os.path.join(HERE, "raw")
    m.OUT = os.path.join(HERE, "submission")
    m.RES = os.path.join(HERE, "results")
    m.WORK = os.path.join(HERE, "work")
    m.ZIP = os.path.join(HERE, "submission_esmfold2_nomsa.zip")
    m.SWAP_FROM_19 = set()
    return m


if __name__ == "__main__":
    sys.exit(load().main())
