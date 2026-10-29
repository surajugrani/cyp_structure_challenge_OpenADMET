"""Glide `inplace` scores for the depth-2048 ESMFold2 set, to make a four-set pool.

Submission 33 tested the pruning direction: drop set 26 from submission 29's three-set
pool and select over the remaining two. It lost 0.0174, and landed below plain ESMFold2.
Set 26 was over-picked relative to chance AND holds the worst LDDT-PLI of the three, and
removing it still cost - so those picks were carrying real value and pool size, not
candidate quality, is what the selector is living on.

This runs the converse. Add a fourth candidate instead of removing one, and make it the
sharpest version of the test available: set 31, msa_max_depth 2048, which scored 0.4116
on its own - the worst ESMFold2 run of the project and 0.0411 below set 23. If pool size
really is the lever, a candidate this bad should STILL help.

Two outcomes, both informative:
  * It helps. Selection within one model scales with the pool, independently of how good
    each candidate is on its own, and the next move is to fold more depth variants purely
    as pool fodder rather than as submissions.
  * It hurts. Submission 33's loss was about set 26 specifically, not about pool size,
    and the selector does depend on candidate quality after all.

`inplace` only. `mininplace` is the mode submission 29 identified as the confound behind
every failed Glide selection in this project, and analysis 28 has the full both-mode
comparison already; there is nothing left to learn from running it again here.

The scoring code is analysis 28's, imported and repointed rather than copied, so the
four sets cannot be scored by two subtly different protocols. Same PrepWizard flags,
same Fe-ligand bond break, same grid generation, same Glide input.

NOTE on OCNT-2313038: set 31 ships submission 19's structure for this ligand, exactly as
set 23 does, because ESMFold2 builds its ring non-flat and PoseBusters rejects it. So for
that one ligand the "new" candidate is identical to one already in the pool and cannot
change the pick. 19 of 20 ligands get a genuine fourth option.
"""
import importlib.util, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "scripts",
                   "01_score_both.py")

spec = importlib.util.spec_from_file_location("a28", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

m.WORK = os.path.join(HERE, "work")
m.OUT_CSV = os.path.join(HERE, "results", "inplace_depth2048.csv")
m.SETS = {"31_depth2048": "31_esmfold2_depth2048"}
m.MODES = ["inplace"]

if __name__ == "__main__":
    sys.exit(m.main())
