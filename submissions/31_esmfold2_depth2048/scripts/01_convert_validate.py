"""Submission 31: ESMFold2 with the MSA subsample depth doubled, 1024 -> 2048.

One changed field against submission 23 (0.4527, still our best single model):
`FoldingConfig.msa_max_depth`. Everything else is submission 23 exactly - same
`fold_all_atom` endpoint, same esmfold2-2026-05 weights, same sequence, same chain
layout, same 14,301-sequence ColabFold MSA, every other config field left at its API
default.

Why this field and not another. Analysis 31 probed the whole config against the server:
`num_sampling_steps` is capped at 100 and `num_loops` at its default, so neither compute
knob can be raised - the public endpoint already runs at maximum effort. What remains
are the three noise fields, and depth is the one that can move in both directions.

Why it is worth a submission. Analysis 31 also ran the control this project had never
done: a second run at settings identical to submission 23. The two agree to a median
0.31 A pocket-superposed ligand RMSD, so ESMFold2 is near-reproducible despite
`lm_dropout=0.3` and a fresh random MSA draw every loop. Against that 0.31 A floor,
depth 2048 sits 1.47 A from submission 23 and 1.48 A from the replicate - the same
distance from both, which is what a real effect looks like rather than noise. 13-14 of
the 20 ligands move more than 1 A.

So this asks a question the project has not been able to ask before: the poses have
genuinely moved, and moved by a controlled amount from a single dial. Better or worse
is what the leaderboard decides. If better, the depth series below 1024 and above 2048
becomes the diversity generator for Glide `inplace` selection, the one rule measured to
work within a model (+0.0111, submission 29). If worse, 1024 is already the model's
optimum and the dial only degrades it.

Conversion, PrepWizard grafting and validation are submission 23's, imported and
repointed rather than copied so the two runs cannot silently diverge.

SWAP_FROM_19 carries OCNT-2313038, exactly as submission 23 does, and for the same
reason: ESMFold2 builds that ligand with a non-flat non-aromatic ring and PoseBusters
rejects it. The failure reproduces here at depth 2048, so it is a property of the
weights rather than of the MSA depth, and PrepWizard does not fix it - submission 23
measured that re-minimising at -rmsd 1.0 moved the ligand 0.19 A mean and left the ring
flagged. Importing the same structure submission 23 imported keeps this a one-variable
comparison: 19 ESMFold2 structures differing from 23 only in msa_max_depth, plus the
identical 20th.
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
    m.ZIP = os.path.join(HERE, "submission_esmfold2_depth2048.zip")
    m.SWAP_FROM_19 = {"OCNT-2313038"}
    return m


if __name__ == "__main__":
    sys.exit(load().main())
