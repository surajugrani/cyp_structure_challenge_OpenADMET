# Submission 31 — ESMFold2 with the MSA subsample depth doubled

One changed field against submission 23, our best single model at 0.4527:
`FoldingConfig.msa_max_depth`, 1024 → 2048. Same `fold_all_atom` endpoint, same
esmfold2-2026-05 weights, same sequence and chain layout, same 14301-sequence
ColabFold MSA, every other field at its API default.

## Why this knob

Analysis 31 probed the whole config against the server. `num_sampling_steps` is
capped at 100 and `num_loops` at its default, so **neither compute knob can be
raised** — the public endpoint already runs at maximum effort. What is left are
the three noise fields, and depth is the only one that moves in both directions.

## Why it is worth a submission

Analysis 31 also ran the control this project had never done: a second run at
settings identical to submission 23. The two agree to a **median 0.31 Å**
ligand RMSD, so ESMFold2 is near-reproducible despite `lm_dropout=0.3` and a
fresh random MSA draw each loop.

Against that floor, depth 2048 sits **1.47 Å** from submission 23 and **1.48 Å**
from the replicate — the same distance from both, which is the signature of a
real effect rather than noise. 14 of 20 ligands move more than 1 Å.

So the poses have genuinely moved, by a controlled amount, from a single dial.
Better or worse is what the leaderboard decides.

## Build

Submission 23's pipeline exactly, imported and repointed rather than copied:
mmCIF → PDB, PrepWizard (`-rmsd 0.3`), re-convert to graft the minimised
coordinates, PoseBusters, zip. **20/20 structural, 20/20 PoseBusters.**

PrepWizard repaired four `minimum_distance_to_protein` clashes that were new at
this depth — slightly tighter packing, consistent with poses having moved.

`OCNT-2313038` is taken from submission 19, exactly as submission 23 does and
for the same reason: ESMFold2 builds that ligand with a non-flat non-aromatic
ring, PoseBusters rejects it, and PrepWizard cannot fix it (submission 23
measured that re-minimising at `-rmsd 1.0` moved it 0.19 Å mean and left the
ring flagged). The failure reproduces at depth 2048, so it belongs to the
weights, not the depth. Importing the same 20th structure submission 23
imported keeps this a one-variable comparison.

## Fe geometry

18 of 20 Fe-ligand distances within 0.3 Å of submission 23, median 2.29 against
2.24 Å — essentially unchanged while ligand RMSD moved 1.47 Å. Fe distance is
one scalar and cannot see rotation about the coordinating atom; use ligand RMSD
to ask whether predictions differ.

Both long-standing oddities reproduce for a fourth independent configuration:
OCNT-2315379 at 7.16 Å (6.89 in sub 23) and OCNT-2395460 at 1.78 Å (1.79).

## What the result decides

* **Above 0.4527** — the dial improves the model, and the depth series becomes
  the diversity generator for Glide `inplace` selection, the one rule measured
  to work within a model (+0.0111, submission 29).
* **Near 0.4527** — depth moves poses without changing quality, which is ideal
  for selection: free diversity at no cost in accuracy.
* **Below 0.4527** — 1024 is already the model's optimum and the dial only
  degrades it; drop the depth series.
