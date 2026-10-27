# Analysis 31 — what the ESMFold2 config knobs actually do

Submissions 23, 26 and 27 are three ESMFold2 sets that disagree, and submission 29 won
by selecting among them. That raised a question the project had never answered: where
does the disagreement come from, and can we generate more of it on purpose?

## 1. Both compute knobs are already at maximum

Probed directly against the server, one field at a time:

| field | default | tried | server |
|---|---|---|---|
| `num_sampling_steps` | 100 | 200 | rejected — "must be ≤ 100" |
| `num_loops` | 20 | 40, 100 | rejected, both over the cap |
| `lm_dropout` | 0.3 | 0.0 | accepted |
| `msa_max_depth` | 1024 | 2048, None | accepted |

The public `fold_all_atom` endpoint will not let you spend more compute than the default
already does. There is no "run it harder" setting. That also explains what
`fold_max_accuracy` is for — it is a separate endpoint with no knobs at all, because the
knobs on this one are already maxed; its accuracy comes from pinned weights, not from
more sampling.

## 2. ESMFold2 is near-reproducible — the control nobody had run

`lm_dropout` defaults to 0.3 and a fresh random 1024 MSA rows are drawn every loop, so
two runs at identical settings ought to differ. They barely do.

A straight replicate of submission 23, same endpoint, same weights, same MSA, no flags,
sits a **median 0.31 Å** pocket-superposed ligand RMSD from it. 1 of 20 ligands moves
more than 1 Å; none moves more than 1.2 Å.

So the three existing ESM sets differ because they used different *configurations*, not
because the endpoint is noisy. Repeating a run does not generate candidates.

## 3. MSA depth is a real dial

Doubling `msa_max_depth` to 2048, one changed field:

| comparison | median | mean | >1 Å | >2 Å |
|---|---|---|---|---|
| **replicate vs sub 23** (control) | **0.31 Å** | 0.34 | 1/20 | 0/20 |
| depth 2048 vs sub 23 | 1.47 Å | 2.24 | 14/20 | 7/20 |
| depth 2048 vs replicate | 1.48 Å | 2.12 | 13/20 | 7/20 |
| depth 2048 vs no-MSA (27) | 1.71 Å | 2.27 | 14/20 | 9/20 |
| depth 2048 vs max-accuracy (26) | 2.46 Å | 3.18 | 17/20 | 12/20 |

Depth 2048 sits ~1.47 Å from submission 23 and ~1.48 Å from the replicate — the *same*
distance from both. That is what a real effect looks like: the replicate and submission
23 are interchangeable as references, and the knob displaces poses about five times the
0.31 Å noise floor.

Costs nothing in time, either: 29 s per fold, identical to the default. Disabling
subsampling entirely (`None`) costs 48 s and changes the mechanism rather than the
number, since it flips `msa_subsample_at_inference` off — which is why 2048 is the
cleaner single-variable test and `None` was abandoned part-way.

## 4. Fe distance stays blind, again

Against submission 23 the depth-2048 set keeps 18 of 20 Fe-ligand distances within
0.3 Å (median 2.29 against 2.24 Å) while ligand RMSD moves a median 1.47 Å. One scalar
to the iron cannot see rotation about the coordinating atom. Same lesson as analysis 27;
use ligand RMSD to ask whether predictions differ.

Both long-standing oddities reproduce for a fourth independent configuration —
OCNT-2315379 at 7.16 Å (6.89 in sub 23) and OCNT-2395460 at 1.78 Å (1.79) — so they are
properties of the esmfold2-2026-05 weights, not of any MSA treatment.

## What this means for selection

Submission 29 showed Glide `inplace` selection works within a model (+0.0111). The
material it needs is a diverse candidate pool, and section 2 says rerunning will not
produce one. Section 3 says the depth dial will: each value yields a distinct structure
~1.5 Å from its neighbours, reproducibly, at no extra cost.

Whether that is *useful* diversity is what submission 31 tests. Pending: fold the series
below 1024 and above 2048 only if 31 shows the dial does not simply degrade the model.

## Files

* `scripts/01_convert.py` — mmCIF → PDB, importing submission 23's converter
* `scripts/run.sh`, `scripts/run_replicate.sh` — the two RMSD comparison jobs
* `results/vs_*.csv`, `results/replicate_vs_*.csv` — per-ligand RMSDs
* `replicate/` — the control run's 20 mmCIFs and fold metrics
