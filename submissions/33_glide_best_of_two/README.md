# Submission 33 — best Glide `inplace` over TWO ESMFold2 runs (23 + 27)

Submission 29 ran this exact rule over three runs and scored 0.4638, the project
best. This drops one and asks whether the selector was better off without it.

## Which one was dropped, and why

Submission 29's 20 winning picks:

| set | picks | its LDDT-PLI | its LDDT-LP |
|---|---|---|---|
| 23 default + ColabFold MSA | 6 | 0.4527 | 0.7985 |
| **26 max-accuracy** | **9** | **0.4365** (worst) | **0.8063** (best) |
| 27 default, no MSA | 5 | 0.4384 | 0.7908 |

Submission 26 supplies **9 of 20** winners against 6.7 expected by chance, while
having the worst protein–ligand interface of the three — and the best pocket
geometry of any submission in this project.

There is a mechanism for that pairing. Glide gscore measures how well a ligand
fits the receptor it is handed; it cannot judge the receptor. A better-built
pocket therefore earns a better gscore whether or not the ligand is correctly
placed. Submission 30 paid 0.0113 for precisely this blind spot by selecting
across models of unequal pocket quality. Within this pool, submission 26 is the
set with the conspicuously better pocket, so its picks are the ones most likely
to reflect pocket quality rather than pose accuracy.

Dropping it leaves the cleanest pair available: 23 and 27 share the
`fold_all_atom` endpoint, the esmfold2-2026-05 weights, the sequence and the
chain layout, and differ only in whether our ColabFold MSA was supplied.
Comparable candidate quality is the precondition submissions 30 and 31 both
established.

## What it selected

| set | ligands won |
|---|---|
| 23 default + MSA | 9 |
| 27 default, no MSA | 11 |

* changes **9 of 20** against submission 29
* differs from submission 23 alone on **11 of 20**
* gives up a median **0.57** and at most **1.08** gscore units on the ones it
  changes — by construction, since the dropped set held the better score there

## What the result decides

* **At or above 0.4638** — Glide was being fooled by pocket quality, and
  restricting the pool to comparable receptors is the better rule.
* **Below 0.4638** — the pool mattered more than the bias. Removing a third of
  the candidates from a rule that was already working cost more than the bias
  did, and submission 29 stands.

The honest risk is the second one: submission 29 *won* with those nine picks
from set 26, so whatever bias exists was not fatal.

## Build

Copy-only. Every structure is copied unchanged from a submission that already
passed the full gate 20/20, verified by sha256. PoseBusters is not re-run — a
byte-identical copy cannot return a different verdict.

* `results/selection.csv` — per-ligand gscores for both kept sets and the
  dropped one, the pick, and what changed against submission 29
* `visualize/33_glide_best_of_two.pse` — PyMOL session, 20 objects
