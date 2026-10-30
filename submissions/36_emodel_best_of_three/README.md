# Submission 36 - best `emodel` over the three ESMFold2 runs

Same pool as submission 29, same structures, **different ranking term.**

## Why change the term rather than the pool

Submission 29 selected on `gscore` under `inplace` and scored 0.4638, the project best.
Both attempts to improve on it by changing the pool have failed, in opposite directions:

| | change | result |
|---|---|---|
| 33 | remove set 26 (3 -> 2 runs) | -0.0174 |
| 29 | the pool as built (3 runs) | **0.4638** |
| 35 | add set 31 (3 -> 4 runs) | -0.0126 |

The pool is a local optimum. What has never been varied is the quantity being ranked.

## Why `emodel` specifically

`gscore` is Glide's estimate of binding affinity - it answers "how well does this ligand
bind this receptor". `emodel` is Glide's **pose-ranking** term, built to choose among
alternative poses of the *same* ligand, weighting force-field interaction energy far
more heavily.

Choosing among alternative structures of one ligand is exactly our task. `emodel` is the
term designed for it, and this project has used it for nothing.

It was also paid for and thrown away. The original reason for preferring
`DOCKING_METHOD mininplace` over `inplace` was that `inplace` returns `emodel` as 0.00 -
and then every selection built on `mininplace` ranked by `gscore` anyway. Analysis 28
scored all 60 structures in both modes, so these numbers have been on disk since.

## The objection, and why it no longer holds

`emodel` describes the **minimised** pose, while the structure shipped is the
unminimised one. The quantity ranked is not a property of the file in the zip. That
mismatch is what made submission 22 uninterpretable: it ranked minimised poses and
shipped raw ones.

Submission 34 settled it. Every ligand in submission 29 was replaced by its minimised
pose - median 0.497 A displacement, max 2.503 A - and the result was 0.4637 against
0.4638, with LDDT-LP identical to four decimals. **As far as the metric is concerned
the minimised and unminimised poses are the same structure**, so ranking one and
shipping the other is no longer a meaningful mismatch.

## What changed

**8 of 20 structures differ from submission 29. The two terms agree on only 12.**

| ligand | submission 29 (`gscore`) | here (`emodel`) |
|---|---|---|
| OCNT-2313397 | 23_default_msa | 26_max_accuracy |
| OCNT-2313429 | 27_default_nomsa | 23_default_msa |
| OCNT-2395328 | 26_max_accuracy | 27_default_nomsa |
| OCNT-2395342 | 26_max_accuracy | 27_default_nomsa |
| OCNT-2395460 | 27_default_nomsa | 23_default_msa |
| OCNT-2395464 | 23_default_msa | 27_default_nomsa |
| OCNT-2395503 | 26_max_accuracy | 23_default_msa |
| OCNT-2395645 | 27_default_nomsa | 26_max_accuracy |

Picks come out **23x7, 26x8, 27x5** - more even than submission 29's 6/9/5, and far
more even than submission 35's skew toward its weakest member. Set 26's
over-representation, the thing submission 33 was built to test, is slightly smaller
here without being targeted.

No quality claim is made by any `emodel` value above. It is a scoring term; the
leaderboard is the only measurement of quality this project has.

## Reading the result

- **Above 0.4638** - the selector, not the pool, was the limiting factor all along, and
  `gscore` was simply the wrong term for a pose-choice problem. That reopens the whole
  selection direction with a better instrument.
- **Below 0.4638** - `gscore` is the better of the two for this task, and three
  independent routes to improving on submission 29 (prune, extend, re-rank) have failed.
- **Indistinguishable** - the two terms disagree on 8 structures but those disagreements
  are not worth anything, which would itself say the signal is thin.

Standing caveat: the leaderboard scores the same 10 of our 20 structures every time and
we do not know which. If the 8 changes fall mostly in the unscored half, nothing
follows.

## Validation

**COPY-ONLY.** Every structure is copied unchanged from a submission that already
passed the full gate 20/20, verified by sha256: 20/20 matched. PoseBusters is not
re-run - a byte-identical copy cannot return a different verdict.
