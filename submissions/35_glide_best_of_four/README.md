# Submission 35 - best Glide `inplace` over FOUR ESMFold2 runs

The converse of submission 33. That one removed a candidate from submission 29's
three-run pool and lost 0.0174; this adds one.

## The pool

| set | run | own LDDT-PLI | picks here |
|---|---|---|---|
| 23 | default endpoint + ColabFold MSA | 0.4527 | 3 |
| 26 | max_accuracy endpoint, own MSA | 0.4365 | 6 |
| 27 | default endpoint, no MSA | 0.4384 | 4 |
| **31** | default endpoint, `msa_max_depth` 2048 | **0.4116** | **7** |

Set 31 is the worst ESMFold2 run of the project, 0.0411 below set 23. It was added
precisely because it is the worst: if pool size is what the selector is living on, a
candidate this weak should still help.

## What changed

**7 of 20 structures differ from submission 29, and set 31 won every one of them.**

| ligand | submission 29 | here | gscore gained |
|---|---|---|---|
| OCNT-2395328 | 26_max_accuracy | 31_depth2048 | +1.502 |
| OCNT-2395460 | 27_default_nomsa | 31_depth2048 | +1.285 |
| OCNT-2312966 | 26_max_accuracy | 31_depth2048 | +0.619 |
| OCNT-2313397 | 23_default_msa | 31_depth2048 | +0.469 |
| OCNT-2395464 | 23_default_msa | 31_depth2048 | +0.316 |
| OCNT-2395436 | 23_default_msa | 31_depth2048 | +0.269 |
| OCNT-2535519 | 26_max_accuracy | 31_depth2048 | +0.109 |

Set 31 takes 7 of 20 against the 5 a fourth equal candidate would take by chance. That
is the same over-selection pattern set 26 showed in submission 29 (9 of 20), and
submission 33 established that removing such a set costs rather than helps.

No claim is made that these 7 structures are better. gscore ranks candidates within one
receptor; whether that tracks LDDT-PLI is the question the leaderboard answers.

## What this cannot settle

The leaderboard scores the same 10 of our 20 structures every time, and we do not know
which 10. If the 7 changed picks fall mostly in the unscored half, the result comes back
close to submission 29's and means nothing either way. A null result here is
uninformative, not support for the hypothesis.

More generally, every delta this project has reasoned from is a mean over those same 10
ligands. Differences of +-0.01 rest on ten structures.

## Reading the result

- **Above 0.4638** - pool size is the lever. Candidate quality barely matters, and the
  next move is to fold more configuration variants purely as pool fodder rather than
  spending submissions on them individually.
- **Below 0.4638** - adding a weak candidate hurts, so submission 33's loss was about
  set 26 specifically rather than about pool size, and the pool should be built from
  comparable-quality runs only.
- **Indistinguishable from 0.4638** - most likely the changes landed in the unscored
  half; nothing follows.

## Note on OCNT-2313038

Sets 23 and 31 both ship submission 19's structure for this ligand, because ESMFold2
builds its ring non-flat and PoseBusters rejects it. For that one ligand the fourth
candidate duplicates one already in the pool and cannot change the pick, so 19 of 20
ligands received a genuine fourth option.

## Validation

**COPY-ONLY.** Every structure is copied unchanged from a submission that already
passed the full gate 20/20, verified by sha256: 20/20 matched. PoseBusters is not
re-run - a byte-identical copy cannot return a different verdict.

## A note on the run

The last two ligands could not be scored for about an hour: the shared Schrodinger
`GLIDE_MAIN` pool was exhausted by other users on the cluster and Glide could not check
out a token. The grids were already built, so the work itself was two minutes. The build
script refuses to assemble on partial scores rather than silently selecting from three
sets for the missing ligands - the failure mode that produced a wrong hybrid on the
first attempt at submission 30.
