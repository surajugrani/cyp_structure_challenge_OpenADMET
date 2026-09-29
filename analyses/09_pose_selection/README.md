# 09 — Which scoring function identifies the native pose?

**Question:** For the pose-prediction competition, the metric that matters is whether the
*top-ranked* pose is native-like. Analyses 01 and 08 showed no scoring function predicts
affinity; this asks the different and more relevant question — can any of them **pick the
right pose**?

**Status:** Steps 1-3 complete. No single ICM score, and no pooled linear combination, beats ICM's default.

## Subset: the 9 most rigid ligands

The 22 affinity-curated structures ranked by rotatable bonds give
3, 4, 4, 5, 5, 7, 8, **10, 10**, 11, 12, … — a natural gap between 8 and 10, so the 9 most
rigid need no tie-breaking:

| comp_id | pdb_id | RotB |
|---|---|---|
| MYT | 6MA6 | 3 |
| CL6 | 8SPD | 4 |
| YNV | 7LXL | 4 |
| 08Y | 3UA1 | 5 |
| TPF | 6MA7 | 5 |
| KLN | 2V0M | 7 |
| PK9 | 4D6Z | 8 |
| J9K | 4D78 | 10 |
| PKT | 4D7D | 10 |

Rigid ligands because analysis 07 found success falls off steeply with flexibility, and a
pose selector cannot help where sampling never produced a near-native pose.

Poses are reused unchanged: **Glide** from analysis 05 (docked from SMILES — the realistic
blind-set protocol), **ICM** from analysis 06 (ICM's own prep, `confs=10`).

## Step 1 — baseline success rates (each program's own scoring function)

| Program | top-1 < 2 Å | best-of-N < 2 Å | top-1 < 2.5 Å | best-of-N < 2.5 Å |
|---|---|---|---|---|
| Glide | **2/9 (22 %)** | 3/9 (33 %) | 2/9 (22 %) | 5/9 (56 %) |
| ICM | **3/9 (33 %)** | 5/9 (56 %) | 3/9 (33 %) | 6/9 (67 %) |

`best-of-N` is the oracle ceiling — what a perfect selector could reach from the poses
already generated. Pose counts are near-identical between programs here, so correcting for
pose-count parity (`best-of-k`, k = min) changes nothing.

### Per structure

| comp_id | pdb_id | RotB | Glide top-1 | Glide best (rank) | poses < 2 Å | ICM top-1 | ICM best (rank) | poses < 2 Å |
|---|---|---|---|---|---|---|---|---|
| MYT | 6MA6 | 3 | 3.48 | 0.27 (2) | 2/10 | 4.61 | 1.85 (3) | 1/10 |
| CL6 | 8SPD | 4 | 5.33 | 3.34 (8) | 0/10 | **0.40** | 0.40 (1) | 1/10 |
| YNV | 7LXL | 4 | 6.07 | 2.25 (8) | 0/10 | 4.37 | 1.19 (2) | 1/10 |
| 08Y | 3UA1 | 5 | **0.56** | 0.53 (2) | 13/13 | **0.66** | 0.66 (1) | 2/20 |
| TPF | 6MA7 | 5 | 4.68 | 3.40 (2) | 0/9 | 4.05 | 2.21 (7) | 0/10 |
| KLN | 2V0M | 7 | **0.82** | 0.58 (2) | 18/19 | **0.53** | 0.33 (2) | 3/20 |
| PK9 | 4D6Z | 8 | 5.84 | 5.00 (4) | 0/10 | 8.38 | 2.81 (5) | 0/10 |
| J9K | 4D78 | 10 | 5.89 | 2.06 (9) | 0/10 | 4.12 | 3.82 (7) | 0/10 |
| PKT | 4D7D | 10 | 6.64 | 5.41 (6) | 0/10 | 5.98 | 5.98 (1) | 0/10 |

RMSD in Å; `(rank)` is the score rank of the lowest-RMSD pose.

## What this says about where a selector can help

**The two programs fail differently, and only one failure mode is fixable by rescoring.**

- **Glide is bimodal.** When it finds the site it finds it overwhelmingly — 13/13 and 18/19
  poses under 2 Å — and its own score already ranks one first. When it misses, **zero** of
  its poses are near-native. There is essentially no case where Glide generated a good pose
  and mis-ranked it, so a better selector on Glide's own poses buys at most 1/9.
- **ICM is the opposite.** It generates a near-native pose in 5/9 cases but only 1–3 of 10
  poses are near-native, and it ranks that pose first in only 2 of those 5. The oracle
  ceiling is 5/9 vs 3/9 delivered — **a perfect selector on ICM's poses gains 2/9**, and the
  mis-ranked cases (6MA6 rank 3, 7LXL rank 2) are near-misses.

So ICM is the right pose source for the selector work, as planned, and the realistic prize
is 33 % → 56 %. With n = 9 that is a 2-structure difference and cannot be called
significant; it sets direction, not proof.

**Caveat on 4D7D and 4D78:** both RotB = 10 and neither program comes close. They will
dominate any "failure" count while telling us nothing about scoring.

## Step 2a — can RTCNN rescue ICM's mis-ranked poses?

**No runs were needed.** `_dockScan` had already written `RTCNN` and `mfScore` into every
pose in analysis 06's `<pdb>_hits.sdf`; analysis 06's parser dropped RTCNN because it
filtered on a tag called `RTCNNscore` when the real SD tag is `RTCNN`. `scripts/02_rtcnn_rescore.py`
re-reads those files and re-ranks the existing poses.

Set: all RotB <= 8 structures ICM docked (33), less 3 whose pose SDFs were lost when their
structure directories were rebuilt (5A1R, 1W0F, 5A1P) -> **n = 30**.

| Ranking score | top-1 < 2 Å | top-1 < 2.5 Å |
|---|---|---|
| ICM `Score` (default) | **8/30 (27 %)** | 8/30 (27 %) |
| `RTCNN` | 6/30 (20 %) | 7/30 (23 %) |
| `mfScore` | 3/30 (10 %) | 3/30 (10 %) |
| oracle (best pose available) | 20/30 (67 %) | 21/30 (70 %) |

**RTCNN is not an improvement — it is slightly worse.** Head to head against `Score` at
2 Å it **saves 1 and breaks 3**:

| | pdb | ligand | RotB | Score top-1 | RTCNN top-1 | oracle |
|---|---|---|---|---|---|---|
| saved | 6MA8 | PMF | 2 | 7.22 | **0.53** | 0.53 |
| broken | 5VCG | 08Y | 5 | **1.02** | 3.44 | 1.02 |
| broken | 5TE8 | 08J | 1 | **0.53** | 3.65 | 0.52 |
| broken | 4NY4 | 2QH | 4 | **0.92** | 8.01 | 0.92 |

`mfScore` is clearly worst overall, but it is the only one of the three that gets 7LXL
(1.19 Å, ranked first) and 9BV6 (1.56 Å) right. So the three scores fail on *different*
structures — which is the argument for a combination rather than for swapping one in.

The gap that matters is unchanged: **8/30 delivered against a 20/30 ceiling.** ICM samples
a near-native pose two-thirds of the time and ranks it first a quarter of the time, and
none of its three built-in scores closes that gap.

### Caveat

`Score`, `RTCNN` and `mfScore` were all compared on the same fixed pose sets, so this is a
fair ranking comparison. But with n = 30 the difference between 8 and 6 successes is two
structures and is not significant; the safe reading is "RTCNN is not better", not "RTCNN
is worse".

## Step 3 — a pooled linear combination of the three scores

Idea: since the three scores succeed on different structures, fit
`predicted_RMSD = w1*Score + w2*RTCNN + w3*mfScore + b` over all 469 poses and rank by the
prediction. Evaluated **leave-one-structure-out** (fit on 29 structures, rank the held-out
one's poses) because poses of one ligand are highly correlated — effective n is 30
structures, not 469 poses.

Two training sets were tried, since 10 of the 30 structures have no pose under 2 A and so
contribute only negative examples. The **test set is always all 30** — on the blind set you
cannot know which ligands have a good pose available, so restricting evaluation would
inflate the result (oracle would be 20/20 by construction).

| Ranker | top-1 < 2 Å |
|---|---|
| ICM `Score` alone | **8/30 (27 %)** |
| `RTCNN` alone | 6/30 |
| `mfScore` alone | 3/30 |
| linear, Score+RTCNN+mfScore, train all 30 | 6/30 |
| linear, Score+RTCNN+mfScore, train solvable 20 | 5/30 |
| linear, Score+RTCNN, train all 30 | 6/30 |
| linear, Score+RTCNN, train solvable 20 | 6/30 |
| oracle | 20/30 (67 %) |
| random pose | 2.0/30 (7 %) |

**Every fitted model is worse than simply using ICM's Score.** Dropping `mfScore` does not
help: its fitted weight is −0.0019, two orders of magnitude below the others, so the
pooled fit was already ignoring it.

### Why: Simpson's paradox

The fitted weight on `Score` is **negative** (−0.056) — the model concludes that a *better*
ICM score means a *worse* pose. Checking the correlations directly explains it:

| feature | pooled ρ vs RMSD | median within-ligand ρ |
|---|---|---|
| `Score` | **−0.202** | +0.017 |
| `RTCNN` | +0.136 | **+0.334** (positive in 20/30) |
| `mfScore` | +0.023 | −0.065 |

`Score`'s pooled and within-ligand correlations have **opposite signs**. The cause is a
between-ligand confound: across structures, mean `Score` vs oracle RMSD gives ρ = −0.305,
i.e. ligands that score *less* favourably actually dock *better* (they are the smaller,
more rigid ones). That between-ligand trend is an order of magnitude stronger than any
within-ligand signal, so a pooled regression fits it and inverts the weight.

The model is not picking at random — trained on all 30 it saves 4 structures (1W0G
8.18→1.29, 6MA8 7.22→0.53, 7LXL 4.37→1.19, 8SO2 3.69→0.54, all landing exactly on the
oracle pose) — but it breaks 6 that `Score` already had right.

**Conclusion: pooling across ligands cannot work for this task.** The selector must be
within-ligand (features and target standardised inside each structure), so the
ligand-level offset is removed and only pose-to-pose differences remain. Note also that
`RTCNN` has by far the best within-ligand signal (+0.334) despite a worse top-1 than
`Score`, so it is the term to build on.

## Step 4 (in progress)

Generating a `_dockScan -n` near-native pose for all 30 structures, to give every ligand a
positive example. **Recorded confound:** `-n` scores the same geometry 0.25–13.4 Score
units better than `_dockScan` does (3UA1: 9.3 units better at 0.66 Å separation), and the
offset is not constant, so a model pooling docked and `-n` rows can identify the `-n` row
from its score alone. See `scripts/04_native_poses.py`.



Test whether any individual score term — Glide's or ICM's — ranks ICM's 10 poses better
than ICM's total score. Must be a **within-ligand** analysis with a structure-level split.
See [[pose-selector-idea]].

## Rerun

```bash
$SCHRODINGER/run python3 scripts/01_baseline_success.py
```

Light: pandas only, no licence, no docking. Outputs `results/rigid9_baseline.csv` and
`results/rigid9_summary.csv`.
