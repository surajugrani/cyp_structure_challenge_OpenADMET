# 09 — Which scoring function identifies the native pose?

**Question:** For the pose-prediction competition, the metric that matters is whether the
*top-ranked* pose is native-like. Analyses 01 and 08 showed no scoring function predicts
affinity; this asks the different and more relevant question — can any of them **pick the
right pose**?

**Status:** Steps 1-8 complete. No single ICM score, and no pooled linear combination, beats ICM's default. Near-native candidate rows now exist for all 30 structures in both programs; the within-ligand model is the next step.

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

## Step 4 — a near-native candidate for every ligand

In 10 of the 30 structures docking never sampled a pose under 2 Å, so those ligands
contributed only negative examples. `_dockScan -n` minimises the deposited ligand in place
and scores it through ICM, giving every ligand an 11th, near-native candidate.

Job 47800305, 30/30, no failures. **All 30 `-n` poses are genuinely near-native:** RMSD to
the deposited pose 0.26–1.42 Å, median 0.56 Å, 30/30 under 2 Å. The pool is now
**499 rows over 30 structures** with an oracle of 30/30.

The Glide counterpart is step 5.

### Confound, recorded — do not pool these rows naively

`-n` scores are not on the same scale as `_dockScan` docked-pose scores. On the 7
structures where both already existed, `-n` scores the SAME geometry **0.25–13.4 Score
units better** (3UA1: 9.3 units better at only 0.66 Å separation), against a within-ligand
pose-score spread of only ~10–20 units. The offset is not constant, so it cannot be
subtracted away.

A model trained on a pooled set of docked + `-n` rows can therefore identify the `-n` row
from its score alone, score 30/30 in training, and be worthless on a blind set where no
row carries that signature. Any use of these rows has to control for it — by z-scoring
within ligand, by holding the `-n` row out of evaluation while keeping it in training, or
at minimum by checking whether the fitted model simply always picks the `-n` row.

## Step 5 — Glide crystal-pose scores, and why Glide has no route artefact

Glide `mininplace` on all 30 crystal poses (job 47800325, 30/30, all `Done`), the Glide
counterpart of step 4. See `scripts/05_glide_crystal.py` for why `mininplace` and not
`inplace`, why not both, and why all 30 were re-run rather than the 23 analysis 01 lacked.

**The confound that threatens the ICM pool does not appear for Glide.** Offset between the
crystal-pose score and the closest docked pose's score, both on the same grid:

| | median | range |
|---|---|---|
| all 30 | −0.00 | −4.25 to +4.28 |
| where a docked pose is < 1.0 Å away (n = 17) | **−0.06** | −1.51 to +1.04 |

Centred on zero and **not directional** — 9 of 17 negative, 8 positive, 12 of 17 within
±0.3 against a within-ligand spread of ~0.77. Contrast ICM, where `-n` scored the same
geometry better in 7 cases out of 7, by up to 13.4 units. Running all 30 on one grid was
what made this clean; the earlier cross-grid estimate (0.10–0.18, analysis 01 grids vs
analysis 05 poses) was tighter-looking but confounded.

So Glide's crystal rows can be pooled with its docked rows directly, while ICM's `-n` rows
cannot be pooled without controlling for the offset.

**Not verified:** `CALC_INPUT_RMS` was requested but Glide did not emit `rmsd_to_input` in
the results CSV, so unlike the ICM side there is no direct measure of how far `mininplace`
moved each crystal pose. The near-zero offsets are indirect evidence it stayed close, not
a measurement of it.



Test whether any individual score term — Glide's or ICM's — ranks ICM's 10 poses better
than ICM's total score. Must be a **within-ligand** analysis with a structure-level split.
See [[pose-selector-idea]].

## Rerun

```bash
$SCHRODINGER/run python3 scripts/01_baseline_success.py
```

Light: pandas only, no licence, no docking. Outputs `results/rigid9_baseline.csv` and
`results/rigid9_summary.csv`.


## Step 6 — fitting the 11-row ICM pool, and measuring the route artefact

Step 4 added a `-n` near-native row per ligand but recorded a risk: because `-n` scores the
same geometry better than `_dockScan`, a model might just learn "pick the suspiciously
good row", score well, and be useless blind. This step fits the pooled model on the
499-row pool and measures that directly.

| | result |
|---|---|
| top-1 < 2 Å on the full 11-row pool | **4/30** |
| picked the `-n` row | **1/30** — artefact meter |
| top-1 < 2 Å, `-n` row removed at test time | 4/30 |
| ICM `Score` baseline | 8/30 |
| oracle, docked only | 20/30 |

**The artefact did not materialise.** The model picks the `-n` row once in 30, so there is
no inflation to discount. But the pool got *worse*: 4/30, down from 6/30 on 469 rows.

### Why, and a correction to step 4's stated risk

The `-n` rows have better scores (mean −16.40 vs −11.57) *and* lower RMSD (0.61 vs 5.87),
which should push the `Score` weight positive and undo the inversion. It moves it — −0.0562
to −0.0371 — but not across zero. The weight stays negative, so the model still believes a
better score means a worse pose and therefore **actively avoids the near-native rows**.
30 correcting points cannot outweigh 469 confounded ones.

**Correction:** step 4 described the 0.25–13.4 unit offset as comparable to the whole
signal, implying the `-n` row would be trivially identifiable. Measured: the `-n` row is
the best-scoring row in its ligand only **9/30** times, median rank 4. The offset is real
but usually smaller than the within-ligand pose spread.

**Conclusion: the blocker is the pooling, not the training data.** Adding positive examples
cannot help a model whose sign is inverted. The within-ligand model is the next step, and
the `-n` rows should help there, since z-scoring inside each ligand removes exactly the
between-ligand confound holding the weight negative.


## Step 7 — the within-ligand model

One model, not 30. The within-ligand step is a **feature transform**: z-score each ligand's
rows against each other (using only that ligand's own scores, no labels), then fit a single
global set of weights. The ligand-level offset that inverted the pooled fit is removed;
the pose-to-pose differences that ranking depends on survive.

Set: 11 rows per ligand — the best 10 docked poses by Score plus the `-n` row (329 rows;
8SO2 has only 9 docked poses). Capping at 10 docked rows matters: uncapped, ligands with
20–30 poses dominated and performance halved. The `-n` row is included both in training and
as a test candidate.

| Model (LOSO) | top-1 < 2 Å | picks `-n` |
|---|---|---|
| Score + RTCNN, target = z-RMSD | **11/30** | 8/30 |
| Score + RTCNN + mfScore, target = raw | 10/30 | 5/30 |
| Score + RTCNN + mfScore, target = z | 10/30 | 6/30 |
| Score + RTCNN, target = raw | 6/30 | 1/30 |

**Baselines on the same 11 candidates:** ICM `Score` **12/30**, `RTCNN` 11/30,
`mfScore` 7/30, oracle 30/30.

### The transform worked; the model still does not beat the baseline

Within-ligand Spearman vs RMSD is now positive for all three terms — RTCNN **+0.459**
(was +0.334), Score +0.127 (was +0.017), mfScore +0.082 (was −0.065) — and the pooled
models' inversion is gone. Performance rose from 4/30 (pooled, 499 rows) to 11/30. Step 3's
diagnosis was correct.

But the fair baseline is **12/30, not the 8/30 quoted earlier in this README**: adding the
`-n` row as a candidate lifts plain `Score` from 8/30 to 12/30 by itself, because it is the
best-scoring row in 9/30 ligands and is always near-native. Most of the 8 → 11 improvement
is the extra candidate, not the model. The best model is still one structure behind Score.

The fitted `Score` weight stays slightly negative (−0.027) even z-scored, while RTCNN
carries the signal (+0.288) — the model is effectively "RTCNN with a small Score penalty",
which is why it lands next to RTCNN-alone's 11/30.

**Target choice matters more than feature count.** Score+RTCNN gets 6/30 with raw RMSD and
11/30 with z-scored RMSD — same features, same rows. Raw RMSD lets ligands whose poses are
all 5–11 Å dominate the least-squares residuals.

**Nothing tried so far — three features, two targets, two pooling schemes — beats ICM's own
Score.**


## Step 8 — is it the model or the features?

| | top-1 < 2 Å |
|---|---|
| **linear ceiling, global weights chosen with hindsight (in-sample)** | **14/30** |
| linear, LOSO | 11/30 |
| GradientBoosting, LOSO | 12/30 |
| RandomForest, LOSO | 9/30 |
| ICM `Score` baseline | 12/30 |
| oracle | 30/30 |

**It is a feature limitation.** A global linear combination of the three ICM scores reaches
only **14/30 even when the weights are fitted on the test data**. The cross-validated model
already gets 11/30, within 3 structures of that cheating ceiling — the missing 16 are not
lost to a bad fit, the features do not contain the information. Non-linear models agree:
GradientBoosting 12/30 (ties the baseline, below the linear ceiling), RandomForest 9/30.

**5/30 are provably unreachable** by any model of these features — in 1W0G, 4D75, 4NY4,
6OOB and 9BVC the near-native pose ranks 4th–11th on *all three* terms at once.

### Which additional ICM terms are worth adding

All per-pose tags in the hitlist, ranked by within-ligand Spearman vs RMSD:

| term | ρ | positive in | top-1 alone | top-1, docked rows only |
|---|---|---|---|---|
| RTCNN | +0.459 | 23/30 | 11/30 | 7/30 |
| **Hbond** | **+0.345** | 22/29 | 9/30 | 6/30 |
| Score | +0.127 | 18/30 | 12/30 | 8/30 |
| SolEl | +0.118 | 18/30 | 5/30 | — |
| Dsolv | +0.100 | 17/30 | 7/30 | — |
| mfScore | +0.082 | 16/30 | 7/30 | 5/30 |
| Eintl | −0.027 | 14/30 | 14/30 ⚠ | 2/30 |
| VwInt | −0.053 | 14/30 | 9/30 | 7/30 |
| Hphob | −0.091 | 12/30 | 4/30 | — |

**Add `Hbond`** — second-best signal and mechanistically independent of RTCNN.

**`VwInt` has no within-ligand signal** (ρ = −0.053, positive in 14/30). The same was true
of Glide's vdW term, so this is consistent across both programs.

**`Eintl` is a route detector, not a pose discriminator.** It looks like the best single
term at 14/30 but it picks the `-n` row **14 times out of 14**, and scores 2/30 on docked
rows alone — the worst of any term. `-n` rows have systematically lower internal strain
(mean 2.80 vs 5.17) because they were minimised from the crystal geometry. This is the
step-4 artefact finally appearing, concentrated in Eintl rather than in Score.

`Natom` and `Nflex` are constant within a ligand and carry no within-ligand information.

### Caveats

With n = 30, GradBoost 12, Score 12 and linear 11 are within noise of each other and
should not be ranked. The 14/30 ceiling is for *global* weights; per-ligand adaptive
weighting could in principle approach the 15/30 single-feature-oracle bound, but on 30
structures that is overfitting.
