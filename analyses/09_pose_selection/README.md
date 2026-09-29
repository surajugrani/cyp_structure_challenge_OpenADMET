# 09 — Which scoring function identifies the native pose?

**Question:** For the pose-prediction competition, the metric that matters is whether the
*top-ranked* pose is native-like. Analyses 01 and 08 showed no scoring function predicts
affinity; this asks the different and more relevant question — can any of them **pick the
right pose**?

**Status:** Step 1 (baselines) complete. Step 2 (cross-scoring / selector) not started.

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

## Step 2 (pending user go-ahead)

Test whether any individual score term — Glide's or ICM's — ranks ICM's 10 poses better
than ICM's total score. Must be a **within-ligand** analysis with a structure-level split.
See [[pose-selector-idea]].

## Rerun

```bash
$SCHRODINGER/run python3 scripts/01_baseline_success.py
```

Light: pandas only, no licence, no docking. Outputs `results/rigid9_baseline.csv` and
`results/rigid9_summary.csv`.
