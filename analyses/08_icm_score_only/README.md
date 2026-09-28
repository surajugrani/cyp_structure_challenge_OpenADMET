# 08 — ICM scoring of crystal poses vs IC50

**Question:** Analysis 01 found Glide's score of the crystal pose does not track pIC50
(rho = −0.27). Do ICM's scores — in particular the neural-network **RTCNN** rescoring — do
better on the same poses?

**Status:** Complete. **Answer: no.**

## Method

`_dockScan -n` ("generate and score near-native conformation as provided in input") scores
the supplied crystal pose through ICM's normal VLS pipeline, returning `Score`, `RTCNN` and
`mfScore` in the hitlist. Reuses analysis 06's ICM-prepared receptors.

**Verified it scores the given pose rather than re-docking:** each output contains exactly
one pose, 0.28–1.01 Å from the crystal coordinates (a re-dock gives 3–7 Å).

**It minimises locally first**, so the like-for-like Glide comparison is `mininplace`
(which also minimises), **not** `inplace`.

### Why not `scanScoreExternal`

Two reasons, both worth recording:

1. **It is broken in icm/3.9-4.** Its own macro calls `minimize cartesian 1000 "vw,..."`,
   a syntax this binary rejects (it accepts `minimize cartesian|mmff [type] [charge]`).
   Fails identically via `-s script`, via stdin, and with or without `convertObject`.
   (Note: `convertObject a_2.` with *defaults* does fix the separate "not a true ICM-object"
   error — the all-`no` PXR-style call does not, since it converts nothing.)
2. Even working, ICM's docs state its scores are not comparable to VLS docking scores —
   VLS uses the lowest-energy MC conformer as the ligand free state, the interactive score
   only minimises locally.

## Results

| Score | All 20 | Rigid only (RotB ≤ 8, n = 7) |
|---|---|---|
| ICM Score | −0.093 (p = 0.70) | −0.464 (p = 0.29) |
| **ICM RTCNN** | **+0.086** (p = 0.72) | **−0.857 (p = 0.014)** |
| ICM mfScore | +0.081 (p = 0.73) | −0.429 (p = 0.34) |
| Glide `inplace` | −0.265 (p = 0.26) | −0.643 (p = 0.12) |
| Glide `mininplace` | −0.093 (p = 0.70) | −0.643 (p = 0.12) |

**Like-for-like** (both minimise the pose): ICM Score −0.093 vs Glide `mininplace` −0.093.
Identical, and both null.

## The RTCNN rigid-subset result is almost certainly spurious

It is the only nominally significant correlation in this entire project, and it should not
be believed:

- **It flips sign between subsets** — +0.086 across all 20, −0.857 on the 7 rigid ligands.
  A real scoring function does not reverse direction depending on which ligands are
  included. That is the signature of small-sample noise.
- **Multiple testing.** Roughly 30 score-vs-pIC50 correlations have now been computed over
  these same ligands (Glide inplace/mininplace/evdw/ecoul; ICM Score/RTCNN/mfScore; on
  crystal, top and best poses; across two subsets). One or two crossing p < 0.05 by chance
  is expected. Bonferroni on p = 0.014 across 30 tests gives **p ≈ 0.42**.
- **The subset was chosen post-hoc**, after the flexibility relationship in analysis 07.

## Conclusion

**No scoring function tested — Glide or ICM, physics-based or neural — correlates reliably
with pIC50 on CYP3A4 crystal poses.** Analysis 01's negative result is not a Glide-specific
artefact; it reproduces in a second program with a different scoring model.

This strengthens rather than weakens the case for the planned pose selector: the right
target is **pose quality (RMSD)**, where there are ~1000 labelled examples and a proper
held-out split is possible, not **affinity**, where n = 20 with heterogeneous BindingDB
values cannot support any conclusion. See [[pose-selector-idea]].
