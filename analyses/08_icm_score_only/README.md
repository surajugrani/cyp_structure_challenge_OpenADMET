# 08 — ICM scoring of crystal poses vs IC50

**Question:** Analysis 01 found Glide's score of the crystal pose does not track pIC50
(rho = −0.27). Do ICM's scores — in particular the neural-network **RTCNN** rescoring — do
better on the same poses?

**Status:** Complete. **Answer: no.**

## Method

`_dockScan -n` ("generate and score near-native conformation as provided in input") scores
the supplied crystal pose through ICM's normal VLS pipeline, returning `Score`, `RTCNN` and
`mfScore` in the hitlist. Reuses analysis 06's ICM-prepared receptors.

**It scores a minimised pose, not the deposited one.** Each output contains exactly one
pose, but that pose is **0.28–1.42 Å from the crystal coordinates** (n = 22, median
0.86 Å). It is near-native — a re-dock gives 3–7 Å — but it is not the crystal pose. The
like-for-like Glide comparison is therefore `mininplace` (which also minimises), **not**
`inplace`.

### `scanScoreExternal` — the second route

**Earlier versions of this README claimed `scanScoreExternal` was broken in icm/3.9-4.
That was wrong.** The failures came from passing the receptor and ligand as two separate
objects. The working pattern (from Molsoft usage) merges the ligand into the receptor
object and passes *molecule selections*:

```
read pdb "${PDB}_rec.pdb" name=rec
read mol "${PDB}_crystal_v2000.sdf" name=lig
move a_2. a_1.          # merge ligand into the receptor object
convertObject a_1.      # default options
receptor = a_2. - a_2.m
scanScoreExternal receptor a_2.m
```

See `scripts/03_sse_rigid7.sh`. It is run non-interactively via a heredoc to `icm64`.

**It also minimises** — its macro runs a Cartesian minimisation before scoring — so it too
is not scoring the deposited coordinates. It writes no pose, so unlike the `-n` route the
displacement cannot be measured. ICM's docs additionally note that interactive scores are
not on the same scale as VLS docking scores: VLS uses the lowest-energy MC conformer as
the ligand free state, the interactive score only minimises locally.

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

### `scanScoreExternal` on the same 7 rigid ligands

| comp_id | pdb_id | RotB | Score | RTCNN | pIC50 |
|---|---|---|---|---|---|
| KLN | 2V0M | 7 | −7.657 | −44.367 | 7.426 |
| CL6 | 8SPD | 4 | −1.018 | −16.463 | 7.237 |
| YNV | 7LXL | 4 | +29.594 | −23.483 | 6.097 |
| 08Y | 3UA1 | 5 | +11.918 | −47.947 | 5.523 |
| MYT | 6MA6 | 3 | +8.706 | −14.477 | 5.415 |
| TPF | 6MA7 | 5 | +42.098 | −2.538 | 4.620 |
| PK9 | 4D6Z | 8 | +37.830 | −11.656 | 4.125 |

Spearman vs pIC50: **RTCNN −0.714 (p = 0.071)**, **Score −0.821 (p = 0.023)**.
Plot: `results/sse_rtcnn_rigid7.png`.

Compared with the `-n` route on the same ligands (RTCNN −0.857, Score −0.464), the two
methods **disagree about which score is the better predictor** — and the sign of that
disagreement reverses. This is not two routes scoring one pose: each minimises the ligand
independently, so they are scoring different (though both near-native) geometries. With
n = 7 the spread across routes is itself evidence that none of these estimates is stable.

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
