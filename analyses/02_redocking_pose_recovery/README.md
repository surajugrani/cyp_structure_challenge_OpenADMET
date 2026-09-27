# 02 — Can Glide reproduce the crystal pose by redocking?

**Question:** Given the correct receptor *and* a grid centred on the true binding site, can Glide regenerate the experimental pose? Analysis 01 tested *scoring given a perfect pose*; this tests *pose generation given a perfect receptor*.

**Status:** Complete.

**Answer: mostly no.** Top-1 success is **4/22 (18 %)** at the standard 2 Å cutoff, and even the best of 10 poses reaches only **6/22 (27 %)**. For roughly three-quarters of these ligands Glide never generates a near-native pose at all.

## Pipeline

| Step | Script | Runs on | Output |
|---|---|---|---|
| 1 | `scripts/01_redock.py` | sbatch | Glide SP `confgen`, 10 poses/ligand |
| 2 | `scripts/02_rmsd.py` | sbatch | `results/pose_rmsd.csv`, `results/pose_rmsd_all_poses.csv` |

```bash
sbatch scripts/run_redock.sh    # redock + RMSD (all 22)
sbatch scripts/run_rmsd.sh      # recompute RMSD only, no re-docking
```

Reuses the prepped receptors and grids from analysis 01 — nothing is re-prepared, so the
whole run takes minutes. The scripts refuse to launch Glide outside SLURM.

## Method notes

**RMSD source.** Glide reports `r_i_glide_rmsd_to_input` per pose (via `CALC_INPUT_RMS`).
Since the ligand file *is* the crystal pose, that is exactly the quantity we want, and it
is symmetry-corrected. Measured **in place** — no superposition — because the receptor is
identical in both cases, so aligning first would hide the displacement being measured.

**Cross-check.** Glide preserves input atom order, so a naive index-wise heavy-atom RMSD
uses the true atom correspondence. Symmetry correction can only *lower* RMSD, so a valid
Glide value must satisfy `glide <= naive`. Across all 212 poses: **0 violations**, max
legitimate symmetry gain 2.83 Å.

> An earlier version cross-checked with `ConformerRmsd`. It was dropped because its
> SMARTS-based renumbering picked degenerate atom mappings for some ligands and
> *inflated* the RMSD — for 4K9T pose 1 it reported 5.85 Å where the naive index-wise
> value and Glide agree on 1.708 Å. The headline numbers were never affected (Glide's
> value was always primary), but it is a good argument for cross-checking the one number
> an analysis actually rests on.

## Results

| Metric | Top-1 | Best-of-10 |
|---|---|---|
| RMSD ≤ 2.0 Å | **4/22 (18 %)** | **6/22 (27 %)** |
| RMSD ≤ 1.0 Å | 2/22 (9 %) | 2/22 (9 %) |
| Median RMSD | 4.85 Å | 3.19 Å |

Published redocking benchmarks on well-behaved targets typically reach 50–80 % top-1
success, so 18 % is poor.

**The successes** (top-1 ≤ 2 Å): 2V0M (0.70 Å), 3UA1 (0.61 Å), 7LXL (1.32 Å), 4K9T (1.71 Å).

### This is a sampling failure, not only a scoring failure

Best-of-10 (27 %) is barely better than top-1 (18 %). Had sampling been working and
ranking been the problem, best-of-10 would be much higher. Instead, a near-native pose is
usually **absent from the pose list entirely**.

The score does have weak, real ability to rank poses *within* a ligand — pooled
Spearman rho = +0.16 (p = 0.02, n = 212 poses); within-ligand median rho = +0.19, and the
sign is correct (better score ↔ lower RMSD). It is just far too weak to matter: the
top-ranked pose is the best available one in only **6/22** cases.

### Likely mechanism: the metal term is inactive

Analysis 01 found that Glide's metal term is identically 0.00 under the chemically
correct (dianionic) heme. **16 of 20 ligands coordinate the heme iron at 1.9–2.4 Å**, yet
the scoring function gives that interaction no credit. Glide therefore has no energetic
reason to place the coordinating nitrogen on the iron, and CYP3A4's active site is a
large, hydrophobic, relatively featureless cavity — so poses drift.

This single observation plausibly explains **both** the failed affinity correlation
(analysis 01) and the poor pose recovery here.

### Nothing else explains the failures

Top-1 RMSD is uncorrelated with the obvious structural suspects:

| | Spearman rho | p |
|---|---|---|
| resolution | +0.10 | 0.67 |
| ligand MW | +0.03 | 0.91 |
| Fe–ligand distance | −0.06 | 0.80 |

**Weak hint:** both multi-copy structures succeeded at rank 1 (2V0M 0.70 Å, 4K9T 1.71 Å).
The second ligand copy occupies pocket volume and restricts placement, which fits the
"under-constrained cavity" explanation — but n = 2, so this is a hint, not a finding.

## Caveats

- **The grid is centred on the crystal ligand**, so this measures pose reproduction given
  the correct site, not site identification. That is the standard redocking setup, but it
  is an optimistic assumption — a prospective run would be harder still.
- **Ligand input geometry** comes from the crystal structure. `confgen` regenerates
  torsions, but bond lengths and angles are crystal-derived, which is a mild advantage
  over a real prospective case. A LigPrep-from-SMILES variant is a planned follow-up.
- Rigid receptor. CYP3A4 is known to reshape its pocket around different ligands, so some
  of this failure is induced fit that rigid docking cannot capture.

## Implication for the challenge

Rigid-receptor Glide docking is not a viable route to CYP3A4 poses on this evidence.
Combined with analysis 01 (the score does not track affinity even for perfect poses),
this supports a **co-folding-first strategy** — matching the PXR challenge result, where
AlphaFold3 co-folding beat docking and Glide post-processing degraded good poses.

**Next test:** re-dock with a Glide metal constraint requiring coordination to the heme
iron. If success rates jump, the mechanism above is confirmed and the fix is directly
useful.
