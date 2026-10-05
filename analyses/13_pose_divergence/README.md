# Analysis 13 — how far apart are AF3's own samples?

**Question:** submission 09 is submission 08 with one line changed, 3 of 20 structures
differ, and LDDT-PLI fell 0.0258 while BiSyRMSD rose 0.17 Å. Is the selection rule worse,
or did one structure go badly wrong?

**Answer: one structure.** And the measurement that fell out of it is more useful than the
answer — AF3's five samples per ligand are genuinely different binding modes, not
variations on one.

## Method

For each pair of PoseBusters-passing, PrepWizard-minimised samples of the same ligand:

1. define the pocket as the protein residues with a heavy atom within 10 Å of the first
   sample's ligand (this mirrors how BiSyRMSD is scored — binding-site superposition, not
   global);
2. Kabsch-superpose the second sample's pocket Cα onto the first's;
3. apply that transform to its ligand and compute the symmetry-corrected RMSD
   (`rdMolAlign.CalcRMS` after `assign_bond_orders_from_template`, with a name-matched
   fallback).

Scripts: `scripts/01_pair_rmsd.py` (38FW-template run, submission 08) and
`scripts/02_pool_rmsd_notemplate.py` (no-template run, submission 04). Both are pure CPU
over files already on disk — 5 and 7 seconds respectively.

## Result 1 — the 09 regression is one ligand

| ligand | 08 → 09 pose shift | pool mean | pool max |
|---|---|---|---|
| OCNT-2313300 | 0.50 Å | 1.83 | 2.50 |
| OCNT-2395328 | 0.66 Å | 1.52 | 2.00 |
| **OCNT-2395415** | **3.51 Å** | 3.53 | 4.86 |

Two swaps barely moved the ligand. The third reoriented it by 3.5 Å.

The arithmetic closes: the leaderboard scores ~10 complexes, so one bad structure raises
mean BiSyRMSD by its own error ÷ 10. Observed +0.173 Å implies that structure got ~1.7 Å
worse — consistent with a 3.5 Å reorientation landing further from the answer. Likewise
−0.0258 × 10 = −0.26 LDDT-PLI concentrated in one complex. LDDT-LP moved +0.0032, i.e.
not at all, confirming the protein is identical and only the ligand moved.

**So 09 does not show that PL-iPTM is a worse rule than mean(PL,HL).** It shows that with
~10 scored structures, a single unlucky swap dominates the score. The same logic retro-
actively deflates 05 → 06: that A/B moved 2 structures for +0.0012 LDDT-PLI, which against
this noise floor is nothing.

## Result 2 — the pools are wide, and this is the important number

| run | pairs | median | mean | max |
|---|---|---|---|---|
| 08, 38FW template | 117 | 1.81 Å | 2.01 Å | 6.42 Å |
| 04, no template | 181 | 1.87 Å | 2.25 Å | 6.85 Å |

Both runs, same picture. Per-ligand mean spread ranges from **0.21 Å** (OCNT-2395464) to
**5.48 Å** (OCNT-2313429) in the no-template run.

Set that against the confidence metrics measured over the same samples: within one
ligand's pool, best-to-worst PL-iPTM is typically **0.01–0.03**. The scoring signal is flat
along an axis the poses vary by ~2 Å on. That is a cleaner statement of the same wall
analyses 09 and 10 hit — no score predicts RMSD (oracle 20/30 vs ~12/30), no AF3 metric
predicts PoseBusters validity (14/20 vs 11.4/20 random) — and it explains *why* every
selection-metric A/B has returned noise.

## What follows

**There is real headroom.** A pool whose members differ by ~2 Å contains both better and
worse poses than the one we ship; the oracle is meaningfully above what we achieve. The
problem is entirely that we cannot tell them apart.

**Ensemble spread is a new, orthogonal signal.** It is not a score of a pose, it is a
measure of how much AF3 agrees with itself, and it is not derivable from iPTM. Tight pools
(OCNT-2395464 at 0.21 Å, OCNT-2313038 at 0.49 Å) are plausibly the reliable predictions;
wide ones (OCNT-2313429 at 5.48 Å, OCNT-2395328 at 4.55 Å) are cases where AF3 is simply
undecided and any single pick is close to a coin flip.

**Caveat on consensus selection.** Picking the ensemble centroid needs the ensemble to
have one. For the four no-template ligands with mean spread above 3.4 Å — OCNT-2313429,
OCNT-2395328, OCNT-2395342, OCNT-2395415 — the centroid is weakly defined and the method
has little to offer. Consensus is worth testing, but it is not a general fix.

## Files

| file | contents |
|---|---|
| `results/swap_rmsd.csv` | per ligand: 08's pick, 09's pick, the shift between them, pool spread |
| `results/pool_pairwise_rmsd.csv` | all 117 within-ligand pairs, template run |
| `results/notemplate_pool_summary.csv` | per ligand: n, mean / median / min / max spread |
| `results/notemplate_pairwise_rmsd.csv` | all 181 within-ligand pairs, no-template run |
