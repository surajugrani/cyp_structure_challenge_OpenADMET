# Project Log

Chronological record of what was done, why, and what came out of it. Newest entries at the bottom. Each analysis has its own folder under `analyses/` with a detailed README.

---

## 2026-09-26: Repo setup

- Created the repo with the challenge inputs (`input/`): CYP3A4 sequence, 20 blinded test ligands, and an RCSB custom report of 122 CYP3A4 PDB entries.
- Added a `.gitignore` to keep large model outputs and docking intermediates out of git.

## 2026-09-26: Analysis 01, Glide score vs IC50 (steps 1–2)

**Folder:** `analyses/01_glide_score_vs_ic50/`

**Goal:** Check whether Glide scores of crystal poses correlate with measured IC50, as a baseline for trusting Glide during pose ranking.

**Done:**
- Fetched affinities for all 122 entries from the RCSB Data API: 70 values (58 IC50, 10 Kd, 2 Ki), all from BindingDB, across 26 entries.
- Curated to one pIC50 and one structure per ligand: **22 ligands, 20 in the main set**.

**Decisions:** IC50 only; median pIC50 per ligand; best-resolution structure per ligand; Ir complexes removed; qualified values dropped.

**Open issues for step 3:**
- 2V0M has 2 ketoconazole copies per chain, and 4K9T also has 2 copies. We need to decide whether to score one copy (with the other kept as part of the receptor) or to remove the second one.
- 3 stereoisomers share one 310 nM value; only one is kept.

**Next:** Glide steps 3–5 (PrepWizard → grid → `inplace` and `mininplace` scoring → correlation), submitted as an sbatch job.

## 2026-09-26: Analysis 01, Glide pipeline implemented (steps 3-4)

**Done:** Wrote structure prep (download, chain/ligand selection, PrepWizard) and Glide
grid + scoring (`inplace` and `mininplace`). Smoke-tested end-to-end on 8SPD via sbatch
(job 47770421, COMPLETED): docking score `inplace` -7.91, `mininplace` -6.86.

**Decisions:**
- Option A for multi-copy structures (2V0M, 4K9T): score the copy closest to the heme Fe,
  keep the other copy as part of the receptor.
- Deprotonate heme propionates after PrepWizard (heme charge +1 -> -1).

**Three schrodinger/2025-3 incompatibilities found vs the PXR scripts** (see analysis
README): `-f OPLS_2005` not `-f 2005`; grids are written as `<job>-gridgen.zip`; Glide
writes a results CSV, not a `.rept` file.

**Note:** steps 3-4 must be submitted with sbatch, never run on the login node.

**Observation to watch in the full run:** for 8SPD, `mininplace` scored *worse* than
`inplace` (-6.86 vs -7.91) with a much larger evdw term. If this is systematic it may
reflect how the minimizer treats the Fe coordination.

**Next:** submit the full 22-structure run, then write step 5 (correlation).

## 2026-09-26: Analysis 01 complete (first term pair) — negative result

**Full run:** job 47770497, COMPLETED in 1h16m, all 22 structures prepped and scored.

**Heme inconsistency found and fixed.** PrepWizard was not deterministic: 21 structures
came out as the porphyrin dianion with deprotonated propionates (heme charge −2), but
8SPD kept an extra pyrrole N−H (charge −1). Since the heme sits directly under the
ligand, that makes scores non-comparable. `normalize_heme()` now *enforces* the −2 state
rather than trusting PrepWizard; re-normalized and re-scored in job 47772332, all 22
verified at −2.

**Bug fixed:** step 4 with `--only` overwrote the whole results CSV with just the re-run
subset; it now merges.

**Result (job 47772502):** no usable correlation between Glide score and pIC50.
`inplace` rho = −0.265 (p = 0.26, 95 % CI [−0.72, +0.28]); `mininplace` rho = −0.093.
All four sensitivity checks agree. The score largely tracks molecular size
(MW vs score rho = −0.47) while size does not track potency (MW vs pIC50 rho = 0.12).

**Implications for the challenge:**
- Glide score is not a usable affinity proxy for CYP3A4 even with a perfect pose — an
  upper bound, since docking a predicted pose only adds error.
- `mininplace` degrades poses relative to `inplace`, echoing the PXR result that Glide
  post-processing worsened good AF3 poses.

**Stage 2 (jobs 47772548, 47772583):** `ecoul` shows no signal either (inplace
rho = +0.15, wrong sign; mininplace rho = −0.08). Pearson added alongside Spearman for
all terms.

**Glide's metal term is identically 0.00 for all 20 structures** under the chemically
correct dianionic porphyrin, so it was dropped from the correlation table. Three 8SPD
heme states pin down the cause: propionates COOH + one pyrrole N−H (heme +1) -> metal
−2.30; propionates deprotonated + one pyrrole N−H (heme −1) -> metal −2.30 (unchanged, so
the carboxylates are irrelevant); both pyrrole N⁻ (heme −2) -> metal 0.00. Not a net-charge
effect — a 2-unit charge change did nothing while a 1-unit change zeroed it. It is the
iron's first coordination shell specifically. Decision: keep the correct chemistry and accept that the
metal term is inactive. This is the most likely explanation for the failed correlation —
16/20 ligands coordinate the Fe at 1.9–2.4 Å and the score gives that no credit.

**Next:** analysis 01 is complete. Possible follow-ups (not yet run): ligand efficiency
(score / heavy atoms) and MW-controlled partial correlation, to test whether removing the
size bias recovers any signal.

## 2026-09-26: Analysis 02 — redocking pose recovery

**Folder:** `analyses/02_redocking_pose_recovery/`

**Goal:** Given the correct receptor and a grid on the true site, can Glide regenerate the
crystal pose? (01 tested scoring given a perfect pose; this tests pose generation.)

**Run:** jobs 47772605 (smoke test, 8SPD), 47772655 (full, 22 x 10 poses), 47772828
(RMSD recompute). Reused analysis 01's receptors and grids, so the full run took minutes.

**Result:** top-1 RMSD <= 2 A in **4/22 (18 %)**; best-of-10 **6/22 (27 %)**; median top-1
4.85 A. Typical benchmarks reach 50-80 %, so this is poor.

**Key point — this is a sampling failure, not only scoring.** Best-of-10 (27 %) is barely
above top-1 (18 %), so a near-native pose is usually absent from the pose list entirely.
The score does rank poses within a ligand slightly better than chance (pooled rho = +0.16,
p = 0.02, n = 212; correct sign) but the top pose is the best available in only 6/22 cases.
This revised an earlier expectation based on 8SPD alone, where sampling had succeeded
(best 1.71 A at rank 9) and only ranking failed — 8SPD was not representative.

**Mechanism linking 01 and 02:** the Glide metal term is identically 0.00 under the correct
heme, yet 16/20 ligands coordinate the Fe at 1.9-2.4 A. With no reward for that
interaction and a large hydrophobic cavity, poses drift. Explains both results.

**Method note:** the `ConformerRmsd` cross-check was dropped — its SMARTS renumbering
inflated RMSD for some ligands (4K9T pose 1: 5.85 A vs the correct 1.708 A). Replaced with
a naive index-wise RMSD plus the invariant `glide <= naive`; 0 violations across 212 poses.
Headline numbers were never affected since Glide's value was always primary.

**Safeguard added:** `require_compute_node()` in the prep/grid/dock scripts now refuses to
launch PrepWizard or Glide outside SLURM, after two accidental login-node runs.

**Follow-up measurement:** Glide displaces ligands off the Fe — median crystal distance
2.33 A vs 3.24 A for top-1 poses; 20/22 move outward, median +0.73 A, Wilcoxon p = 0.0003.
The two exceptions are the two non-coordinating ligands.

**CORRECTION (2026-09-27).** I first explained this as Glide penalising the crystal geometry
as a steric clash. That was wrong. (a) A dative N->Fe bond at 2.0-2.4 A is normal
coordination chemistry, not an error in 22 crystal structures. (b) inplace_evdw is
favourable for all 22 (-14 to -65). There IS localised strain at short Fe contacts (vdW gain
on minimization vs Fe distance: rho = -0.86, p < 0.0001; 21-31 kcal at 1.94-2.01 A vs 1-4
kcal for non-coordinating), but GlideScore weights evdw at only ~0.046 (fitted by least
squares over 212 poses; lipo 1.19, rewards 1.05, ecoul 0.19, esite 2.51). So that strain is
worth ~1 GlideScore unit while a -2.30 metal term enters at weight ~1.0. The correct
mechanism is the *absence of a reward*, not the presence of a penalty — and the metal term
is decisive at this scale, not marginal. The positive-scoring crystal poses (4K9V, 4D6Z)
come from the `rewards` term, not vdW.

**Next:** analysis 03 — try to reactivate the metal term without breaking the chemistry.
Primary lever: set Fe to the physiological ferric Fe(III) state (PrepWizard leaves it +2).

## 2026-09-27: Analysis 03 — heme representation and the metal term

**Folder:** `analyses/03_heme_representation/`

**Goal:** Can Glide's metal term be reactivated without misrepresenting the chemistry?
Gate test on 3 structures (8SPD, 4K9X, 2V0M), 4 heme variants, grids rebuilt per variant,
all 13 Glide terms recorded. Jobs 47778194 and 47778201 failed on my own path bugs
(SDF carries no chain/resnum; Glide resolves LIGANDFILE relative to cwd); 47778301 ran.

**Results — identical across all three structures:**

| Variant | Fe | Porphyrin | Heme | Metal |
|---|---|---|---|---|
| fe2_dianion (01/02 baseline) | +2 | dianion | -2 | 0.00 |
| fe3_dianion (chemically correct) | +3 | dianion | -1 | **0.00** |
| fe2_monoanion | +2 | monoanion | -1 | 0.00 |
| fe3_monoanion | +3 | monoanion | 0 | **-2.30** |

1. **Fe(III) alone does not work.** The physiological ferric resting state with the correct
   dianionic porphyrin still gives metal = 0.00.
2. **Not a net-charge effect.** fe3_dianion and fe2_monoanion are both net -1 and give 0.00,
   while 8SPD "state B" (also -1, but with a pyrrole N-H) gave -2.30. Glide's perception
   depends on something more specific; left unresolved rather than guessed at.
3. **The metal term is a fixed reward, not distance-graded** — exactly -2.30 in all five
   activated cases seen so far, across Fe distances 1.94-2.10 A and different ligands.
   So even an active metal term would NOT rescue analysis 01: a constant applied to every
   coordinating ligand adds no discrimination. It should still help pose recovery.
4. **Heme representation shifts scores more than the metal term does.** 8SPD ecoul swings
   -1.00 -> +7.02 -> -7.12 across variants (14 raw units ~ 2.7 GlideScore units at weight
   0.19). fe3_monoanion's better total comes from metal AND electrostatics.

**Decision:** keep the chemically correct heme; do not adopt the monoanion. Leave analyses
01/02 intact as the Fe(II) baseline — the sensitivity is a more useful result than any one
set of numbers.

**Next:** Glide `metal_cons` at grid generation — imposes coordination geometrically
without breaking the chemistry. Prediction recorded in advance: improves pose RMSD while
leaving the metal term at 0.00 (a geometric fix, not a scoring fix).

---

## Planned work

**CYP inhibition dataset (`input/cyp-challenge_cyp-inhibition-data.csv`).** 4,905 compounds
with pIC50 for CYP1A2/2C9/2D6/3A4 plus confidence bounds and std. This belongs to a
different track of the same challenge and will be used later.

Roughly **15 of the 20 blinded structure-track compounds also appear in this dataset**.
There are no experimental structures for them, so unlike analysis 01 (which scores
crystal poses in place) that work will require **docking and then scoring** — pose
generation first, affinity correlation second. Analysis 01 sets the baseline for how far
a Glide score can be trusted before adding pose-prediction error on top.
