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
correct dianionic porphyrin, so it was dropped from the correlation table. It was −2.30
for 8SPD under PrepWizard's original (chemically wrong) mono-protonated porphyrin, so the
term is sensitive to that choice. Decision: keep the correct chemistry and accept that the
metal term is inactive. This is the most likely explanation for the failed correlation —
16/20 ligands coordinate the Fe at 1.9–2.4 Å and the score gives that no credit.

**Next:** analysis 01 is complete. Possible follow-ups (not yet run): ligand efficiency
(score / heavy atoms) and MW-controlled partial correlation, to test whether removing the
size bias recovers any signal.

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
