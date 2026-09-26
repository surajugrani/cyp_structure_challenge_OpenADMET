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
