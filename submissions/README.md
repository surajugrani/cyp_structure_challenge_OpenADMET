# Submissions

Contest deliverables for the **OpenADMET CYP3A4 Structure Prediction Track** (track 3 only).

This folder is kept separate from [`analyses/`](../analyses/), which is a trial environment
for working out what does and does not work. Nothing here is exploratory: each subfolder is
a candidate submission, reproducible end to end.

## Task

Predict the CYP3A4 complex for each of the 20 blinded ligands in
[`input/cyp-challenge-TEST-BLINDED_structures.csv`](../input/cyp-challenge-TEST-BLINDED_structures.csv).

- **Submission format:** a single `.zip` of all 20 predicted complexes
- **Primary metric:** LDDT-PLI (OpenStructure), against blinded ground truth
  (20 cryo-EM structures)
- **Secondary:** BiSyRMSD — binding-site superposed, symmetry-corrected pose RMSD
- **Gating:** PoseBusters physical-plausibility checks. Failing OpenStructure validation
  scores LDDT-PLI 0 and a 20 Å BiSyRMSD penalty; failing PoseBusters also scores 0.

## Results

See [`LEADERBOARD.md`](LEADERBOARD.md) for scores returned by the organizers.

| # | LDDT-PLI | BiSyRMSD | LDDT-LP | Coverage | Rank |
|---|---|---|---|---|---|
| 03 | 0.3921 | 4.1179 | 0.7962 | 1 | 5 / 9 |

## Attempts

| # | Method | Settings | Status |
|---|---|---|---|
| 01 | [Boltz-2](01_boltz2_default/) | defaults + heme, MSA server | complete — 20/20 valid, PoseBusters 17/20 |
| 02 | [Boltz-2 + ICM tethered min](02_boltz2_icm_tethered/) | ICM `minimize cartesian`, self-tether | abandoned — corrupted structures |
| 03 | [Boltz-2 + PrepWizard](03_boltz2_prepwizard/) | `-rmsd 0.3` restrained minimisation | **submitted** — 20/20 valid, 20/20 PoseBusters, LDDT-PLI 0.3921, rank 5/9 |
| 04 | [AlphaFold3 default](04_af3_default/) | ColabFold MSA, heme, no templates, 5 samples | superseded by 05 — 13/20 PoseBusters before prep |
| 05 | [AF3 + PrepWizard](05_af3_prepwizard/) | PB-filtered, best mean(PL, haem-L) iPTM | **submitted** — LDDT-PLI 0.4279, BiSyRMSD 4.0653 |
| 06 | [AF3, PL-iPTM only](06_af3_pliptm/) | as 05 but selection on protein-ligand iPTM alone (2 structures differ) | **ready** — 20/20 valid, 20/20 PoseBusters |

## What analyses/ says that matters here

- Docking into crystal receptors recovers the pose in only ~27 % of cases top-1
  (analysis 09), so co-folding is not obviously worse than docking for this target.
- The blinded ligands share **no Murcko scaffold** with the CYP3A4 PDB ligands and have a
  maximum Tanimoto of 0.34 (analysis 07) — template-based approaches are weakly supported.
- Glide's metal term never fires under correct heme chemistry (analysis 03), and docking
  systematically pushes ligands off the heme iron (analysis 04). Whether a co-folding model
  reproduces Fe coordination is worth checking in the output.
