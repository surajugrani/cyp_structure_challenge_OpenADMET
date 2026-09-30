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

## Attempts

| # | Method | Settings | Status |
|---|---|---|---|
| 01 | [Boltz-2](01_boltz2_default/) | defaults + heme, MSA server | complete — 20/20 valid, PoseBusters 17/20 |
| 02 | [Boltz-2 + ICM tethered min](02_boltz2_icm_tethered/) | ICM `minimize cartesian`, self-tether | abandoned — corrupted structures |
| 03 | [Boltz-2 + PrepWizard](03_boltz2_prepwizard/) | `-rmsd 0.3` restrained minimisation | **complete — 20/20 valid, 20/20 PoseBusters** |

## What analyses/ says that matters here

- Docking into crystal receptors recovers the pose in only ~27 % of cases top-1
  (analysis 09), so co-folding is not obviously worse than docking for this target.
- The blinded ligands share **no Murcko scaffold** with the CYP3A4 PDB ligands and have a
  maximum Tanimoto of 0.34 (analysis 07) — template-based approaches are weakly supported.
- Glide's metal term never fires under correct heme chemistry (analysis 03), and docking
  systematically pushes ligands off the heme iron (analysis 04). Whether a co-folding model
  reproduces Fe coordination is worth checking in the output.
