# 07 — AF3 with cryo-EM templates (38FY + 38FU)

**Status:** ready. 20/20 structural, 20/20 PoseBusters.
`submission_af3_templates.zip`

Identical to [05](../05_af3_prepwizard/) in every respect — same ColabFold MSA, same
haem/ligand setup, same PrepWizard treatment, same selection rule — **except that the
protein chain now carries two templates**. Templates are the only variable.

## Why these two templates

From the challenge's own cryo-EM preprint (`input/`, see
[analysis 11](../../analyses/11_preprint_insights/)): the **F/G loop** reshapes the CYP3A4
pocket per ligand and is unmodelled in 40 % of PDB entries, and the authors name this as a
likely cause of co-folding failure on this target. The five new cryo-EM structures have it
resolved.

Chosen on **model quality, not ligand similarity** — none of the four cryo-EM ligands
resembles the blinded set (max Tanimoto 0.218, **0/80 pairs ≥ 0.4**), so per-ligand
matching has no basis:

| | FSC 0.5 (Å) | Q score | MolProbity | clash | F/G loop |
|---|---|---|---|---|---|
| **38FY** vardenafil | **3.5** | **0.478** | 1.18 | 3.95 | complete |
| **38FU** apo | **3.5** | 0.416 | **1.00** | **2.24** | complete |
| 38FV ritonavir | 3.7 | 0.423 | 1.02 | 2.42 | complete |
| 38FW ketoconazole | 3.6 | 0.351 | 1.11 | 3.14 | complete |
| 38FX azamulin | 3.6 | 0.393 | 1.25 | 4.82 | **214–217 missing** |

38FX is excluded deliberately: worst clash score and MolProbity, and the only structure
missing part of the F/G loop — the very thing the templates are for. 38FV holds the pocket
open for a 721 Da ligand, against 240–412 Da for the blinded set.

## Three AF3 requirements that must be met, all verified against the AF3 source

1. **Single chain.** `get_polymer_features` raises *"The structure must be filtered to a
   single chain"*; both cryo-EM entries are **C3 trimers** (chains A/B/C). Chain A is
   extracted.
2. **Release date.** It also raises *"The structure must have a release date"*, read from
   `_pdbx_audit_revision_history.revision_date`. A naive extraction drops that category and
   fails with a confusing, unrelated-looking error — so AF3's own
   `Structure.filter()`/`to_mmcif()` is used, and the script re-reads what it wrote to
   confirm both properties survived.
3. **`templateIndices` are not an offset.** They index the template's *full* sequence
   including unmodelled residues. The sequences do not line up arithmetically — 38FU
   chain A is 487 residues and 38FY 485, against a 503-residue query and a nominal Δ3-22
   length of 483 (His tag; 38FY also lacks the leading `MA`). The map comes from a real
   pairwise alignment, keeping only positions where residues match.

Resulting maps: 38FY_A 481 mapped positions (query 22–502), 38FU_A 483 (query 0–502).

## Result: templates barely moved anything

| | 05 no templates | 07 templates |
|---|---|---|
| samples passing PoseBusters | 95/100 | 95/100 |
| ligands with ≥1 passing sample | 20/20 | 20/20 |
| median PL-iPTM over 100 samples | 0.880 | 0.885 |
| mean(PL,HL) of the 20 chosen | 0.769 | 0.774 (up 14/20, down 5/20) |

No change in physical validity and a marginal confidence gain. Expect a small leaderboard
effect at best — and PXR's template runs scored *below* its plain AF3 run (0.4735 vs
0.5055), so the sign is not guaranteed.

## Rerun

```bash
module load alphafold/3.0 && python scripts/01_make_templates.py
python3 scripts/02_make_inputs.py
AF3=$(sbatch --parsable scripts/03_run_af3.sh)
PREP=$(sbatch --dependency=afterok:$AF3 --parsable scripts/run_prep.sh)
sbatch --dependency=afterok:$PREP scripts/run_pb.sh
```

## Leaderboard result: templates hurt

| | 05 no templates | 07 templates | change |
|---|---|---|---|
| LDDT-PLI | 0.4279 | **0.3853** | **−0.0426** |
| BiSyRMSD | 4.0653 | 4.2613 | +0.20 Å |
| LDDT-LP | 0.6977 | 0.6762 | −0.0215 |

All three metrics worse, by ~35x the 05-vs-06 selection difference, so this is a real
effect.

**CORRECTION to an earlier claim in this file:** it previously said PXR's template runs
also scored below plain AF3. That was a misreading of the PXR README — 0.4735 was its best
*Boltz2* submission, not a template run. PXR's actual template results mostly BEAT plain
AF3 (0.5087 template search, 0.5104 +9fzj, 0.5132 +7axe best-iPTM, **0.5285** +9fzj
best-iPTM, against 0.5055 plain), and its best AF3 submission used a template.

**What distinguishes this run is TWO templates, one apo (38FU) and one holo (38FY).**
Every successful PXR template run used exactly one. Apo and holo pockets differ precisely
in the F/G loop region, so supplying both may give AF3 contradictory restraints on the one
feature templates were meant to fix.

A secondary factor may be ligand mismatch: template ligands are 479–721 Da against 240–412
Da for the blinded set, max Tanimoto 0.218, 0/80 pairs >= 0.4.

**Do not conclude "templates don't work" from this — test a single holo template next.**
