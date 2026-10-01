# 05 — AlphaFold3 + PrepWizard, selected on mean(protein–ligand, haem–ligand) iPTM

**Status:** ready. **20/20 structural checks, 20/20 PoseBusters.**
`submission_af3_prepwizard.zip`

## Pipeline

| step | where | cost |
|---|---|---|
| AF3 co-folding, 5 samples/ligand, ColabFold MSA, no templates | [04](../04_af3_default/) | 1 h |
| PrepWizard on all 100 samples | [04/glide_rescoring](../04_af3_default/glide_rescoring/) | 9 h |
| PoseBusters on all 100 prepped | same | 45 min |
| select + strip H + validate | here | 5 min |

**MSA searched zero times.** All 20 ligands share one CYP3A4 sequence, so the ColabFold
MSA that Boltz already fetched in submission 01 was reused via `unpairedMsaPath`
(14,301 sequences). One gotcha: MMseqs2 terminates its a3m with a NUL byte, which Boltz
tolerates and AF3 rejects with "Unknown residues in MSA" — it must be stripped.

## Selection

Among each ligand's **PoseBusters-passing** prepped samples, the highest
`mean(chain_pair_iptm[0][2], chain_pair_iptm[1][2])` — protein–ligand and haem–ligand
interface confidence averaged. The first scores placement in the pocket, the second the
coordination geometry, which is what matters for a type-II CYP binder. Global `iptm` is
not used: it averages in the near-perfect protein–haem interface (~0.97) and says little
about the blinded ligand.

Chosen samples average 0.769. **Haem–ligand confidence is consistently much lower than
protein–ligand (0.51–0.78 vs 0.81–0.95)** — AF3 is far less sure about coordination
geometry than about pocket placement. Including it changed the pick for only 2 of 20
versus `pl_iptm` alone.

## Why PoseBusters is the filter and not a confidence metric

[Analysis 10](../../analyses/10_af3_selection_vs_posebusters/) tested every AF3 metric
against PoseBusters outcome on all 100 samples. None predicts it — the best gets 14/20
against an 11.4/20 random baseline, and `has_clash` is **identically 0.0 across all 100
samples** while 43 % fail a distance check.

Glide was tested too ([04/glide_rescoring](../04_af3_default/glide_rescoring/), 97/100
scored): selecting on any Glide term gives **10–12/20, at or below random**. Median gscore
is −5.68 for PB-passing poses and −6.26 for failing ones (MWU p = 0.042) — Glide mildly
*prefers* the poses PoseBusters rejects, because it rewards the close contacts PoseBusters
flags. `docking_score` is identical to `gscore` here (all 97 samples, max difference
0.0000), as expected with no conformational search.

## What PrepWizard did

| | before | after |
|---|---|---|
| samples passing PoseBusters | 57/100 | **95/100** (fixed 38, broke 0) |
| ligands with ≥1 passing sample | 19/20 | **20/20** |

It rescued OCNT-2535519, whose 5 raw samples all failed, by moving the ligand–Fe contact
from 2.03 Å to 2.10–2.15 Å.

**Fe distance is not a clean threshold.** Among the 5 samples still failing, two have Fe at
3.66 and 3.82 Å (ordinary protein clashes), and the minimum Fe distance among *passing*
structures is 2.02 Å. It is the specific contact geometry, not the Fe distance alone.

**Covalent Fe–ligand bonds matter.** 49/100 samples place a ligand atom ≤2.2 Å from the
iron and PrepWizard (`-rehtreat`) bonds them — 58/100 had a bond to break. Splitting the
complex across such a bond would have handed Glide a fragment with a dangling valence, so
the rescoring script deletes the Fe–ligand bond before extracting.

Hydrogens are stripped from the submitted files: PrepWizard adds them (HEM 43 → 73 atoms)
while submission 03 and the organizers' own samples are heavy-atom only, and LDDT-PLI and
BiSyRMSD are heavy-atom metrics.

## Rerun

```bash
python3 scripts/01_select_and_package.py
sbatch scripts/run_validate.sh      # add --skip-pb to the python call to skip PoseBusters
```
