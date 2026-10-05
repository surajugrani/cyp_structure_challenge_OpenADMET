# Submission 08 — AF3 with a single holo cryo-EM template (38FW)

**Result: 0.4192 LDDT-PLI / 4.0543 Å BiSyRMSD / 0.6894 LDDT-LP.**
Better than the two-template run (07) on every metric, still slightly below no template (05).

## Why this run exists

Submission 07 supplied AF3 with two cryo-EM templates — **38FY** (holo, vardenafil) and
**38FU** (apo) — and every metric fell. The hypothesis was that an apo/holo pair gives AF3
contradictory restraints on exactly the F/G-loop region templates were meant to fix, and
that templates per se are not the problem: on PXR, single templates mostly *helped*, and
PXR's best AF3 submission used one. This run isolates that variable — one template, holo.

## Why 38FW and not 38FY

Chosen on **binding mode**, not map resolution. 12 of the 20 blinded ligands coordinate the
haem iron through an imidazole / pyridine / triazole nitrogen at 2.06–2.32 Å — classic
type II binding, which is exactly what **ketoconazole (38FW)** does. Vardenafil (38FY) is
described in the cryo-EM preprint as neither type I nor type II, so its pocket geometry is
the less representative of the two for this ligand set.

## Pipeline

Identical to submission 05 apart from the template block, so the comparison is clean.

| step | script | notes |
|---|---|---|
| 1 | `scripts/01_make_templates.py` | extracts chain A of 38FW with AF3's own `Structure.filter()` / `to_mmcif()`, which preserves `_pdbx_audit_revision_history` (AF3 rejects mmCIFs without a release date); builds the query→template index map from a Bio.pairwise2 alignment |
| 2 | `scripts/02_make_inputs.py` | AF3 JSON per ligand; MSA reused from the Boltz-2 MMseqs2 run (uniref + bfd a3m, trailing `\x00` stripped) |
| 3 | job 47869043 | AF3 inference, 20 ligands × 5 samples, 59:50 |
| 4 | job 47869044 | PrepWizard on all 100 samples, 6:56:38 |
| 5 | job 47869045 | PoseBusters + selection + packaging, 17:39 |

Template build: `38FW_A.cif`, 1 chain, 487 residues, 483 mapped positions over query 0–502.

## Numbers from the run

- **80/100** samples survived PrepWizard (05 managed 95/100).
- **77/80** passed PoseBusters.
- All 20 ligands had at least one passing sample, so **20/20 structural, 20/20 PoseBusters** —
  nothing was zeroed.
- Selection: highest `mean(chain_pair_iptm[0][2], chain_pair_iptm[1][2])` among passing
  samples, the same rule as 05. Mean selected PL-iPTM 0.892.
- **16 of 20** chosen samples differ from the ones 05 chose, so the template genuinely moved
  the predictions rather than reproducing them.

## Verdict

| | 05 no template | 07 two templates | 08 one holo template |
|---|---|---|---|
| LDDT-PLI | **0.4279** | 0.3853 | 0.4192 |
| BiSyRMSD | 4.0653 | 4.2613 | **4.0543** |
| LDDT-LP | **0.6977** | 0.6762 | 0.6894 |

Going from two templates to one recovers **+0.0339 LDDT-PLI** and **0.21 Å** of BiSyRMSD —
roughly 80 % of submission 07's regression. The apo/holo conflict was real.

What it does not do is beat having no template at all: 08 is 0.0087 LDDT-PLI and 0.0083
LDDT-LP below 05. Its one win is BiSyRMSD at 4.0543 Å, the best of any attempt so far, so
the template pulls the pose marginally closer while slightly costing interface contacts.
Both gaps are small enough to be run-to-run noise.

**Templates are a dead end on CYP3A4, unlike PXR, and no further variants are worth
running.** The most likely reason is chemical distance: the four ligand-bound cryo-EM
structures carry 479–721 Da ligands against 240–412 Da for the blinded set, max Tanimoto
0.218, with 0 of 80 pairs reaching 0.4. PXR's useful templates were close analogues of its
targets; these are not.

Best attempt overall remains **06** (no template, PL-iPTM selection, 0.4291).
