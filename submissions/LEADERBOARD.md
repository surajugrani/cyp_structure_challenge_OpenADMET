# Leaderboard results

Scores returned by the organizers for each submitted attempt. Evaluated on a subset
(~10 of the 20 structures) while the challenge is open.

| # | Submission | Date | LDDT-PLI | BiSyRMSD (Å) | LDDT-LP | Coverage | Rank |
|---|---|---|---|---|---|---|---|
| 03 | [Boltz-2 defaults + heme, PrepWizard-minimised](03_boltz2_prepwizard/) | 2026-09-30 | **0.3921** | **4.1179** | 0.7962 | 1 | 5 / 9 |
| 05 | [AF3 + PrepWizard, PB-filtered, mean(PL,HL) iPTM](05_af3_prepwizard/) | 2026-10-01 | **0.4279** | **4.0653** | 0.6977 | 1 | — |

Metric notes: **LDDT-PLI** is the primary score (protein-ligand interface, higher better).
**BiSyRMSD** is binding-site-superposed symmetry-corrected pose RMSD (lower better).
**LDDT-LP** scores the protein pocket itself. **Coverage** = fraction of structures that
scored at all; 1 means nothing was zeroed by OpenStructure or PoseBusters.

## 03 vs 05

| | 03 Boltz-2 | 05 AF3 | change |
|---|---|---|---|
| LDDT-PLI | 0.3921 | **0.4279** | **+0.0358** |
| BiSyRMSD | 4.1179 | **4.0653** | −0.05 Å |
| LDDT-LP | 0.7962 | 0.6977 | **−0.0985** |
| Coverage | 1 | 1 | — |

**AF3 wins on the metric that counts.** LDDT-PLI is the primary score and it rose 9 %
relative, while BiSyRMSD barely moved — so the gain is in interface contacts, not in
getting the pose fundamentally more right. Both still sit ~4 Å, well outside the ~2 Å
usually called correct.

**LDDT-LP got worse, and that is the surprise.** AF3's protein pocket scores 0.0985 lower
than Boltz's. So AF3 is not simply the better model here: it builds a worse pocket and
still places the ligand better within it. Worth remembering before assuming AF3 dominates
— PXR saw AF3 ahead on both.

That split also means the two methods are not redundant. A hybrid - Boltz protein with
AF3's ligand placement, or consensus between them - is now motivated by data rather than
by hope.

## Reading the first result

**Coverage 1** — every structure scored. The PoseBusters work in submission 03 did its job;
nothing was zeroed on a technicality.

**LDDT-LP 0.7962 vs LDDT-PLI 0.3921** — the protein pocket is predicted well and the
*ligand* is what is wrong. That is the diagnostic that matters: effort spent on protein
quality is wasted, effort on ligand placement is not.

**BiSyRMSD 4.12 Å** — well outside the 2 Å usually called correct, and close to the
~4 Å best-pose values analyses 05/06 saw for structures where docking failed to sample a
near-native pose. Consistent with co-folding putting the ligand in roughly the right
cavity but the wrong orientation.

**Rank 5/9** — mid-field, from an essentially default Boltz-2 run. The gap to the leaders
is therefore unlikely to be closed by cleanup; it needs a different pose strategy.

## Candidate next attempts

Ordered by expected value, not by effort:

1. **`--diffusion_samples 5`** and select among them. Analysis 09 established the selection
   problem directly: oracle 20/30 versus ~12/30 delivered by any single score. More
   candidates plus a selector is the one lever measured to help.
2. **Fe-coordination restraint or filter.** 13/20 blinded ligands carry a classical
   Fe-ligating motif but Boltz coordinated only 11, and BiSyRMSD 4.12 Å suggests
   orientation errors a coordination constraint would penalise.
3. **AlphaFold3** as an independent co-folder, for consensus rather than replacement.
4. **Docking into the Boltz pocket** — but analysis 09 measured only ~27 % top-1 success
   for Glide/ICM on crystal receptors, so this is the weakest of the four.
