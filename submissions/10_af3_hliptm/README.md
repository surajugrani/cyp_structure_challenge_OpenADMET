# Submission 10 — AF3 no-template run, selected on haem-ligand iPTM alone

**Result: 0.4221 LDDT-PLI / 4.0842 Å BiSyRMSD / 0.6953 LDDT-LP.** Slightly the weakest of
the three rules, but the whole three-way spread is 0.0070 LDDT-PLI — see the verdict at
the bottom.

Third point on the selection-metric A/B over the same no-template predictions
(submission 04), completing the set:

| # | selection key | LDDT-PLI |
|---|---|---|
| 05 | mean(`chain_pair_iptm[0][2]`, `chain_pair_iptm[1][2]`) | 0.4279 |
| 06 | `chain_pair_iptm[0][2]` — protein-ligand | 0.4291 |
| **10** | **`chain_pair_iptm[1][2]` — haem-ligand** | **0.4221** |

## Why HL-iPTM is worth a slot

It scores the **coordination geometry** rather than placement in the pocket. 12 of the 20
blinded ligands coordinate Fe through an imidazole / pyridine / triazole nitrogen at
2.06–2.32 Å — classic type II binding — so for most of the set this is the chemically
decisive contact.

It also has by far the widest dynamic range of the three. Best-of-5 across the 20 ligands:

| metric | range |
|---|---|
| global iptm | 0.91 – 0.97 |
| PL-iPTM | 0.81 – 0.95 |
| **HL-iPTM** | **0.51 – 0.78** |

Global iptm is nearly saturated because it averages in the near-perfect protein-haem
interface, which says nothing about the blinded ligand. HL-iPTM is the only one of the
three with room to discriminate.

## What changed

Nothing was recomputed — AF3 inference, PrepWizard on all 100 samples and PoseBusters all
come from submission 04. 95 PB-passing prepped samples, 20/20 ligands covered.
20/20 structural, 20/20 PoseBusters on the final set.

**Differs from 06 (PL-only) in 4 of 20, from 05 (mean) in 2 of 20** — identical to the
counts over all five samples, so the PoseBusters filter did not change which samples win.

| ligand | HL pick | vs 06 PL-only | vs 05 mean |
|---|---|---|---|
| OCNT-2312554 | sample-2 | differs | same |
| OCNT-2395342 | sample-2 | differs | differs |
| OCNT-2395503 | sample-0 | differs | differs |
| OCNT-2535422 | sample-1 | differs | same |

Only 2 of 20 differ from 05 because the mean is half HL by construction.

## Expected noise

[Analysis 13](../../analyses/13_pose_divergence/) measured the no-template pool spread for
these ligands. **OCNT-2395342 has a mean pairwise spread of 3.53 Å and OCNT-2395503
3.83 Å** — both are cases where AF3's five samples disagree strongly, so changing the
metric can move the pose by several Å.

That is exactly the situation that produced the submission 09 regression, where one 3.5 Å
swap cost 0.026 LDDT-PLI on its own. Expect this comparison to move more than 05 → 06 did
(+0.0012), in either direction, and **do not read the sign of the change as a verdict on
HL-iPTM** without checking which ligands fell in the scored subset.

## Result and verdict

| # | selection key | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|---|
| 06 | PL-iPTM alone | **0.4291** | **4.0647** | **0.6985** |
| 05 | mean(PL, HL) | 0.4279 | 4.0653 | 0.6977 |
| 10 | HL-iPTM alone | 0.4221 | 4.0842 | 0.6953 |

HL-only comes last on all three metrics, so the wider dynamic range did not translate into
a better choice. But the **total spread across all three rules is 0.0070 LDDT-PLI** —
against submission 09, where one ligand's 3.5 Å pose shift moved the score by 0.0258 on
its own. The metric choice is worth roughly a quarter of one unlucky pose.

The ordering is monotone in PL weight (PL > mean > HL), a weak hint that protein-ligand
confidence carries more information than the haem-ligand term. Not enough to build on.

Worth noting against the prediction made above: the 2 ligands that differ from 05,
OCNT-2395342 and OCNT-2395503, are precisely the high-spread pools (3.53 Å and 3.83 Å)
that analysis 13 flagged as risky — yet the score moved only 0.0058. Either neither fell
in the ~10-complex scored subset, or both swaps were benign. **The variance lives in which
ligands get scored, not in the rule.**

**Selection-metric variants are now exhausted.** Three rules over identical predictions,
0.007 apart, all inside single-structure noise. Further A/Bs of this kind are not worth a
submission slot.
