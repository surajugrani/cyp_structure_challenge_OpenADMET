# Submission 10 — AF3 no-template run, selected on haem-ligand iPTM alone

Third point on the selection-metric A/B over the same no-template predictions
(submission 04), completing the set:

| # | selection key | LDDT-PLI |
|---|---|---|
| 05 | mean(`chain_pair_iptm[0][2]`, `chain_pair_iptm[1][2]`) | 0.4279 |
| 06 | `chain_pair_iptm[0][2]` — protein-ligand | 0.4291 |
| **10** | **`chain_pair_iptm[1][2]` — haem-ligand** | pending |

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
