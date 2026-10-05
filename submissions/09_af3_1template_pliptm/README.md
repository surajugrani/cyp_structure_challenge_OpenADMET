# Submission 09 — as 08, but selected on PL-iPTM alone

**Result: 0.3934 LDDT-PLI / 4.2269 Å BiSyRMSD / 0.6926 LDDT-LP.**

Submission 08 with exactly one line changed: among a ligand's PoseBusters-passing,
PrepWizard-minimised samples, take the highest `chain_pair_iptm[0][2]` rather than its
mean with the haem-ligand term. The same A/B that 05 → 06 ran on the no-template
predictions.

Nothing was recomputed. AF3 inference (job 47869043), PrepWizard on all 100 samples
(job 47869044) and PoseBusters on the 80 that survived prep (job 47869045) all live in
submission 08; `scripts/01_select_and_package.py` re-reads `08/results/sample_pb.csv`,
re-selects, strips hydrogens, re-runs the organizers' structural checks plus PoseBusters
on the final 20, and zips. 20/20 structural, 20/20 PoseBusters.

## What changed

**3 of 20** structures: OCNT-2313300 (→ sample-4), OCNT-2395328 (→ sample-0),
OCNT-2395415 (→ sample-3). In two of the three the competing samples were tied in PL-iPTM
to two decimals, so the switch was decided by the absence of the haem term rather than by
any real PL preference.

| | 08 mean(PL,HL) | 09 PL-iPTM only | change |
|---|---|---|---|
| LDDT-PLI | **0.4192** | 0.3934 | **−0.0258** |
| BiSyRMSD | **4.0543** | 4.2269 | +0.17 Å |
| LDDT-LP | 0.6894 | 0.6926 | +0.0032 |

## Why this is NOT evidence that PL-only is a worse rule

[Analysis 13](../../analyses/13_pose_divergence/) measured the actual pose shift for each
of the three swaps:

| ligand | pose shift |
|---|---|
| OCNT-2313300 | 0.50 Å |
| OCNT-2395328 | 0.66 Å |
| **OCNT-2395415** | **3.51 Å** |

Two barely moved. One reoriented by 3.5 Å. With ~10 complexes scored, a single structure
that gets ~1.7 Å worse raises mean BiSyRMSD by exactly the +0.17 Å observed, and
−0.0258 × 10 = −0.26 LDDT-PLI is one complex going from decent to bad. LDDT-LP did not
move, confirming the protein is identical and only the ligand placement changed.

**The whole regression is one ligand.** The selection rule is not the variable that moved
the score; luck is.

That cuts both ways — it also means the 05 → 06 result (+0.0012 LDDT-PLI from a
2-structure change, which got 06 logged as "best so far" and prompted "use PL-iPTM going
forward") was well inside the same noise floor and should not have been read as a signal.
mean(PL,HL) and PL-only are not meaningfully different rules.

## Takeaway

Selection-metric A/Bs on a ~10-structure scored subset cannot resolve differences this
small, and are not worth further submission slots. The measurable problem is upstream:
within one ligand's pool the poses differ by ~2 Å median while PL-iPTM varies by 0.01–0.03.
Any improvement has to come from a different *kind* of signal, not a different weighting
of the same confidences.
