# Submission 14 — AF3 samples selected by ICM RTCNN

**Result: 0.3949 LDDT-PLI / 4.4057 Å BiSyRMSD / 0.7004 LDDT-LP** against submission 06's
0.4291 / 4.0647 / 0.6985.

Better than Glide selection (13, 0.3730) by 0.022, still 0.034 below AF3's own PL-iPTM.

## Why RTCNN was the one to try

Analysis 09 ranked pose-scoring terms by within-ligand Spearman against RMSD — the
quantity a selector needs:

| term | ρ |
|---|---|
| **ICM RTCNN** | **+0.334** |
| Glide lipo | +0.115 |
| Glide gscore | +0.086 |
| Glide emodel | +0.057 |

Submission 13 tested the Glide end and lost 0.056. RTCNN was the one selector with real
measured skill behind it.

## Outcome

**The ranking transferred; the level did not.** RTCNN beat Glide by exactly the direction
analysis 09 predicted, but both lost to AF3's own confidence:

| selector | LDDT-PLI |
|---|---|
| PL-iPTM (06) | **0.4291** |
| mean(PL,HL) (05) | 0.4279 |
| HL-iPTM (10) | 0.4221 |
| ICM RTCNN (14) | 0.3949 |
| Glide gscore (13) | 0.3730 |

16 of 20 structures differ from 06, yet the Fe chemistry barely moves — 11 of 15 capable
ligands coordinate against 06's 12, median 2.34 Å against 2.25 Å. The swaps are lateral:
a different sample in the same binding mode, losing 0.034 LDDT-PLI with nothing visible
changing.

`RTCNN` and analysis 09's preferred `RTCNN+Hbond` picked identically on all 20 ligands, so
the Hbond term added nothing and the simpler rule was used.

## Build notes

Shipped coordinates are the **original PrepWizard poses** — ICM is used only to rank.
`_dockScan -n` minimises before scoring (~0.3–1.4 Å), so its own output coordinates were
deliberately not used.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: chosen sample, RTCNN, PL/HL iPTM, whether it differs from 05 and 06 |
