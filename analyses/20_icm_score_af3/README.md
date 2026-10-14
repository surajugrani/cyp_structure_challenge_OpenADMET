# Analysis 20 — ICM scoring of the AF3 samples, and the end of the selector question

Scored all 100 PrepWizard-minimised AF3 samples in place with ICM (`_dockScan -n` on each
sample's own receptor conformation), then compared what each ICM term would select
against submission 06.

Together with submission 13 (Glide) this closes the question: **can any external scoring
function pick AF3's best sample better than AF3's own confidence?** No.

## The answer

| selector | LDDT-PLI | submission |
|---|---|---|
| PL-iPTM | **0.4291** | 06 |
| mean(PL, HL) | 0.4279 | 05 |
| HL-iPTM | 0.4221 | 10 |
| **ICM RTCNN** | **0.3949** | 14 |
| **Glide gscore** | **0.3730** | 13 |

Two findings.

**AF3's own confidence wins.** The three iPTM rules sit within 0.0070 of each other; both
external scores sit 0.03-0.06 *below* all of them. A physics or ML score computed on the
final pose carries less information about which sample is right than the model's own
confidence did at generation time.

**The analysis 09 ranking transferred.** That analysis measured within-ligand Spearman
against RMSD on crystal receptors — RTCNN +0.334, Glide gscore +0.086 — and predicted
RTCNN should be the better selector. It is, by 0.022 LDDT-PLI. The relative ordering held
when moving from crystal to co-folded receptors; the absolute level was simply too low to
beat AF3.

## What the rules actually change

94 of 100 samples were both scored and PoseBusters-passing; 1 failed ICM scoring.

| rule | differs from 06 | coordinating (of 15 capable) | median Fe |
|---|---|---|---|
| 06 PL-iPTM | — | **12/15** | 2.25 Å |
| ICM `Score` | 13/20 | 11/15 | 2.30 Å |
| ICM `RTCNN` | 16/20 | 11/15 | 2.34 Å |
| `RTCNN+Hbond` | 16/20 | 11/15 | 2.34 Å |
| `mfScore` | 17/20 | 11/15 | 2.28 Å |

Every rule reshuffles most of the set while leaving the Fe chemistry essentially
untouched. The swaps are **lateral** — a different sample with the same binding mode —
which is why they can lose 0.03-0.06 LDDT-PLI without any visible change in coordination.

`RTCNN` and analysis 09's preferred `RTCNN+Hbond` pick **identically on all 20 ligands**,
so the Hbond term contributes nothing here and the simpler rule was submitted.

## Method notes

* `_dockScan -n` scores the supplied pose through the normal VLS pipeline, but **minimises
  first** — the scored pose sits ~0.3-1.4 Å from the input (median 0.86 Å, analysis 08).
  Submission 14 ships the ORIGINAL PrepWizard coordinates; ICM only ranks.
* Each of the 100 samples is a different receptor conformation, so each needs its own
  `_dockBatch` project. One at a time, 60 minutes total.
* ICM writes the SD tag `RTCNN`, not `RTCNNscore` — analysis 06 had this wrong and
  silently produced an empty column.

| file | contents |
|---|---|
| `results/icm_sample_scores.csv` | all 100 samples: Score, RTCNN, mfScore, Hbond, Hphob, Dsolv, Eintl, Fe distance |
| `results/rule_comparison.csv` | per rule and ligand: chosen sample, 06's choice, whether it differs |
