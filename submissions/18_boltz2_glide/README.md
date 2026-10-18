# Submission 18 — Boltz ensemble selected by Glide gscore

**Result: 0.4174 LDDT-PLI / 3.7523 Å BiSyRMSD / 0.8047 LDDT-LP.**

**Best BiSyRMSD of any submission in the project.**

## This contradicted the recommendation

Analysis 23 measured that Glide's picks agree with PL-iPTM at chance (4/20 vs ~4/20
expected), and on that basis — plus submission 13's −0.056 when Glide selected the AF3
ensemble — the advice was not to spend a slot here. **That was wrong.**

| Boltz selection rule | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|
| 15 Boltz `confidence_score` | 0.4063 | 4.0017 | 0.7983 |
| 16 PL-iPTM | **0.4229** | 3.9283 | **0.8060** |
| 17 complex_iPDE | 0.4157 | 3.8021 | 0.7976 |
| **18 Glide gscore** | 0.4174 | **3.7523** | 0.8047 |

Two inferences failed:

1. **Glide's failure on AF3 did not transfer to Boltz.** It cost 0.056 there; here it
   costs 0.0055 against the best rule and wins outright on BiSyRMSD.
2. **Chance-level agreement with PL-iPTM does not imply poor selection.** Two rules can
   disagree on *which* sample to take and both land on good poses — agreement with a
   working rule is not a proxy for performance. The leaderboard was the only way to know.

## The two metrics reward different things

BiSyRMSD improves monotonically across the four rules — 4.0017 → 3.9283 → 3.8021 →
3.7523 — while LDDT-PLI peaks at submission 16. Analysis of the Fe distances shows the
mechanism: although Glide's picks differ from iPDE's on 16 of 20 samples, only **6 of 20**
move the ligand more than 0.05 Å relative to the haem, and all six are **non-coordinating
(type I)** ligands that Glide consistently pulls closer:

| ligand | 17 iPDE | 18 Glide |
|---|---|---|
| OCNT-2395328 | 4.31 Å | 3.01 Å |
| OCNT-2395342 | 4.72 Å | 3.62 Å |
| OCNT-2312966 | 4.24 Å | 3.79 Å |
| OCNT-2315379 | 4.23 Å | 3.81 Å |
| OCNT-2395712 | 4.11 Å | 3.80 Å |
| OCNT-2395645 | 4.15 Å | 3.93 Å |

Coordination is untouched — 11 of 15 capable ligands in both, median 2.29 Å. Glide never
makes or breaks a metal contact; it repositions the type I ligands, and that is where the
BiSyRMSD gain comes from.

This also qualifies submission 11's conclusion. There, moving three ligands closer to the
iron by hand *hurt*. Here, a score that happens to prefer closer type I poses *helps*. The
difference is that Glide is optimising a physical interaction model, not a distance
heuristic.

## Build

No new compute: analysis 23 had already PrepWizard-ed all 100 Boltz samples, so all 20
prepped structures were reused. Selection is the most negative `gscore` per ligand.

Differs from submission 15 in 14/20, from 16 in 16/20, from 17 in 16/20.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: chosen model, gscore, number scored, gscore spread |
| `raw_sel/` | the selected Boltz model per ligand, haem renamed to HEM |
