# Submission 16 — Boltz-2 ensemble selected on protein–ligand iPTM

**Result: 0.4229 LDDT-PLI / 3.9283 Å BiSyRMSD / 0.8060 LDDT-LP.**

**First submission under 4 Å BiSyRMSD, and the best LDDT-LP of any attempt.** Submission
06 still leads the primary metric by 0.0062.

Same 5-sample Boltz ensemble as submission 15; only the choice of sample changes.

| | 15 Boltz's own pick | **16 PL-iPTM** | change |
|---|---|---|---|
| LDDT-PLI | 0.4063 | **0.4229** | **+0.0166** |
| BiSyRMSD | 4.0017 | **3.9283** | **−0.073 Å** |
| LDDT-LP | 0.7983 | **0.8060** | +0.0077 |

## The first selection rule that actually worked

Every previous attempt to improve selection failed or backfired:

| rule | LDDT-PLI | vs its baseline |
|---|---|---|
| AF3 PL-iPTM / mean / HL-iPTM (06/05/10) | 0.4291 / 0.4279 / 0.4221 | spread 0.0070 — noise |
| AF3 + Glide gscore (13) | 0.3730 | −0.056 |
| AF3 + ICM RTCNN (14) | 0.3949 | −0.034 |
| **Boltz + PL-iPTM (16)** | **0.4229** | **+0.0166** |

**Why this one worked: dynamic range.** A selector can only help if it varies across the
candidates. Median spread over the 5 samples of a ligand:

| metric | spread |
|---|---|
| Boltz `confidence_score` (its own ranking) | 0.0038 |
| **Boltz `pair_chains_iptm[0][2]`** | **0.079** |
| AF3 PL-iPTM across its 5 samples | 0.01–0.03 |

20× more range than the metric Boltz ranks by. The rule changes **14 of 20** structures
away from Boltz's own `model_0`.

## Which number, and why the matrix is asymmetric

Chains are indexed `0` = protein (chain A), `1` = haem (B), `2` = ligand (C). The same
protein–ligand interface therefore appears twice, and the two entries disagree:
`[0][2]` = 0.698 against `[2][0]` = 0.956 for the first ligand.

ipTM is a TM-score-style sum whose distance cutoff `d0` is set by the size of the set
being scored. `[0][2]` normalises over the ~27 ligand atoms — a tight tolerance, so errors
count and the score spreads out. `[2][0]` normalises over ~500 protein residues — a loose
tolerance, pinning every sample near 0.95 and discriminating nothing.

**`[0][2]` is the usable direction**, and it is also the direct analogue of AF3's
`chain_pair_iptm[0][2]`, the metric that produced submission 06. Averaging both directions
would change the pick on 4 of 20; the AF3-consistent choice was kept.

## Cumulative effect

From default Boltz (submission 03, 0.3921):

| step | LDDT-PLI | gain |
|---|---|---|
| 03 defaults | 0.3921 | — |
| 15 tuned sampling (recycling 5, sampling 300, 5 diffusion samples) | 0.4063 | +0.0142 |
| **16 + PL-iPTM selection** | **0.4229** | **+0.0166** |
| | | **+0.0308 total** |

Both halves matter, and neither alone would have been enough. Generating more diverse
candidates is only useful with a selector that can tell them apart, and a selector is only
useful when the candidates differ.

## Build notes

* 6 of 20 PrepWizard results were reused from submission 15 (those ligands where PL-iPTM
  also picked `model_0`); the other 14 were prepped fresh with identical options.
* Boltz labels **every** het chain `LIG`, including the haem — chain B is renamed to `HEM`
  before PrepWizard. Skipping this failed validation 0/20 in submission 15's first attempt.

| file | contents |
|---|---|
| `raw_sel/` | the PL-iPTM-selected Boltz model per ligand, haem renamed |
| `../15_boltz2_tuned/results/pl_iptm_picks.csv` | per ligand: chosen model, PL-iPTM, spread vs confidence_score |
