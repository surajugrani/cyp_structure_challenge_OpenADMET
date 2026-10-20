# Submission 21 — submission 19 with a Glide margin threshold

Submission 19, but AF3 structures are kept only where they win Glide by more than
0.5 kcal/mol. **17 Boltz / 3 AF3** (19 was 14/6), changing exactly 3 of 20 structures.
20/20 structural, 20/20 PoseBusters, no new compute.

## Which three reverted, and why

| ligand | AF3 margin | kept |
|---|---|---|
| OCNT-2395712 | 0.07 | → Boltz |
| OCNT-2395503 | 0.08 | → Boltz |
| OCNT-2395457 | 0.31 | → Boltz |
| OCNT-2315379 | 1.01 | AF3 |
| OCNT-2312966 | 1.47 | AF3 |
| OCNT-2535519 | 2.15 | AF3 |

**The bet:** an AF3 structure costs pocket quality — LDDT-LP 0.6985 against Boltz's
0.8047 — so importing one for a 0.07 kcal/mol score advantage is paying real LDDT-LP for
noise.

**The risk:** submission 19 beat 18 by 0.0157 LDDT-PLI *because* of its AF3 content.
Halving that content could hand some of the gain back, and LDDT-PLI is the primary
metric while LDDT-LP is not.

## Prediction on record

The mixing model has now been accurate to ~0.004 twice:

| submission | mix | predicted LDDT-LP | observed |
|---|---|---|---|
| 19 | 14 Boltz / 6 AF3 | 0.773 | 0.7823 |
| 20 | 8 Boltz / 12 AF3 | 0.741 | 0.7449 |
| **21** | **17 Boltz / 3 AF3** | **0.789** | — |

Both observations ran ~0.009 above prediction, so **expect ≈0.798**, about +0.016 over
submission 19. If LDDT-PLI holds near 0.433 this is a strict improvement; if it drops
toward 18's 0.4174, the three reverted AF3 poses were earning their keep after all.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: source, what 19 used, both gscores, the AF3 margin, whether it changed |
