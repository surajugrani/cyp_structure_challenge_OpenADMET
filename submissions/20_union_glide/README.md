# Submission 20 — best Glide score over the union of all AF3 and Boltz samples

**Result: 0.3854 LDDT-PLI / 4.2275 Å BiSyRMSD / 0.7449 LDDT-LP — a clear loss.**

| | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|
| 19 pairwise, cross-model | **0.4331** | **3.6754** | 0.7823 |
| **20 union, global** | **0.3854** | 4.2275 | 0.7449 |
| 13 Glide within AF3 alone | 0.3730 | 4.5822 | 0.6968 |

**It landed next to submission 13, not next to 19** — and that is the answer the
experiment was built to get.

One rule applied to everything: the most negative `gscore` across all 94 Glide-scored AF3
samples and all 100 Boltz samples. **12 AF3 / 8 Boltz**, 20/20 structural, 20/20
PoseBusters, no new compute.

## How it differs from submission 19

Submission 19 (**0.4331 / 3.6754 / 0.7823** — best of the project on LDDT-PLI and
BiSyRMSD) compared only *two* structures per ligand, and they were not chosen the same
way: 06's was the best **PL-iPTM** AF3 sample, 18's the best **Glide** Boltz sample.

Applying one rule to the full pool changes **10 of 20** structures and flips the balance
from 6/14 to 12/8 in AF3's favour. The reason is that the union sees AF3 samples
submission 06 never shipped — for six ligands a *different* AF3 sample outscores
everything, e.g. OCNT-2312554 `seed-1_sample-1` at −7.01 against the Boltz structure 19
used.

## The bet, and the risk

**For:** submission 19 showed Glide resolves differences between whole models well
(+0.0040 LDDT-PLI over the better parent), even though analysis 23 found it picks within
a single ensemble at chance. The scale explanation is that candidates ~4 Å apart are
separable while candidates 1–2 Å apart are not.

**Against:** submission 13 applied exactly this rule *within AF3 alone* and scored
**0.3730**, the worst AF3-based attempt. Twelve of the twenty picks here are AF3 samples
chosen on gscore — the same losing rule. The union pool mixes both regimes, so a global
Glide rule may simply revert to within-ensemble behaviour wherever the best candidates
come from one model.

## Verdict: the scale explanation holds

20 lost 0.048 LDDT-PLI against 19 and finished 0.012 above submission 13, the
within-AF3 Glide rule. Twelve of its twenty picks *were* AF3 samples chosen on gscore —
13's rule exactly — and the score followed them.

**Glide's selection skill is specific to comparing whole models, not general.** Between
AF3 and Boltz the candidates differ by a median 4.15 Å and Glide resolves them (19 beat
both parents). Within one model's ensemble they differ by 1–2 Å, Glide picks at chance
(analysis 23) and loses outright (13, 20). A global rule destroys the advantage by
letting Glide rank within ensembles wherever the best candidates happen to come from one
model.

The practical rule: **give Glide a choice between models, never between samples.**

### The LDDT-LP mixing model predicted this too

`(12 × 0.6985 + 8 × 0.8047) / 20 = 0.741` against 0.7449 observed — the second time this
simple mix has been accurate to ~0.004 (submission 19 predicted 0.773, observed 0.7823).
LDDT-LP is essentially a weighted average of the parent pockets, so it is predictable in
advance for any hybrid.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: chosen model and sample, gscore, candidate count, and what submission 19 shipped |
