# Submission 20 — best Glide score over the union of all AF3 and Boltz samples

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

If 20 beats 19, the scale story is wrong and Glide is just a good selector given enough
candidates. If it lands below, the scale story holds and cross-model comparison is the
specific thing Glide is good at.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: chosen model and sample, gscore, candidate count, and what submission 19 shipped |
