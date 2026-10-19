# Submission 19 — AF3/Boltz hybrid, picked by Glide score

**Result: 0.4331 LDDT-PLI / 3.6754 Å BiSyRMSD / 0.7823 LDDT-LP.**

**Best of the project on both LDDT-PLI and BiSyRMSD**, beating both parents on the
primary metric (06: 0.4291, 18: 0.4174).

| submission | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|
| 06 AF3 + PL-iPTM | 0.4291 | 4.0647 | 0.6985 |
| 16 Boltz + PL-iPTM | 0.4229 | 3.9283 | **0.8060** |
| 17 Boltz + iPDE | 0.4157 | 3.8021 | 0.7976 |
| 18 Boltz + Glide | 0.4174 | 3.7523 | 0.8047 |
| **19 hybrid** | **0.4331** | **3.6754** | 0.7823 |

## Glide can choose between models even though it cannot rank poses within one

This is the result the submission was built to test. Analysis 23 showed Glide picks
samples from a single ensemble at chance, and submission 13 showed it *loses* 0.056 doing
so on AF3's. Yet choosing between two whole co-folding models it beats both — +0.0040
over the better parent on LDDT-PLI and −0.077 Å on BiSyRMSD.

The likely reason is scale. Within one model's five samples the candidates differ by
~1–2 Å and Glide's score differences are noise; between AF3 and Boltz they differ by a
median 4.15 Å, which is a large enough gap for a physics score to resolve.

## The LDDT-LP drop is the pocket, not the ligand

LDDT-LP scores the **protein binding site**, not ligand placement. The value here is
fully explained by mixing in AF3's weaker pocket:

```
(14 × 0.8047  +  6 × 0.6985) / 20  =  0.773      observed 0.7823
```

So nothing about submission 19's ligands is implicated — the cost is simply the price of
importing 6 AF3 structures for their better ligand poses. An all-Boltz variant would
recover ~0.02 LDDT-LP but give up whatever those 6 AF3 poses contributed to LDDT-PLI.

For each ligand, whichever of submission 06 or 18 Glide scores better.

| source | | count |
|---|---|---|
| 18 | tuned Boltz + Glide selection — 0.4174 / **3.7523** / 0.8047 | **14/20** |
| 06 | AF3 + PL-iPTM selection — **0.4291** / 4.0647 / 0.6985 | 6/20 |

AF3 wins on OCNT-2312966, -2315379, -2395457, -2395503, -2395712, -2535519.
20/20 structural, 20/20 PoseBusters. No new compute — both parents were already
PrepWizard-minimised and validated, so the structures are copied verbatim.

**OCNT-2395789 takes the Boltz structure by instruction**: its chosen AF3 sample is one
of three that carried Glide's `10000.0` "no score" sentinel, so no comparable score
exists.

## Two caveats

**The two gscores come from different receptors.** 06's is computed on AF3's protein,
18's on Boltz's, each with its own grid. Glide scores are not calibrated for
cross-receptor comparison, so "better gscore" is a weaker claim than it appears.

**The parents disagree enormously and Glide does not track the disagreement.** Median
ligand RMSD between them is 4.15 Å (max 7.71 Å, 17 of 20 beyond 2 Å), yet the size of
that disagreement is unrelated to which pose Glide prefers — Pearson **+0.110**, Spearman
**+0.018** over n=19. The rule is choosing between genuinely different answers with a
score that is blind to the difference.

Two coordinating ligands make the point: OCNT-2313300 is 6.21 Å apart between the two
submissions while both sit at 2.13 Å from the iron, and OCNT-2535519 is 4.90 Å apart with
both at ~2.10 Å. Same metal contact, the molecule pivoted around its coordinating
nitrogen.

## Why run it anyway

Glide selection *within* the Boltz ensemble produced the best BiSyRMSD of the project
(submission 18), overturning analysis 23's recommendation against trying it. Predictions
about what Glide can and cannot do here have been wrong in both directions, so the
leaderboard is the only arbiter.

What makes this a genuinely different bet from 18: it mixes two **co-folding models**
rather than re-ranking one model's samples. Since the parents agree on almost nothing yet
both score ~0.42 LDDT-PLI, roughly half of each one's apparent accuracy must come from
different structures — and if Glide can tell which half, the hybrid should beat both.
If it lands between them, Glide cannot.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: source submission, both gscores, their difference |
| `../../figures/rmsd_vs_delta_glide.svg` | the pose disagreement vs score difference scatter |
