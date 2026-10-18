# Submission 19 — AF3/Boltz hybrid, picked by Glide score

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
