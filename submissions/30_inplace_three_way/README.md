# Submission 30 — the three-model hybrid, redone with Glide `inplace`

Submission 24 picked, for each of the 20 ligands, whichever of AF3 (sub 06),
Boltz-2 (sub 18) and ESMFold2 (sub 23) scored best on Glide gscore, and lost:
0.4348 against 0.4527 for shipping ESMFold2 on its own. That looked like
evidence against cross-model selection.

It was not. Every losing Glide selection in this project — 13, 14, 20, 24 — used
`DOCKING_METHOD mininplace`, which minimises the pose before `POSTDOCK` rescores
it. Analysis 28 measured what that does: it moves the pose off the co-folded
geometry and gives the *worse* gscore on 38 of 60 structures (median +0.177).
The one selection run under plain `inplace`, submission 29, is the best result
of the project at 0.4638.

So this is submission 24's rule with the confound removed: same three sets, same
best-gscore rule, scored `inplace`.

## Result of the rebuild

Scores come from analysis 30 (AF3 and Boltz-2, 40 structures, job 48278134) and
analysis 28 (ESMFold2, already scored both ways). All 20 ligands have all three
sets scored.

| | AF3 (06) | Boltz-2 (18) | ESMFold2 (23) |
|---|---|---|---|
| ligands won | 3 | 8 | 9 |

* **5 of 20 picks change** versus submission 24's `mininplace` picks.
* **11 of 20 differ** from shipping ESMFold2 alone, so this is a genuinely
  different set from submission 23, not a near-copy.

Switching mode moved Boltz-2 from 5 wins to 8 and AF3 from 5 to 3. Three of the
five changes (OCNT-2312554, OCNT-2315379, OCNT-2535519) move *to* ESMFold2.

OCNT-2395789 is worth noting: under `mininplace` neither AF3 nor Boltz-2
produced a usable number for it, so it could not be selected at all. Under
`inplace` all three score, and ESMFold2 wins it by a wide margin (−7.903 against
−4.645 and −5.879).

## What this submission tests

Submission 29 showed `inplace` selection works *within* one model's three runs,
where the candidates are close together. This asks the harder question: does it
work *across* models, where the candidates are far apart and the models differ
in overall quality?

The two outcomes are informative either way:

* **Above 0.4527** — cross-model selection works once the scoring mode is right,
  and submission 24's loss was entirely the mode.
* **Below 0.4527** — the mode was not the whole story, and the cross-model
  problem is real: Glide cannot compare poses built on different receptor
  conformations, even scored as given.

The honest prior is caution. AF3 and Boltz-2 ship at roughly 0.42–0.43 on their
own, below ESMFold2's 0.4527, so 11 of the 20 structures here come from weaker
parents. The rule has to be better than good to overcome that.

## Build

Copy-only, by `scripts/01_build.py`. Every structure is a byte-identical copy of
one already validated in its own submission, verified by sha256 against the
source; PoseBusters is not re-run, since it scores each complex in isolation and
an identical file cannot score differently.

* `results/selection.csv` — per-ligand gscores from all three sets, the pick, and
  whether it changed from submission 24
* `visualize/30_inplace_three_way.pse` — PyMOL session, 20 objects

## Result — 0.4414 / 3.7178 / 0.7808

Scored 2026-10-27. The second of the two outcomes above.

`inplace` is worth **+0.0066** over submission 24's `mininplace` version of the
same rule, so analysis 28's finding transfers. But the set still lands 0.0113
below shipping ESMFold2 alone and 0.0224 below submission 29, which applied the
same rule inside ESMFold2's three runs. The mode was a confound; it was not the
only one.

LDDT-LP is the diagnostic: 0.7808 here against 0.7985–0.8063 for every ESM-only
submission. That metric scores the protein pocket, which belongs to the model
that built it, not to the pose. Eleven of these 20 structures come from AF3 and
Boltz-2, and Glide's gscore cannot see that their pockets are worse — it reports
how well a ligand fits the receptor handed to it, so the best fit inside a worse
pocket still wins the comparison and still loses on the metric.

Conclusion for the rest of the challenge: **select within the best model, not
across models.** Within one model the candidates share a pocket quality and
gscore compares like with like, which is why submission 29 worked.
