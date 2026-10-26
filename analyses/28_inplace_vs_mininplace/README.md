# Analysis 28 — `inplace` vs `mininplace`, and a third chance-level selector

**Question.** Every Glide-based selection shipped here (submissions 13, 18, 19, 20, 21,
24) used `DOCKING_METHOD mininplace`, on one argument: plain `inplace` returns emodel and
einternal as 0.00, and emodel is Glide's own pose-ranking term. That had never been
tested against a leaderboard metric, and testing it by shipping two submissions would
cost two attempts.

Instead: score the three ESMFold2 sets both ways and ask whether the modes would even
pick different structures. 60 structures × 2 modes, one grid per structure so the mode is
the only thing that varies.

## The modes are not interchangeable

| | |
|---|---|
| `inplace` emodel | 0.0 for every structure — the original argument holds |
| Pearson r between the two modes' gscore | **0.791** |
| median `mininplace` − `inplace` gscore | +0.177 (range −0.93 to +1.92) |
| **selections that differ** | **7 of 20** |

So the mode materially changes which structure gets picked. A comparison of the two as
submissions would be measuring something real.

## But neither mode can select, and neither can RTCNN

The three sets have known scores, so "did the selector pick well?" is answerable without
submitting anything:

| set | LDDT-PLI |
|---|---|
| 23 default endpoint + ColabFold MSA | **0.4527** |
| 26 max_accuracy, own MSA | 0.4365 |
| 27 default endpoint, no MSA | 0.4384 |

| selector | picks 23 | picks 26 | picks 27 |
|---|---|---|---|
| Glide `mininplace` | 7 | 7 | 6 |
| Glide `inplace` | 6 | 9 | 5 |
| ICM RTCNN (analysis 29) | 6 | 6 | 8 |

Chance is 6.67 each. **All three land on it.** Submission 23 is the better set by about
0.015 LDDT-PLI and no selector prefers it above chance.

## Conclusion

The `inplace` vs `mininplace` question is real but not worth two attempts: it changes
7 of 20 picks, and both modes select at chance, so a head-to-head would compare two
coins. Keep `mininplace` — the emodel argument is the only discriminating reason either
way.

More importantly this is the third independent selector to fail on these three sets,
after analysis 23 (Glide on the Boltz ensemble, agreement with PL-iPTM at chance) and
submission 24 (Glide across three models, 0.0179 below its best parent). The sets sit a
median 1.53 Å apart by pocket-superposed ligand RMSD (analysis 27), the within-ensemble
regime where every selector tried here has failed.

## Reproduce

```bash
sbatch analyses/28_inplace_vs_mininplace/scripts/run.sh   # resumable
```

Ran in 4 h + 2 h 47 m serially on 1 CPU, peak 2.1 GB. Results in
`results/scores_both_modes.csv`.
