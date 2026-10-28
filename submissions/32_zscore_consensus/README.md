# Submission 32 — z-score consensus of Glide `inplace` and ICM RTCNN

A third selection rule over the same three-ESMFold2-run pool used by submissions
28 and 29.

| | rule | LDDT-PLI |
|---|---|---|
| 23 | best single parent, no selection | 0.4527 |
| 28 | best ICM RTCNN | 0.4411 |
| 29 | best Glide `inplace` gscore | **0.4638** |
| 32 | consensus of the two | this one |

## Why

The two scorers are nearly decorrelated — analysis 29 found they pick the same
winning structure for only **10 of 20** ligands, Pearson +0.373. Decorrelated
errors are the textbook precondition for consensus beating either input. And
RTCNN is not a junk input: analysis 09 measured its within-ligand Spearman
against true RMSD on crystals at **+0.334**, about four times Glide gscore's
+0.086.

## Why it is probably wrong anyway

The two have already met prospectively on this exact pool and Glide won by
0.0227. Every structure consensus changes relative to submission 29 is a case
of RTCNN overruling Glide — so this imports precisely the disagreements that
submission 28 lost on. The evidence leans against it; it is being run because
an unused submission is worth nothing and this is the version that returns a
readable answer.

## Method

gscore and RTCNN are different units on different scales, so a raw sum means
nothing. Each is z-scored *within* a ligand's three candidates — subtract that
ligand's mean, divide by its population standard deviation — then added. Both
are "more negative is better", so the minimum summed z wins.

With three candidates the standard deviation comes from three points and is
noisy. That is a genuine weakness. Rank-averaging was tried first and has no
resolution at n=3 (ranks 0/1/2 tie constantly), collapsing onto Glide for 18 of
20 — it would have tested nothing.

A safer variant was also computed: Glide primary, RTCNN breaking ties only
below a 0.25 gscore margin. It changes 3 of 20, all with Glide margins of
0.00–0.12 where Glide is genuinely indifferent. It was rejected because ~1–2 of
those fall inside the ~10 scored structures, which is unreadable. Between a
safe unreadable test and a risky readable one, this is the readable one.

## What it selected

| set | ligands won |
|---|---|
| 23 default + ColabFold MSA | 7 |
| 26 max-accuracy endpoint | 5 |
| 27 default, no MSA | 8 |

* scorers agreed outright on **10/20**
* differs from submission 29 (Glide) on **7/20**
* differs from submission 23 alone on **13/20**
* chose a set **neither scorer picked alone** on **3/20** — OCNT-2315379,
  OCNT-2395457, OCNT-2535422. These are the only places in the project where a
  combined rule produced a new choice rather than siding with one input.

## Build

Copy-only. Every structure is copied unchanged from a submission that already
passed the full gate 20/20, verified by sha256. PoseBusters is not re-run: a
byte-identical copy cannot return a different verdict.

* `results/selection.csv` — per-ligand z-scores, both scorers' picks, and what
  changed against submission 29
* `visualize/32_zscore_consensus.pse` — PyMOL session, 20 objects
