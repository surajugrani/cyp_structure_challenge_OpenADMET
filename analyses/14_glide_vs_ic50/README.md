# Analysis 14 — do Glide scores of the AF3 poses predict potency?

**No.** Not the total score, not the vdW term, not a fitted combination of vdW and
electrostatics. Every cross-validated model predicts a held-out ligand worse than guessing
the mean.

This closes a question the user asked to be reminded about since 2026-09-30.

## Setup

15 of the 20 blinded structure ligands also carry a real `CYP3A4_pIC50_direct_inhibition`
value in `input/cyp-challenge_cyp-inhibition-data.csv`, so the test runs directly on the
blinded set. Glide scores come from the rescoring of the 100 PrepWizard-minimised AF3
samples (submission 04).

Two pose-selection rules, because which pose you score is itself a choice:

| rule | which sample |
|---|---|
| 1 | best (most negative) Glide value for that ligand |
| 2 | the sample submission 06 ships — best PL-iPTM among PoseBusters-passing |

## Result 1 — correlations

| x | rule | Spearman | p | Pearson | p |
|---|---|---|---|---|---|
| gscore | 1 best score | +0.268 | 0.33 | +0.118 | 0.67 |
| gscore | 2 best PL-iPTM | +0.032 | 0.91 | −0.108 | 0.70 |
| evdw | 1 best evdw | −0.029 | 0.92 | +0.054 | 0.85 |
| evdw | 2 best PL-iPTM | −0.068 | 0.81 | −0.035 | 0.90 |

**The gscore sign is backwards.** Glide is negative-is-better, so a real affinity signal
would give a negative Spearman. Both gscore rows are positive — better-scoring poses
belong to the slightly *less* potent compounds. Same direction analysis 01 found on
crystal poses and analysis 08 found for ICM.

**Rule 1 is not a fair comparison.** The two rules pick the same sample for only 1 of 15
ligands on gscore, and rule 1 takes a *minimum over ~5 samples* — an order statistic whose
downward bias depends on how spread that ligand's pool is, and analysis 13 showed those
spreads vary from 0.21 to 5.48 Å. Rule 1's marginally larger rho is more likely that
selection artefact than signal. Rule 2 is the honest number: +0.03.

**evdw is mostly molecular weight.** It spans −31 to −47 kcal/mol against gscore's −5.3 to
−7.9; the big-|evdw| compounds simply bury more atoms. Over 240–412 Da and ~2 log units of
potency, size is not tracking potency — which also explains gscore's spurious positive
rho, since gscore contains evdw. Consistent with analysis 09, evdw barely varies within a
ligand's pool, so the two rules agree far more often here (5/15 vs 1/15).

## Result 2 — fitted models, and why the fit is a mirage

pIC50 ~ f(evdw, ecoul), fitted by least squares, scored by leave-one-out Q².

| rule | model | params | R² in-sample | **Q² (LOO)** | RMSE (LOO) |
|---|---|---|---|---|---|
| best PL-iPTM | mean only | 1 | — | −0.148 | 0.539 |
| best PL-iPTM | linear | 3 | 0.004 | **−0.484** | 0.613 |
| best PL-iPTM | quadratic | 6 | 0.551 | **−0.100** | 0.527 |
| best gscore | linear | 3 | 0.010 | −0.377 | 0.590 |
| best gscore | quadratic | 6 | 0.305 | −1.249 | 0.754 |

Q² ≤ 0 means the model predicts a held-out ligand no better than the mean of the other 14.
**Every Q² is negative.**

The quadratic is the instructive case: R² = 0.551 scored on its own training points, and
Q² = −0.100 the moment each ligand is predicted from a fit that excluded it. That gap is
6 parameters memorising 15 points. Its in-sample R² also swings 0.551 → 0.305 depending
only on which pose is fed in, while both Q² stay negative — a real relationship would not
be that sensitive to pose selection.

The linear model is at least honest about itself (R² = 0.004). Its coefficients are tiny
against the data: evdw −0.0023 and ecoul −0.0161 per kcal/mol, so the full ~16 kcal/mol
evdw range moves predicted pIC50 by ~0.04 log units against a pIC50 sd of 0.52. It
outputs 5.52–5.86 for ligands whose true values run 4.96–6.94 — a flat line at the mean.

## Verdict

Treat Glide-terms-vs-potency as closed. Across analyses 01, 08 and 14, ~30 correlations
plus these five models have been tested against the same 15 pIC50 values over ~2 log
units, with nothing surviving cross-validation. The binding constraint is n and dynamic
range, and no change of term, model class or pose-selection rule fixes that — only more
compounds or a wider potency spread would.

## Notes for reuse

* `scripts/02_corr_term.py` takes any GlideScore term as `argv[1]` (lipo, hbond, …).
* `scripts/03_fit_model.py` fits mean / linear / quadratic with LOO Q².
* scipy and sklearn are **not** in `cyp_posebusters_env`; Spearman, Pearson, the t-test
  tail and all model fitting are implemented on numpy in `01_corr.py`.
* Glide writes `10000.0` as a "no score" sentinel and can pair it with `status ok` — one
  OCNT-2395789 sample does. Filter on magnitude, not just status.

| file | contents |
|---|---|
| `results/glide_vs_ic50.csv` | per ligand: pIC50, both rules' sample and gscore |
| `results/term_vs_ic50_evdw.csv` | same for evdw |
| `results/model_fit.csv` | per ligand in-sample and LOO predictions for each model |
