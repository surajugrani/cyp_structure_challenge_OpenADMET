# 12 — Does a Glide metal constraint fix Fe coordination in redocking?

**Answer: no.** Not a single structure reached 2 Å with or without it, and applied blindly
the constraint destroys cases Glide already got right.

## Motivation

Analysis 04 found Glide docking systematically pushes ligands off the haem iron (median
+0.73 Å vs crystal, 20/22 structures, p = 0.0003), and its metal reward only fires below
~2.6 Å, so just 13 of 214 docked poses earned any. Glide's metal constraint turns that
weak reward into a hard filter.

## Design

10 structures, 2 arms, control built in. The unconstrained arm reuses analysis 05's
redocking (same grids, same LigPrep inputs), so only the constrained arm was run.

| arm | structures | crystal Fe distance |
|---|---|---|
| coordinating | 8SPD, 4I4G, 6MA7, 6DAJ, 6UNI | 1.94–2.40 Å |
| non-coordinating | 6UNJ, 3UA1, 7LXL, 6OO9, 4D6Z | 2.82–6.09 Å |

## Results

| pdb | arm | uncons top-1 | cons top-1 | uncons best | cons best |
|---|---|---|---|---|---|
| 8SPD | coord | 5.33 | 3.95 | 3.34 | 1.53 |
| 4I4G | coord | 3.88 | 5.68 | 2.97 | 3.45 |
| 6MA7 | coord | 4.68 | 5.12 | 3.40 | 2.85 |
| 6DAJ | coord | 6.98 | 5.43 | 4.33 | 2.13 |
| 6UNI | coord | 7.91 | 4.68 | 3.04 | 4.39 |
| 6UNJ | non | 6.60 | 6.10 | 6.44 | 6.10 |
| **3UA1** | non | **0.56** | **no poses** | 0.53 | — |
| 7LXL | non | 6.07 | 4.90 | 2.25 | 4.86 |
| **6OO9** | non | 2.97 | **no poses** | 2.44 | — |
| 4D6Z | non | 5.84 | 3.54 | 5.00 | 3.54 |

Median top-1: coordinating 5.33 → 5.12, non-coordinating 6.07 → 4.90.
**Under 2 Å: 0/10 either way.**

## Why this is a "no"

**The control arm fails.** The constraint helped non-coordinating ligands as much as
coordinating ones — the "it only ever helps" signature indicating it shrinks the search
space rather than capturing Fe chemistry.

**It destroys correct answers.** 3UA1 was the one structure Glide got right unconstrained
(0.56 Å) and the constraint rejected it outright, as it did 6OO9 (2.97 Å). Two of five
non-coordinating structures produced no poses — honest self-selection, but blind
application converts successes into total failures.

**It discards good poses where it does run.** Best-of-N worsened for 7LXL (2.25 → 4.86 Å)
and 6UNI (3.04 → 4.39 Å).

Fe coordination was not the bottleneck — consistent with analysis 09, where Glide's
weakness on CYP3A4 is pose generation and ranking generally.

## The syntax, since it took six attempts

Grid side works as documented:
`generate_glide_grids -rec_file rec.mae -lig_asl "res. LIG" -metal_cons "fe#<atomnum>" -j grid`

Docking side, from the official Glide keyword reference:

> `USE_CONS`: "Each constraint is specified in the format `title[:n]` ... **If n is not
> given, the default features are used.**"
> `[CONSTRAINT_GROUP:m]`: "**All constraints must be defined as part of a group.**"

So: a bare label, inside a group, with **no** FEATURE block.

| attempt | result |
|---|---|
| `USE_CONS fe:1`, no FEATURE block | `feature '1' is not defined` |
| bare `USE_CONS fe` outside a group | silently zero poses |
| group + `fe:1` + hand-written FEATURE | `list index out of range` |
| group + `fe(metal):1` + FEATURE | `list index out of range` |
| xglide `DOCK_USE_CONS fe` | zero poses |
| **group + bare `fe`, no FEATURE** | **works** |

The Maestro GUI always emits explicit FEATURE blocks because it supports non-default
features; reading its serialiser led me to think they were mandatory. Schrodinger notes
the feature file "is not designed to be written from scratch". The keyword reference
answers this in one sentence — it should have been the first thing consulted.

`GLIDECONSUSEMET True` sets the cutoff to r(metal)+r(ligand)+0.4 Å rather than a flat
2.4 Å, relevant since crystal Fe distances here span 1.94–2.40 Å.
