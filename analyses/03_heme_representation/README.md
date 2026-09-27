# 03 — Heme representation, Glide's metal term, and score sensitivity

**Question:** Analyses 01 and 02 found Glide's metal term is identically 0.00 under the
chemically correct heme, while 16/20 ligands coordinate the iron. Can the term be
reactivated *without* misrepresenting the chemistry — and how much does heme
representation affect docking scores generally?

**Status:** Gate test complete on 3 structures (8SPD, 4K9X, 2V0M). `metal_cons` route pending.

## Why this matters

GlideScore term weights, fitted by least squares over 212 redocked poses from analysis 02:

| Term | Weight |
|---|---|
| esite | +2.51 |
| lipo | +1.19 |
| rewards | +1.05 |
| hbond | +0.41 |
| ecoul | +0.19 |
| evdw | +0.046 |

The metal term enters at ~1.0 (unestimable here — zero variance in the baseline data).
It is therefore **decisive at the relevant scale**: for 8SPD, drifting off the iron gains
~21 kcal of raw vdW = ~1.0 GlideScore unit, while a −2.30 metal term would outweigh it.

## Method

Each variant is applied to the already-prepped complex from analysis 01 (no PrepWizard
re-run), the **grid is regenerated** (metal perception is fixed at grid generation), and
the crystal pose is rescored in place. **All 13 Glide terms are recorded**, so a change in
the total score can be attributed rather than guessed at.

```bash
sbatch scripts/run_variants.sh --pdbs 8SPD,4K9X,2V0M
```

## Results

Identical across all three structures:

| Variant | Fe | Porphyrin | Heme charge | **Metal** | 8SPD score | 8SPD ecoul |
|---|---|---|---|---|---|---|
| `fe2_dianion` — analyses 01/02 baseline | +2 | dianion | −2 | 0.00 | −4.78 | −1.00 |
| `fe3_dianion` — **chemically correct** | +3 | dianion | −1 | **0.00** | −3.37 | **+7.02** |
| `fe2_monoanion` | +2 | monoanion | −1 | 0.00 | −4.52 | −1.83 |
| `fe3_monoanion` | +3 | monoanion | 0 | **−2.30** | −7.62 | **−7.12** |

### 1. The chemically correct heme does not activate the metal term

CYP3A4's resting state is ferric: Fe(III) + porphyrin dianion + 2 deprotonated propionates
= net −1. That is `fe3_dianion`, and it gives **metal = 0.00**. Setting the physiological
oxidation state does not help.

### 2. Activation is not explained by net formal charge

`fe3_dianion` (−1) and `fe2_monoanion` (−1) both give 0.00, yet 8SPD "state B" from
analysis 01 — also net −1, but with a pyrrole **N–H** present — gave −2.30. Same net
charge, different outcome. Glide's metal-site perception depends on something more
specific, which we could not determine from outside the program. **Left unresolved rather
than guessed at.**

### 3. The metal term looks like a fixed reward, not a distance-graded one

Every activated case gives **exactly −2.30**: 8SPD states A and B from analysis 01, and
all three `fe3_monoanion` structures here — across Fe distances of 1.94, 1.97 and 2.10 Å
and three different ligands.

**Consequence:** even a working metal term would *not* rescue analysis 01. A constant
−2.30 applied to every coordinating ligand shifts all scores equally and adds no
discrimination between them. It would likely help **pose recovery** (analysis 02) by
anchoring the coordinating atom, but not **affinity ranking**.

### 4. Heme representation shifts scores by more than the metal term itself

`ecoul` for 8SPD swings from −1.00 → +7.02 → −7.12 across variants: a 14-unit raw range,
≈2.7 GlideScore units at weight 0.19. So `fe3_monoanion`'s better total (−7.62 vs −4.78)
comes from **both** the metal term (−2.30) and electrostatics (~−1.2), not the metal term
alone.

This is arguably the most transferable result here: **docking scores for heme proteins are
strongly sensitive to cofactor representation**, on a scale comparable to the entire
dynamic range being correlated against affinity. Pipelines rarely report which
representation they used.

## Decision

Keep the **chemically correct** heme. The only activation route found so far
(`fe3_monoanion`) requires a monoanionic porphyrin, which is wrong for a metalated
porphyrin — and per point 3 it would not improve affinity ranking anyway.

Analyses 01 and 02 are deliberately **left intact** as the Fe(II) dianion baseline rather
than being re-run under a different representation; the sensitivity shown here is a more
useful result than any single set of numbers.

## Next

Test Glide's `metal_cons` constraint at grid generation. That imposes coordination
geometrically **without** misrepresenting the chemistry, and is the remaining legitimate
route to better pose recovery. Prediction recorded in advance: it should improve pose
RMSD while leaving the metal term at 0.00 — i.e. a geometric fix, not a scoring fix.
