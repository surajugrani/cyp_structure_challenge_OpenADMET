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

| Variant | Fe | Pyrrole N | N–H | Heme charge | **Metal** | 8SPD score | 8SPD ecoul |
|---|---|---|---|---|---|---|---|
| `fe2_dianion` — analyses 01/02 baseline | +2 | 2 × N⁻ | 0 | −2 | 0.00 | −4.78 | −1.00 |
| `fe2_monoanion` | +2 | 1 × N⁻ | 0 | −1 | 0.00 | −4.52 | −1.83 |
| **`fe2_nh`** — reproduces "state B" | +2 | 1 × N⁻ | **1** | −1 | **−2.30** | −7.07 | −4.22 |
| `fe3_dianion` — **chemically correct** | +3 | 2 × N⁻ | 0 | −1 | 0.00 | −3.37 | +7.02 |
| `fe3_monoanion` | +3 | 1 × N⁻ | 0 | 0 | **−2.30** | −7.62 | −7.12 |
| `fe3_nh` | +3 | 1 × N⁻ | **1** | 0 | 0.00 | −5.00 | −9.51 |

Identical metal values across all three structures tested.

### 1. The chemically correct heme does not activate the metal term

CYP3A4's resting state is ferric: Fe(III) + porphyrin dianion + 2 deprotonated propionates
= net −1. That is `fe3_dianion`, and it gives **metal = 0.00**. Setting the physiological
oxidation state does not help.

### 2. Activation follows no chemically sensible rule

Two of six representations activate the term, and they share no identifiable feature:

- **Net heme charge does not predict it.** `fe2_monoanion` and `fe2_nh` are both −1, giving
  0.00 and −2.30. `fe3_monoanion` and `fe3_nh` are both 0, giving −2.30 and 0.00.
- **Fe–porphyrin core charge does not predict it.** `fe2_monoanion` and `fe2_nh` are both
  +1; `fe3_monoanion` and `fe3_nh` are both +2. Same split.
- **The N–H does not predict it either** — it *activates* the term at Fe²⁺ and
  *deactivates* it at Fe³⁺.

> An earlier version of this README claimed only `fe3_monoanion` activates the term. That
> was **wrong**: `fe2_nh` also does, at the same −2.30, and it reproduces the "state B"
> already observed in analysis 01. The error came from testing `fe2_monoanion` (charge
> neutralised, no H added) as a stand-in for state B and not noticing the proxy differed in
> exactly the respect that mattered.

**Practical caution:** Glide's metal-site perception flips on single-atom changes in a way
that is not derivable from formal charge or protonation state. For metalloprotein docking,
the metal term should be checked empirically per receptor preparation, never assumed.

### 3. The metal term is saturated at crystal geometries — not fixed

Every activated case here gives **exactly −2.30**: 8SPD states A and B from analysis 01,
and all three `fe3_monoanion` structures — across Fe distances of 1.94, 1.97 and 2.10 Å.

> **Correction (from analysis 04).** This was first written up as the term being a *fixed*
> reward rather than distance-graded. That was an overgeneralisation from a badly chosen
> sample: all the crystal poses sit at 1.94–2.10 Å, i.e. essentially identical geometry, so
> identical values say nothing about distance dependence. Redocked poses show intermediate
> values (e.g. 8SPD: −0.79, −0.53, 0.00 at various geometries), so −2.30 is a **saturation
> ceiling**, not a constant.
>
> More strikingly, the ordering can be inverted: for 8SPD the *near-native* pose
> (RMSD 1.73 Å) gets metal = **0.00** while a displaced pose at RMSD 3.96 Å gets **−0.79**.
> Glide's metal-term geometry window may sit at longer distances than a genuine dative
> bond, in which case it would actively mis-rank type-II poses rather than merely miss
> them. Being quantified across all 22 structures in analysis 04.

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

Keep the **chemically correct** heme for the baseline. Both activation routes
(`fe2_nh`, `fe3_monoanion`) require a chemically wrong porphyrin — a metalated porphyrin
has no pyrrole N–H and is a dianion.

Analysis 04 nonetheless repeats analyses 01 and 02 under **`fe2_nh`**, to test empirically
whether the wrong chemistry docks and scores better — understood going in that any
improvement would be right for the wrong reasons. `fe2_nh` is chosen over `fe3_monoanion`
because it perturbs electrostatics roughly three times less (mean `ecoul` shift −1.7 vs
−4.9; for 4K9X, exactly zero) and leaves Fe at +2 as in the baseline, making it a far
cleaner single-variable comparison.

Analyses 01 and 02 are deliberately **left intact** as the Fe(II) dianion baseline rather
than being re-run under a different representation; the sensitivity shown here is a more
useful result than any single set of numbers.

## Next

Test Glide's `metal_cons` constraint at grid generation. That imposes coordination
geometrically **without** misrepresenting the chemistry, and is the remaining legitimate
route to better pose recovery. Prediction recorded in advance: it should improve pose
RMSD while leaving the metal term at 0.00 — i.e. a geometric fix, not a scoring fix.
