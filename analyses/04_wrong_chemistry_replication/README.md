# 04 — Does the chemically *wrong* heme dock and score better?

**Question:** Analysis 03 found Glide's metal term is 0.00 for every chemically defensible
heme state, and −2.30 for two wrong ones. Does activating it actually improve docking?

**Status:** Complete.

**Answer: no.** Pose recovery and affinity correlation are both unchanged. The metal term
is real and behaves sensibly — but Glide's sampling almost never reaches the geometry
where it fires.

## Design

`fe2_nh` (Fe +2, one pyrrole N neutralised **and protonated**) reproduces "state B" from
analysis 01. Everything chemically correct and *independent of the metal* is held identical
to analyses 01/02 — same PrepWizard run, same deprotonated propionates, same chain and
ligand selection. **Only the porphyrin changes**, by one hydrogen and one formal charge.

`fe2_nh` was chosen over `fe3_monoanion` (the other activating state) because it perturbs
electrostatics far less and leaves Fe at +2 as in the baseline. That choice is vindicated
below: the residual `ecoul` contribution turned out to be −0.11 on average against a metal
contribution of −1.40, so this is close to a pure metal-term test.

Understood going in: **any improvement here would be right for the wrong reasons.** A
metalated porphyrin has no pyrrole N–H.

## Results

| | Baseline (`fe2_dianion`) | Wrong chem (`fe2_nh`) |
|---|---|---|
| Top-1 RMSD ≤ 2 Å | **4/22 (18 %)** | 3/22 (14 %) |
| Top-1 RMSD ≤ 1 Å | 2/22 | 3/22 |
| Best-of-10 ≤ 2 Å | 6/22 (27 %) | 6/22 (27 %) |
| Top-1 median RMSD | 4.85 Å | 4.54 Å |
| Affinity Spearman rho | −0.265 [−0.72, +0.28] | −0.293 [−0.74, +0.24] |
| Affinity Pearson r | −0.363 (p = 0.12) | −0.395 (p = 0.085) |

The correlation moves by 0.028 with near-identical confidence intervals — no meaningful
change. Pose recovery is flat to marginally worse.

### Attribution

Mean per-structure change in docking score: **−1.13**, decomposing as
metal **−1.40**, ecoul **−0.11**, evdw **+0.01** (component weights from analysis 03).
The score shift is essentially all metal term, as intended by the design.

### The metal term is distance-graded, and behaves sensibly

In-place crystal poses span the full range — 6UNJ −0.16, 6UNI −0.92, 6DAG −1.60,
8SPD/2V0M/4D7D/4K9X −2.30 — and the four non-coordinating ligands (6OO9, 7LXL, 3UA1,
4D6Z) all correctly receive **0.00**.

> **Correction.** Analysis 03 first reported the term as a *fixed* −2.30 reward. That came
> from three gate-test structures that all happened to sit at saturation (Fe 1.94–2.10 Å).
> With 22 structures the grading is obvious. Corrected there too.

### Why activating it doesn't help: the reward is unreachable

Metal term versus Fe distance across all 214 redocked poses:

| Fe distance | n poses | mean metal |
|---|---|---|
| 2.2–2.6 Å | 14 | **−0.289** |
| 2.6–3.0 Å | 35 | −0.029 |
| 3.0–3.5 Å | 72 | −0.015 |
| > 3.5 Å | 93 | ~−0.012 |

`rho = +0.271` (p = 6×10⁻⁵) — closer contacts do get more reward, the correct direction.
But **only 13 of 214 poses (6 %) receive any metal reward at all**, and no pose falls below
2.2 Å. Crystal ligands sit at 1.94–2.82 Å; docked poses essentially never get there.

The mechanism is circular: the metal reward would pull a ligand onto the iron, but it only
switches on once the ligand is already there. Glide's sampling does not reach into the
reward's basin, so activating the term changes almost nothing during docking — even though
it clearly rewards the crystal poses when scored in place.

> **Hypothesis refuted.** After the 8SPD smoke test I suggested the metal term might
> *favour displaced poses over near-native ones*, which would mean it actively mis-ranks
> type-II complexes. Across all 22 structures this is not true: near-native poses
> (RMSD ≤ 2 Å, n = 26) average metal −0.043 versus −0.033 for displaced poses (n = 188),
> Mann-Whitney **p = 0.71**. The 8SPD observation was a coincidence.

## Conclusion

Adopting the wrong chemistry buys nothing. Keep the chemically correct heme
(analyses 01/02 stand as the baseline).

Combined with analyses 01–03, the picture for **rigid-receptor Glide docking on CYP3A4**:

- the score does not track affinity even with a perfect pose (01);
- pose recovery is 18 % top-1, and sampling rarely produces a near-native pose at all (02);
- heme representation shifts scores by several GlideScore units and flips the metal term on
  single-atom changes that follow no chemical logic (03);
- and the one term that *should* anchor type-II inhibitors is unreachable from where
  sampling operates (04).

This supports a **co-folding-first strategy** for the challenge, consistent with the PXR
result where AlphaFold3 beat docking and Glide post-processing degraded good poses.
