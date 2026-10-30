# Leaderboard results

Scores returned by the organizers for each submitted attempt. Evaluated on a subset
(~10 of the 20 structures) while the challenge is open.

| # | Submission | Date | LDDT-PLI | BiSyRMSD (Å) | LDDT-LP | Coverage | Rank |
|---|---|---|---|---|---|---|---|
| 03 | [Boltz-2 defaults + heme, PrepWizard-minimised](03_boltz2_prepwizard/) | 2026-09-30 | **0.3921** | **4.1179** | 0.7962 | 1 | 5 / 9 |
| 05 | [AF3 + PrepWizard, PB-filtered, mean(PL,HL) iPTM](05_af3_prepwizard/) | 2026-10-01 | **0.4279** | **4.0653** | 0.6977 | 1 | — |
| 06 | [AF3 + PrepWizard, PL-iPTM only](06_af3_pliptm/) | 2026-10-01 | **0.4291** | **4.0647** | 0.6985 | 1 | — |
| 07 | [AF3 + 38FY/38FU cryo-EM templates](07_af3_templates/) | 2026-10-02 | 0.3853 | 4.2613 | 0.6762 | — | — |
| 08 | [AF3 + single 38FW holo template](08_af3_1template/) | 2026-10-02 | 0.4192 | **4.0543** | 0.6894 | — | — |
| 09 | [as 08 but PL-iPTM selection only](09_af3_1template_pliptm/) | 2026-10-04 | 0.3934 | 4.2269 | 0.6926 | — | — |
| 10 | [no template, HL-iPTM selection only](10_af3_hliptm/) | 2026-10-05 | 0.4221 | 4.0842 | 0.6953 | — | — |
| 11 | [06 with 3 structures from Boltz-2](11_af3_boltz_fe_swap/) | 2026-10-12 | 0.4241 | 4.0565 | 0.6972 | — | — |
| 12 | [ICM redocked poses in AF3's receptor](12_icm_redock_af3/) | 2026-10-13 | not submitted | — | — | — | — |
| 13 | [AF3 samples picked by Glide gscore](13_af3_glide_select/) | 2026-10-13 | 0.3730 | 4.5822 | 0.6968 | — | — |
| 14 | [AF3 samples picked by ICM RTCNN](14_af3_icm_rtcnn/) | 2026-10-14 | 0.3949 | 4.4057 | 0.7004 | — | — |
| 15 | [Boltz-2 tuned sampling, PrepWizard](15_boltz2_tuned/) | 2026-10-15 | 0.4063 | 4.0017 | 0.7983 | — | — |
| 16 | [Boltz-2 tuned, PL-iPTM selection](16_boltz2_pliptm/) | 2026-10-16 | **0.4229** | 3.9283 | **0.8060** | — | — |
| 17 | [Boltz-2 tuned, complex_iPDE selection](17_boltz2_ipde/) | 2026-10-16 | 0.4157 | 3.8021 | 0.7976 | — | — |
| 18 | [Boltz-2 tuned, Glide gscore selection](18_boltz2_glide/) | 2026-10-18 | 0.4174 | **3.7523** | 0.8047 | — | — |
| 19 | [AF3/Boltz hybrid, better Glide score](19_hybrid_glide/) | 2026-10-18 | 0.4331 | **3.6754** | 0.7823 | — | — |
| 20 | [best Glide over all AF3+Boltz samples](20_union_glide/) | 2026-10-19 | 0.3854 | 4.2275 | 0.7449 | — | — |
| 21 | [hybrid 19 with a 0.5 kcal/mol Glide margin](21_hybrid_margin/) | 2026-10-19 | 0.4174 | 3.7523 | 0.8047 | — | — |
| 22 | [19 with 5 Glide-ranked ICM poses](22_icm_fe_swap/) | 2026-10-21 | 0.4160 | 3.7855 | 0.7823 | — | — |
| 23 | [ESMFold2 + ColabFold MSA](23_esmfold2_default/) | 2026-10-22 | **0.4527** | 3.8551 | 0.7985 | — | — |
| 24 | [best Glide over AF3 + Boltz + ESMFold2](24_hybrid_three_way/) | 2026-10-22 | 0.4348 | 3.7261 | 0.7798 | — | — |
| 25 | [23 with 3 structures taken from 19](25_esm_outlier_fix/) | 2026-10-22 | 0.4465 | 3.9500 | 0.7969 | — | — |
| 26 | [ESMFold2 fold_max_accuracy endpoint](26_esmfold2_maxacc/) | 2026-10-24 | 0.4365 | 3.9024 | **0.8063** | — | — |
| 27 | [ESMFold2 default endpoint, NO MSA](27_esmfold2_nomsa/) | 2026-10-24 | 0.4384 | **3.8017** | 0.7908 | — | — |
| 28 | [best ICM RTCNN over the 3 ESMFold2 runs](28_rtcnn_best_of_three/) | 2026-10-25 | 0.4411 | 3.7722 | 0.7965 | — | — |
| 29 | [best Glide `inplace` over the 3 ESMFold2 runs](29_glide_inplace_best_of_three/) | 2026-10-26 | **0.4638** | 3.7277 | 0.8001 | — | — |
| 30 | [best Glide `inplace` over AF3 + Boltz-2 + ESMFold2](30_inplace_three_way/) | 2026-10-27 | 0.4414 | 3.7178 | 0.7808 | — | — |
| 31 | [ESMFold2 with `msa_max_depth` 2048](31_esmfold2_depth2048/) | 2026-10-27 | 0.4116 | 4.1680 | 0.7894 | — | — |
| 32 | [z-score consensus of Glide `inplace` + ICM RTCNN](32_zscore_consensus/) | 2026-10-28 | 0.4438 | 3.7662 | 0.8014 | — | — |
| 33 | [best Glide `inplace` over 2 ESMFold2 runs, 23+27](33_glide_best_of_two/) | 2026-10-28 | 0.4464 | 3.7705 | 0.7952 | — | — |
| 34 | [submission 29 with Glide-minimised ligand poses](34_glide_minimised/) | 2026-10-29 | 0.4637 | 3.7483 | 0.8001 | — | — |
| 35 | [best Glide `inplace` over 4 ESMFold2 runs, +set 31](35_glide_best_of_four/) | 2026-10-29 | 0.4512 | 3.8222 | 0.7975 | — | — |
| 36 | [best Glide `emodel` over 3 ESMFold2 runs](36_emodel_best_of_three/) | 2026-10-30 | pending | pending | pending | — | — |

Metric notes: **LDDT-PLI** is the primary score (protein-ligand interface, higher better).
**BiSyRMSD** is binding-site-superposed symmetry-corrected pose RMSD (lower better).
**LDDT-LP** scores the protein pocket itself. **Coverage** = fraction of structures that
scored at all; 1 means nothing was zeroed by OpenStructure or PoseBusters.

## Submission 30 — the cross-model rule, with the `mininplace` confound removed

Submission 24 applied the best-gscore rule across AF3, Boltz-2 and ESMFold2 and
scored 0.4348, which read as evidence that cross-model selection does not work.
Analysis 28 then found that every losing Glide selection in the project used
`mininplace`, which minimises the pose before rescoring and gives the worse
gscore on 38 of 60 structures. The only `inplace` selection, submission 29, is
the project's best at 0.4638.

Submission 30 is submission 24's exact rule rescored `inplace`. Analysis 30
supplied the 40 AF3/Boltz-2 scores; ESMFold2's already existed from analysis 28.
All 20 ligands have all three sets scored.

Picks: AF3 3, Boltz-2 8, ESMFold2 9. Switching mode moved Boltz-2 from 5 wins to
8 and AF3 from 5 to 3, and changed **5 of the 20 picks**. The set differs from
shipping ESMFold2 alone on **11 of 20**, so it is a real alternative rather than
a near-copy of submission 23.

OCNT-2395789 could not be selected at all under `mininplace` — neither AF3 nor
Boltz-2 returned a usable number. Under `inplace` all three score and ESMFold2
takes it at −7.903 against −4.645 and −5.879.

This is the harder version of the question submission 29 answered. Within one
model's three runs the candidates sit close together and the parents are equally
good. Across models they are far apart and 11 of the 20 come from parents that
ship at roughly 0.42–0.43, below ESMFold2's 0.4527 — so the rule has to beat a
real quality handicap. Above 0.4527 means the mode was the whole story behind
submission 24's loss; below it means Glide genuinely cannot compare poses built
on different receptor conformations, however they are scored.

### Result: 0.4414 / 3.7178 / 0.7808 — the mode was not the whole story

| | 24 `mininplace` | **30 `inplace`** | 23 ESMFold2 alone | 29 within-ESM |
|---|---|---|---|---|
| LDDT-PLI | 0.4348 | **0.4414** | 0.4527 | **0.4638** |
| BiSyRMSD | 3.7261 | **3.7178** | 3.8551 | 3.7277 |
| LDDT-LP | 0.7798 | **0.7808** | 0.7985 | 0.8001 |

Switching to `inplace` is worth **+0.0066** on the cross-model rule. That is a
real effect in the direction analysis 28 predicted, and it is the second-best
BiSyRMSD the project has posted (only submission 19's 3.6754 is lower). But it
closes barely a third of the gap: the set still lands **0.0113 below** simply
shipping ESMFold2 and **0.0224 below** selecting within ESMFold2's own three
runs.

So the two questions separate cleanly:

* **Within one model** (submission 29) `inplace` selection works, +0.0111 over
  the best parent.
* **Across models** it does not, even with the confound removed.

The LDDT-LP column says why. At 0.7808 submission 30 sits next to submission
24's 0.7798 and far below every ESM-only submission (0.7985–0.8063). LDDT-LP
scores the protein pocket, which is a property of the model that built it, not
of the pose Glide scored. Eleven of these 20 structures come from AF3 and
Boltz-2, whose pockets are worse, and the gscore has no way to see that — it
reads how well the ligand fits the receptor it is given, not how good that
receptor is. Picking the best-fitting ligand inside a worse pocket is still a
worse structure.

This also explains why within-model selection is the easier problem: the three
ESMFold2 runs share a pocket quality, so gscore is comparing like with like and
is free to be right about the pose. The practical rule for the rest of the
challenge: **select within the best model, never across models.**

## Submission 31 — MSA depth is a real dial that turns the wrong way

0.4414 was submission 30; this is worse. Doubling `msa_max_depth` from 1024 to
2048, a single changed field against submission 23, cost **0.0411 LDDT-PLI** —
the largest single-knob loss the project has recorded.

| | 23 (depth 1024) | 31 (depth 2048) | delta |
|---|---|---|---|
| LDDT-PLI | 0.4527 | 0.4116 | **-0.0411** |
| BiSyRMSD | 3.8551 | 4.1680 | +0.3129 |
| LDDT-LP | 0.7985 | 0.7894 | -0.0091 |

The metric split is the finding. LDDT-LP moves 0.009, so the protein pocket is
essentially unchanged; the loss is entirely in where the ligand sits. Analysis
31 measured the displacement in advance at a median 1.47 A against a 0.31 A
reproducibility floor, so we knew the poses had moved. We now know they moved
the wrong way, almost all of it error.

Two conclusions.

**1024 is the model's tuned optimum, not an arbitrary default.** Departing from
it degrades ligand placement while leaving the fold intact. The depths below
1024 are not worth folding on this evidence.

**The depth dial cannot generate a selection pool.** Submission 30 established
that selection needs candidates of comparable quality: picking across sets that
differ systematically in quality means the selector must first beat "always use
the best set", and it loses. Depth variants are exactly that case - 1024 is
systematically better by 0.041, four times the spread across the three existing
ESMFold2 configurations (0.4527 / 0.4384 / 0.4365, range 0.016). A pool holding
a 0.4116 set offers the selector mostly traps.

Which leaves submission 29's pool - the three real ESMFold2 configurations - as
the only diverse, comparable-quality candidate set this project has, and it
already produced the best result at 0.4638.

## Submission 32 — consensus averaged its inputs

0.4438 / 3.7662 / 0.8014. The z-score consensus of Glide `inplace` and ICM
RTCNN landed between its two inputs (0.4411 and 0.4638) and slightly below
their mean, dragged toward the weaker one.

So the near-decorrelation between the scorers - they agree on only 10 of 20
winners, Pearson +0.373 - bought nothing. Decorrelated errors are the textbook
precondition for consensus to beat both inputs, and the precondition held here;
the conclusion did not. Even the three ligands where consensus chose a set
neither scorer picked alone (OCNT-2315379, OCNT-2395457, OCNT-2535422) did not
pull it above the mean.

The reading: combining a good ranker with a worse one on the same pool
interpolates between them. It does not synthesise a better one. Glide alone
remains the rule to beat.

## Two templates (apo + holo) made things worse

| | 05 no templates | 07 with templates | change |
|---|---|---|---|
| LDDT-PLI | 0.4279 | **0.3853** | **−0.0426** |
| BiSyRMSD | 4.0653 | 4.2613 | +0.20 Å |
| LDDT-LP | 0.6977 | 0.6762 | −0.0215 |

**Every metric got worse**, and by more than the 05/06 selection difference — a real
effect, not noise.

### CORRECTION: this is NOT evidence that "templates don't work"

An earlier version of this file claimed PXR's template runs scored below its plain AF3 run
(0.4735 vs 0.5055). **That was a misreading** — 0.4735 was PXR's best *Boltz2* submission,
not a template run. The actual PXR numbers say the opposite:

| PXR entry | LDDT-PLI |
|---|---|
| 18 Plain AF3, ColabFold MSA | 0.5055 |
| 23 AF3 template search | 0.5087 |
| 24 AF3 + 9fzj template | 0.5104 |
| 36 AF3 + 7axe template, best-iPTM | 0.5132 |
| **26 AF3 + 9fzj template, best-iPTM** | **0.5285** — best AF3 result on PXR |

Single templates mostly helped there, and PXR's best AF3 submission used one.

**The distinguishing feature of submission 07 is that it used TWO templates, one apo
(38FU) and one holo (38FY).** Every successful PXR template run used exactly one. Apo and
holo pockets differ precisely in the F/G loop region, so supplying both plausibly gives
AF3 contradictory restraints on the one feature templates were meant to fix.

A secondary factor may be ligand mismatch: template ligands are 479–721 Da against 240–412
Da for the blinded set, max Tanimoto 0.218, 0/80 pairs >= 0.4.

**Next test: a single holo template, not an apo/holo pair.** — done, see below.

## One holo template (08): most of the damage was the pair, but templates still do not help

Submission 08 used a single holo template, **38FW (ketoconazole)**, chosen on binding mode
rather than map quality: 12 of the 20 blinded ligands coordinate Fe through an
imidazole/pyridine/triazole nitrogen at 2.06–2.32 Å (type II), exactly what ketoconazole
does, whereas vardenafil (38FY) is reported in the cryo-EM preprint as neither type I nor
type II. Selection used mean(PL,HL) iPTM, identical to 05 — so 05 vs 08 is the clean
apples-to-apples comparison.

| | 05 no template | 07 two templates | 08 one holo template |
|---|---|---|---|
| LDDT-PLI | **0.4279** | 0.3853 | 0.4192 |
| BiSyRMSD | 4.0653 | 4.2613 | **4.0543** |
| LDDT-LP | **0.6977** | 0.6762 | 0.6894 |

**The apo/holo pair was indeed most of the problem.** Going from two templates to one
holo recovered +0.0339 LDDT-PLI and 0.21 Å of BiSyRMSD — about 80 % of the regression.
That supports the conflicting-restraints explanation.

**But a single template still does not beat no template on CYP3A4.** 08 sits 0.0087
LDDT-PLI below 05 and 0.0083 LDDT-LP below it. The one thing it does win is BiSyRMSD —
4.0543 Å, the best value of any attempt so far — so the template nudges the pose
fractionally closer while slightly degrading interface contacts. Both differences are
small enough to be within run-to-run noise.

**Conclusion: templates are a dead end here, unlike PXR.** The likely reason is the one
already noted — the four cryo-EM ligands are 479–721 Da against 240–412 Da for the blinded
set, max Tanimoto 0.218, 0/80 pairs >= 0.4. PXR's useful templates were chemically close
to its targets; these are not. No further template variants are worth running.

**Best attempt remains 06** (no template, PL-iPTM selection, 0.4291).

## 09 vs 08: three structures, and the selection metric is not the story

Submission 09 is submission 08 with one line changed - select on `chain_pair_iptm[0][2]`
alone instead of the mean with the haem-ligand term. Exactly **3 of 20** structures differ
(OCNT-2313300, OCNT-2395328, OCNT-2395415). The score fell hard:

| | 08 mean(PL,HL) | 09 PL-iPTM only | change |
|---|---|---|---|
| LDDT-PLI | **0.4192** | 0.3934 | **-0.0258** |
| BiSyRMSD | **4.0543** | 4.2269 | +0.17 A |
| LDDT-LP | 0.6894 | 0.6926 | +0.0032 |

### What the size of that change implies

The leaderboard scores ~10 of the 20 complexes, so the damage is concentrated in however
many of those 3 structures are in the scored subset. Summed, it is about **-0.26 LDDT-PLI**
and **+1.7 A BiSyRMSD**. One scored structure would mean a single pose lost 0.26 LDDT-PLI
and moved 1.7 A; two would mean ~0.13 and ~0.9 A each. Either way this is a **pose flip**,
not a refinement - and LDDT-LP barely moved, confirming the protein is identical and only
the ligand placement changed.

### Two conclusions

**1. The 5 samples in a ligand's pool contain genuinely different binding modes.** In 2 of
the 3 swapped ligands the competing samples were tied in PL-iPTM to two decimals, yet the
poses differ by ~1.7 A. Within-ligand PL-iPTM spread is typically 0.01-0.03 (see the
best-of-5 scatter), so the metric is flat across exactly the axis that matters.

**2. The 05 -> 06 result was noise, and the conclusion drawn from it was wrong.** That A/B
moved 2 structures for +0.0012 LDDT-PLI and was logged as making 06 the best attempt, with
"use PL-iPTM going forward" as the takeaway. Against a -0.0258 swing from a 3-structure
change here, +0.0012 is indistinguishable from zero. **mean(PL,HL) and PL-only are not
meaningfully different as rules.** What moves the score is which pose is picked, and
neither metric controls that.

This is the same wall analyses 09 and 10 measured: no score predicts pose quality
(oracle 20/30 vs ~12/30 for any score), and no AF3 confidence metric predicts PoseBusters
validity (best 14/20 vs 11.4/20 random). Further selection-metric A/Bs are not worth a
submission slot - the next lever has to be a different *kind* of signal, or a bigger pool.

### Analysis 13 settled it: the regression is ONE ligand

Measuring the actual pose shift for each of the three swaps:

| ligand | 08 -> 09 pose shift |
|---|---|
| OCNT-2313300 | 0.50 A |
| OCNT-2395328 | 0.66 A |
| **OCNT-2395415** | **3.51 A** |

Two barely moved; one reoriented by 3.5 A. With ~10 complexes scored, one structure that
gets ~1.7 A worse raises mean BiSyRMSD by exactly the +0.17 A observed. **The selection
rule is not what moved the score.** See [analysis 13](../analyses/13_pose_divergence/).

### And the pools are wide

| run | within-ligand pairs | median | mean | max |
|---|---|---|---|---|
| 08, 38FW template | 117 | 1.81 A | 2.01 A | 6.42 A |
| 04, no template | 181 | 1.87 A | 2.25 A | 6.85 A |

AF3's five samples are genuinely different binding modes. Within one ligand's pool the
poses differ by ~2 A median while PL-iPTM varies by 0.01-0.03. **The confidence signal is
flat along the axis the poses vary on** - which is the mechanism behind every null result
in analyses 09, 10 and here.

Two consequences. There is real headroom: the pool contains better poses than the one we
ship. And **ensemble spread is a new signal orthogonal to iPTM** - tight pools
(OCNT-2395464 at 0.21 A) plausibly mark reliable predictions, wide ones (OCNT-2313429 at
5.48 A) mark cases where AF3 is undecided and any pick is near a coin flip.

## 03 vs 05

| | 03 Boltz-2 | 05 AF3 | change |
|---|---|---|---|
| LDDT-PLI | 0.3921 | **0.4279** | **+0.0358** |
| BiSyRMSD | 4.1179 | **4.0653** | −0.05 Å |
| LDDT-LP | 0.7962 | 0.6977 | **−0.0985** |
| Coverage | 1 | 1 | — |

**AF3 wins on the metric that counts.** LDDT-PLI is the primary score and it rose 9 %
relative, while BiSyRMSD barely moved — so the gain is in interface contacts, not in
getting the pose fundamentally more right. Both still sit ~4 Å, well outside the ~2 Å
usually called correct.

**LDDT-LP got worse, and that is the surprise.** AF3's protein pocket scores 0.0985 lower
than Boltz's. So AF3 is not simply the better model here: it builds a worse pocket and
still places the ligand better within it. Worth remembering before assuming AF3 dominates
— PXR saw AF3 ahead on both.

That split also means the two methods are not redundant. A hybrid - Boltz protein with
AF3's ligand placement, or consensus between them - is now motivated by data rather than
by hope.

## Reading the first result

**Coverage 1** — every structure scored. The PoseBusters work in submission 03 did its job;
nothing was zeroed on a technicality.

**LDDT-LP 0.7962 vs LDDT-PLI 0.3921** — the protein pocket is predicted well and the
*ligand* is what is wrong. That is the diagnostic that matters: effort spent on protein
quality is wasted, effort on ligand placement is not.

**BiSyRMSD 4.12 Å** — well outside the 2 Å usually called correct, and close to the
~4 Å best-pose values analyses 05/06 saw for structures where docking failed to sample a
near-native pose. Consistent with co-folding putting the ligand in roughly the right
cavity but the wrong orientation.

**Rank 5/9** — mid-field, from an essentially default Boltz-2 run. The gap to the leaders
is therefore unlikely to be closed by cleanup; it needs a different pose strategy.

## Candidate next attempts

Ordered by expected value, not by effort:

1. **`--diffusion_samples 5`** and select among them. Analysis 09 established the selection
   problem directly: oracle 20/30 versus ~12/30 delivered by any single score. More
   candidates plus a selector is the one lever measured to help.
2. **Fe-coordination restraint or filter.** 13/20 blinded ligands carry a classical
   Fe-ligating motif but Boltz coordinated only 11, and BiSyRMSD 4.12 Å suggests
   orientation errors a coordination constraint would penalise.
3. **AlphaFold3** as an independent co-folder, for consensus rather than replacement.
4. **Docking into the Boltz pocket** — but analysis 09 measured only ~27 % top-1 success
   for Glide/ICM on crystal receptors, so this is the weakest of the four.

## The three-way selection A/B is settled: the metric is worth ~0.007

Submissions 05, 06 and 10 are the *same* no-template predictions, differing only in which
confidence term picks the sample from each ligand's PoseBusters-passing pool.

| # | selection key | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|---|
| 06 | PL-iPTM alone | **0.4291** | **4.0647** | **0.6985** |
| 05 | mean(PL, HL) | 0.4279 | 4.0653 | 0.6977 |
| 10 | HL-iPTM alone | 0.4221 | 4.0842 | 0.6953 |

**Total spread across all three rules: 0.0070 LDDT-PLI and 0.02 A BiSyRMSD.**

Set that against submission 09, where a *single* ligand whose pose shifted 3.5 A moved
LDDT-PLI by 0.0258 - nearly four times the entire range of the metric choice. This is the
quantitative version of what analysis 13 predicted: **one unlucky pose outweighs the whole
selection-rule question.**

A weak secondary signal: the ordering is monotone in how much weight PL carries
(PL > mean > HL), so protein-ligand confidence does appear to hold slightly more
information than the haem-ligand term. But 0.007 across three rules is not a result to
build on.

One detail worth noting. HL-only differs from 05 in just 2 ligands, OCNT-2395342 and
OCNT-2395503 - and analysis 13 measured both as *high-spread* pools (3.53 A and 3.83 A
mean pairwise), exactly the risky kind. Yet the score moved only 0.0058. Either neither
fell in the ~10-complex scored subset, or both swaps happened to be benign. Compare
submission 09, where one such swap cost 0.026. **The variance is in which ligands the
organizers happen to score, not in the rule.**

### Where this leaves the strategy

Selection-metric variants are exhausted. Three rules, 0.007 apart, all well inside
single-structure noise. The remaining levers, in order of expected value:

1. **A bigger pool.** 180 PrepWizard-minimised, PoseBusters-passing samples already exist
   across submissions 04 and 08; every attempt so far has chosen within one run's five.
2. **A different kind of signal.** Analysis 13 identified ensemble spread - how much AF3
   agrees with itself - which is orthogonal to iPTM and currently unused.
3. **Accepting that ~12 of 20 pools are simply undecided** and that no selector will fix
   those; effort is better spent where the pool is tight.

## Where the error actually is (analyses 14 and 16)

Two negative results that each close off a direction.

**Glide scores say nothing about potency.** Analysis 14 tested gscore and evdw against the
15 blinded ligands that carry a CYP3A4 pIC50, under two pose-selection rules, plus fitted
linear and quadratic models on (evdw, ecoul). Every leave-one-out Q² is negative - the
models predict a held-out ligand worse than guessing the mean. The quadratic reaches
R² = 0.551 in-sample and Q² = -0.100 cross-validated, which is 6 parameters memorising 15
points. Combined with analyses 01 and 08, ~30 correlations have now been tested against
these same pIC50 values over ~2 log units with no survivor. Closed.

**Neither model misplaces the pocket.** Analysis 16 measured the haem Fe to ligand distance
in submissions 06 and 03. Nothing in either lies beyond 4.7 A of the iron; the distribution
is two clean clusters, 12 type II binders coordinating at 2.05-2.31 A and 8 type I binders
at 3.59-4.66 A. AF3 coordinates MORE ligands than Boltz (12 vs 11), not fewer.

That second result is the more useful one, because it **localises the remaining error**:
both models place the ligand correctly relative to the haem, so the ~4 A BiSyRMSD is
orientation and conformation *within* the right pocket, not gross mislocation. Effort aimed
at pocket identification is wasted. It also removes the stated motivation for a
Boltz-protein / AF3-ligand graft, since the LDDT-LP gap is backbone quality, a separate
thing from where the ligand sits.

One genuine chemical disagreement survives: **OCNT-2535422**, where AF3 coordinates Fe at
2.28 A but the ligand has no imidazole, pyridine or triazole nitrogen - only a sulfonamide
NH and an indoline NH, both poor donors. It should not be type II, so Boltz's
non-coordinating 3.80 A is likely right. Best-justified single swap available; deliberately
held out of submission 11 to keep that a clean test.

### Candidate rules not submitted (analysis 15)

| rule | differs from 06 | median shift | >2 A |
|---|---|---|---|
| best Glide gscore | 16/20 | 2.23 A | 8 |
| consensus centroid, 04 pool | 14/20 | 1.17 A | 4 |
| consensus centroid, 04+08 pooled | 15/20 | 1.68 A | 5 |

Glide is the wildest bet and the worst-motivated - analysis 10 found it selects
PoseBusters-valid poses below random, and analysis 14 found its scores anti-correlate with
potency. High variance with no evidence of skill.

## Submission 11 result: the Fe-distance heuristic does not work

| | 06 | 11 (3 structures from Boltz) | change |
|---|---|---|---|
| LDDT-PLI | **0.4291** | 0.4241 | -0.0050 |
| BiSyRMSD | 4.0647 | **4.0565** | -0.0082 (better) |
| LDDT-LP | **0.6985** | 0.6972 | -0.0013 |

**A prediction recorded here was wrong.** Submission 11's README argued that taking 3 of
20 structures from Boltz should raise LDDT-LP by ~+0.015 "fairly reliably", because
Boltz's protein scored 0.7962 against 06's 0.6985. It fell 0.0013 instead. The 0.098 gap
is an average over 20 structures, not a per-structure guarantee - whichever swapped
structures were scored did not carry it.

**The deltas are non-zero, so at least one swap IS in the ~10-complex scored subset.** So
this is a real negative result, not the uninformative no-change case that was flagged as
likely. Scaling by 10: one scored swap lost ~0.05 LDDT-PLI while gaining ~0.08 A
BiSyRMSD.

That split is the interesting part and it agrees with analysis 16. The Boltz pose sits
*closer to the truth in position* while making *worse interface contacts*. Both models
already have the right pocket, so moving a ligand nearer the iron changes its contacts
without necessarily improving them - "closer to Fe" is not a proxy for "better pose".

Still untested and better justified: OCNT-2535422, held out of 11 deliberately. That swap
rests on chemistry (no imidazole / pyridine / triazole N, so it should not be a type II
binder and AF3's 2.28 A coordination looks forced) rather than on this heuristic.

## Docking into a predicted receptor (analysis 17, submission 12)

All 20 ligands redocked with ICM into both predicted receptors, ligands built from SMILES
so no pose information leaks in. Settings match analysis 06: confs=10, effort=3., seed=42.

| | Boltz (03) receptor | AF3 (06) receptor |
|---|---|---|
| top-1 RMSD to co-folded, median | 6.38 A | **5.04 A** |
| best-of-10 RMSD, median | 2.36 A | **1.69 A** |
| top-1 within 2 A of co-folded | 0/20 | **2/20** |
| best-of-10 within 2 A | 6/20 | **12/20** |
| top-1 coordinating Fe | 0/20 | **4/20** |

**Backbone accuracy and dockability come apart.** Boltz's pocket scores 0.098 higher on
LDDT-LP, yet AF3's receptor docks better on every measure. A prediction made before this
ran - that Boltz's receptor would be the better target *because* of its LDDT-LP - was
wrong. LDDT-LP scores the pocket backbone; it says nothing about whether the side chains
leave a dockable cavity.

**Docking loses the Fe coordination that co-folding gets for free.** 12/20 co-folded poses
coordinate the iron; ICM keeps 4/20, or 0/20 in the Boltz receptor. Analysis 12 already
showed Glide metal constraints cannot restore it either.

**Sampling is fine, ranking is not.** Best-of-10 median 1.69 A with 12/20 pools holding a
near-co-folded pose, against a top-1 median of 5.04 A. Same split as analysis 09 on
crystal receptors, and the same wall as analyses 10 and 13.

Submission 12 ships the AF3-receptor top-1 poses as a deliberate floor measurement
(20/20 PoseBusters). Expected well below 06.

## Settled: AF3's own confidence beats every external scoring function

Five submissions now select from the SAME 100 PrepWizard-minimised AF3 samples, differing
only in the ranking rule.

| selector | LDDT-PLI | BiSyRMSD | submission |
|---|---|---|---|
| PL-iPTM | **0.4291** | **4.0647** | 06 |
| mean(PL, HL) | 0.4279 | 4.0653 | 05 |
| HL-iPTM | 0.4221 | 4.0842 | 10 |
| ICM RTCNN | 0.3949 | 4.4057 | 14 |
| Glide gscore | 0.3730 | 4.5822 | 13 |

**The three confidence rules cluster within 0.0070; both external scores sit 0.03-0.06
below all of them.** A physics or ML score computed on the finished pose carries less
information about which sample is right than the model's own confidence did while
generating it. That is the single clearest result of the project.

**The analysis 09 ranking transferred.** That analysis measured within-ligand Spearman
against RMSD on crystal receptors - RTCNN +0.334, Glide lipo +0.115, gscore +0.086 - and
predicted RTCNN would be the better of the two selectors. It is, by 0.022. The relative
order held moving from crystal to co-folded receptors; the absolute level was simply far
too low to beat AF3.

**The swaps are lateral.** Both external rules change 16 of 20 structures while leaving
the Fe chemistry essentially untouched: of the 15 ligands with a nitrogen donor, 06
coordinates 12, RTCNN 11, Glide 11, with median Fe distances 2.25 / 2.34 / 2.32 A. So a
rule can reshuffle most of the set, keep every binding mode intact, and still lose
0.03-0.06 LDDT-PLI - the damage is in orientation within the pocket, invisible to any
coordination check.

### Docking into a predicted receptor is worse still

Analysis 17 (ICM) and 18 (Glide) redocked all 20 into the predicted receptors. Of the 15
Fe-capable ligands, ICM's top-1 keeps 4 coordinations and Glide's keeps **0**, against
AF3's 12. Submission 12 (ICM redock) was built and validated but not submitted once
submission 13 showed selection alone already costs 0.056; the Glide redock was not
packaged at all.

## Why Boltz beats AF3 on LDDT-LP (analysis 21)

Scored every predicted pocket against 38 ligand-bound CYP3A4 crystals at <=2.5 A, on a
fixed 10 A haem shell, with OpenStructure's lDDT.

| | max lDDT | median vs ensemble |
|---|---|---|
| Boltz-2 (03) | **0.9726** | **0.9513** |
| AF3 (06) | 0.9650 | 0.8955 |

Boltz wins on 19 of 20 ligands, and **the gap is in the median (0.056), not the max
(0.008)**. Both models match some crystal closely; only Boltz resembles the whole
ensemble.

The best-match counts make the mechanism explicit: **AF3 lands on 6DAA for 15 of 20
ligands** - the same pocket conformation almost regardless of ligand - while Boltz
spreads across 7UFB, 4D75, 6BCZ, 6DAJ and 6OOB. AF3 produces a near ligand-independent
pocket; Boltz adapts.

CYP3A4's site is highly plastic and the blinded references will span several states, so a
model locked to one state scores well on the few that match and poorly on the rest. That
is a plausible mechanism for the 0.098 gap.

It also resolves what looked like a contradiction: analysis 17 found AF3's receptor is the
BETTER docking target despite the worse LDDT-LP. A single well-formed induced-fit pocket
can be easy to dock into and still be the wrong conformation for most references.

## Tuning Boltz's sampling is the first change in a long while that moved anything

Submission 15 is submission 03 with three flags changed - `--recycling_steps 5`
(default 3), `--sampling_steps 300` (default 200), `--diffusion_samples 5` (default 1) -
and nothing else.

| | 03 default | 15 tuned | change |
|---|---|---|---|
| LDDT-PLI | 0.3921 | 0.4063 | **+0.0142** |
| BiSyRMSD | 4.1179 | **4.0017** | **-0.116 A** |
| LDDT-LP | 0.7962 | **0.7983** | +0.0021 |

**Best of any submission on BiSyRMSD and LDDT-LP.** Submission 06 still leads LDDT-PLI
at 0.4291.

**The PXR transfer held.** On PXR the same change gave +0.0182 (run 8 -> run 9); here it
gives +0.0142. PXR also showed 6/400 gives back a third of the gain, so 5/300 is the
setting to keep, not a direction to push further.

Worth the perspective: +0.0142 is **double the 0.0070** that separated every
selection-metric variant we tried (05, 06, 10), and it came from generation settings
rather than post-hoc selection. After analyses 09/10/20 and submissions 13/14 all failed
to improve selection, this says the remaining headroom is in how the poses are GENERATED.

`diffusion_samples 5` gave an ensemble we cannot yet exploit: Boltz's own
`confidence_score` varies by only 0.001-0.007 across the five, as flat as AF3's iPTM was.
The four unused samples per ligand are kept on disk. The untried alternative is Boltz's
`pair_chains_iptm[0][2]` - its protein-ligand iPTM, the analogue of the metric that won
for AF3 - which reads ~0.70 where confidence_score reads ~0.95, so it has far more range.

## Submission 16: the first selection rule that worked, and why

Same 5-sample Boltz ensemble as submission 15, different sample chosen.

| | 15 Boltz's own pick | 16 PL-iPTM | change |
|---|---|---|---|
| LDDT-PLI | 0.4063 | **0.4229** | **+0.0166** |
| BiSyRMSD | 4.0017 | **3.9283** | -0.073 A |
| LDDT-LP | 0.7983 | **0.8060** | +0.0077 |

**First submission under 4 A BiSyRMSD and the best LDDT-LP of any attempt** (0.806 against
AF3's 0.699). Submission 06 still leads LDDT-PLI by 0.0062.

Every earlier selection attempt failed: the three AF3 confidence rules spanned 0.0070
(noise), Glide gscore lost 0.056, ICM RTCNN lost 0.034. **The difference here is dynamic
range.** Median spread across a ligand's 5 samples: Boltz `confidence_score` 0.0038,
Boltz `pair_chains_iptm[0][2]` **0.079** - 20x more. A selector can only help if it
varies across the candidates it is ranking.

Note the asymmetry of `pair_chains_iptm`: `[0][2]` = 0.698 while `[2][0]` = 0.956 for the
same interface, because ipTM's distance cutoff is set by the size of the set being scored
- ~27 ligand atoms (strict, discriminating) versus ~500 protein residues (loose, pinned
near the ceiling). `[0][2]` is the usable direction and the analogue of AF3's
`chain_pair_iptm[0][2]`.

### Cumulative, from default Boltz

| step | LDDT-PLI | gain |
|---|---|---|
| 03 defaults | 0.3921 | - |
| 15 tuned sampling | 0.4063 | +0.0142 |
| 16 + PL-iPTM selection | **0.4229** | +0.0166 |
| | | **+0.0308** |

**Generation and selection are complementary.** More diverse candidates are only useful
with a selector that can tell them apart, and a selector is only useful when the
candidates differ. That reframes the earlier null results: AF3's selection experiments
failed partly because AF3's five samples were too similar in confidence to choose between.

## Analysis 23: Glide cannot select Boltz samples either

All 100 Boltz samples scored with Glide `mininplace` (PrepWizard -> own-conformation grid
-> score), 0 failures. Agreement with the PL-iPTM pick that produced submission 16:

| | agrees with PL-iPTM |
|---|---|
| Glide `gscore` | 4/20 |
| Glide `emodel` | 5/20 |
| chance | ~4/20 |

Glide is not flat - median within-ligand spread 1.83 kcal/mol against PL-iPTM's 0.079 -
so it discriminates confidently on something uncorrelated with pose quality.

Three independent lines now agree: analysis 09 measured within-ligand rho = +0.086 on
crystal receptors, submission 13 lost 0.056 selecting the AF3 ensemble with it, and here
it agrees with the working rule at chance. The failure belongs to the scoring function,
not to the candidate pool, which is why changing the ensemble changed nothing.

**No further submission slots on Glide or ICM selection.** What distinguishes PL-iPTM is
not dynamic range - Glide has more - but that it is the model's own estimate of how well
it placed the ligand, formed with access to generation, rather than a physics score
re-judging a finished pose.

## Submission 18 overturns the recommendation against Glide selection

Analysis 23 found Glide's picks agree with PL-iPTM at chance, and submission 13 had lost
0.056 using Glide on the AF3 ensemble. On that basis the advice was to skip it. **The
advice was wrong.**

| Boltz selection rule | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|
| 15 Boltz `confidence_score` | 0.4063 | 4.0017 | 0.7983 |
| 16 PL-iPTM | **0.4229** | 3.9283 | **0.8060** |
| 17 complex_iPDE | 0.4157 | 3.8021 | 0.7976 |
| 18 Glide gscore | 0.4174 | **3.7523** | 0.8047 |

Two reasoning errors worth keeping:

1. **A score's failure on one ensemble does not transfer to another.** Glide cost 0.056 on
   AF3's samples and 0.0055 on Boltz's, while winning BiSyRMSD outright.
2. **Agreement with a working rule is not a proxy for performance.** Two selectors can
   disagree on which sample to take and both land on good poses. Only the leaderboard
   settles it.

**BiSyRMSD improves monotonically across the four rules while LDDT-PLI peaks at 16**, so
the two metrics reward different things. The Fe analysis shows why: Glide's picks differ
from iPDE's on 16 of 20 samples but move the ligand >0.05 A for only 6, all of them
non-coordinating type I ligands that Glide pulls closer to the haem (e.g. OCNT-2395328
4.31 -> 3.01 A). Coordination is identical at 11/15, median 2.29 A.

This qualifies submission 11, where moving three ligands closer to the iron by hand HURT.
A physical interaction model that happens to prefer closer type I poses helps; a bare
distance heuristic does not.

## Submission 19: can Glide choose between two whole models?

Per ligand, whichever of submission 06 (AF3 + PL-iPTM) or 18 (Boltz + Glide) Glide scores
better. 14 from 18, 6 from 06. 20/20 structural and PoseBusters.

**The parents agree on almost nothing**: median ligand RMSD 4.15 A, max 7.71 A, 17 of 20
beyond 2 A - yet both score ~0.42 LDDT-PLI. Two coordinating ligands show the pattern
sharply: OCNT-2313300 is 6.21 A apart with both at 2.13 A from the iron, OCNT-2535519 is
4.90 A apart with both at ~2.10 A. Same metal contact, molecule pivoted around the
coordinating nitrogen.

Since they cannot both be right, roughly half of each submission's apparent accuracy comes
from different structures. If Glide can identify which half, the hybrid beats both parents;
if it lands between them, it cannot.

Caveats recorded in the submission README: the two gscores are computed on DIFFERENT
receptors (AF3's protein vs Boltz's, each with its own grid), and the pose disagreement is
uncorrelated with the score difference (Pearson +0.110, Spearman +0.018, n=19).

## Submission 19 result: Glide CAN choose between models

**0.4331 LDDT-PLI / 3.6754 A BiSyRMSD / 0.7823 LDDT-LP - best of the project on the
primary metric and on BiSyRMSD**, beating both parents (06: 0.4291, 18: 0.4174).

| submission | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|
| 06 AF3 + PL-iPTM | 0.4291 | 4.0647 | 0.6985 |
| 16 Boltz + PL-iPTM | 0.4229 | 3.9283 | **0.8060** |
| 17 Boltz + iPDE | 0.4157 | 3.8021 | 0.7976 |
| 18 Boltz + Glide | 0.4174 | 3.7523 | 0.8047 |
| **19 hybrid** | **0.4331** | **3.6754** | 0.7823 |

**The resolution to the Glide puzzle is scale.** Within one model's five samples the
candidates differ by ~1-2 A, Glide's score differences are noise, and it picks at chance
(analysis 23) or loses outright (submission 13, -0.056). Between AF3 and Boltz the poses
differ by a median 4.15 A - a gap large enough for a physics score to resolve - and it
beats both parents. Selection ability is not a fixed property of a score; it depends on
how far apart the candidates are.

**The LDDT-LP drop is the pocket, not the ligands.** LDDT-LP scores the protein binding
site. Mixing 14 Boltz structures (0.8047) with 6 AF3 ones (0.6985) predicts
(14*0.8047 + 6*0.6985)/20 = 0.773 against the observed 0.7823 - the whole deficit is the
price of importing AF3's weaker pocket for its ligand poses, and says nothing about
ligand placement.

## Submission 20: one Glide rule over the whole pool

Best gscore across all 94 scored AF3 samples and all 100 Boltz samples. **12 AF3 /
8 Boltz**, changing 10 of 20 structures from submission 19 and flipping the balance from
19's 6/14 - the union sees AF3 samples submission 06 never shipped.

It is a clean test of the scale explanation for Glide's behaviour:

* **If 20 beats 19**, the scale story is wrong and Glide is simply a good selector given
  enough candidates.
* **If 20 lands below 19**, the scale story holds: Glide's skill is specifically in
  comparing whole models (~4 A apart), not in ranking within an ensemble (1-2 A apart),
  and a global rule reverts to the losing within-ensemble behaviour.

The risk is concrete: submission 13 applied this exact rule within AF3 alone and scored
0.3730, the worst AF3-based attempt, and 12 of these 20 picks are AF3 samples chosen on
gscore.

## Submission 20 result: Glide's skill is cross-model only

0.3854 / 4.2275 / 0.7449 - a 0.048 LDDT-PLI loss against submission 19, finishing just
above submission 13 (0.3730), the within-AF3 Glide rule.

| | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|
| 19 pairwise, cross-model | **0.4331** | **3.6754** | 0.7823 |
| 20 union, global | 0.3854 | 4.2275 | 0.7449 |
| 13 Glide within AF3 alone | 0.3730 | 4.5822 | 0.6968 |

**The scale explanation is confirmed.** Between AF3 and Boltz, candidates differ by a
median 4.15 A and Glide resolves them - submission 19 beat both its parents. Within one
model's ensemble they differ by 1-2 A, and Glide picks at chance (analysis 23) or loses
(13, 20). A global rule destroys the advantage, because 12 of 20 picks became AF3 samples
ranked on gscore - submission 13's rule exactly.

**Practical rule: give Glide a choice between models, never between samples.**

### LDDT-LP for a hybrid is predictable in advance

(12 x 0.6985 + 8 x 0.8047) / 20 = 0.741 against 0.7449 observed. Submission 19 predicted
0.773 against 0.7823. Twice accurate to ~0.004, so LDDT-LP is essentially a weighted
average of the parent pockets and can be computed before submitting.

## Submission 21: is a 0.07 kcal/mol Glide margin worth an AF3 pocket?

Submission 19 with AF3 kept only where it wins Glide by more than 0.5 kcal/mol:
17 Boltz / 3 AF3, changing 3 structures (OCNT-2395457, -2395503, -2395712, whose margins
were 0.31, 0.08 and 0.07).

**Prediction on record.** The mixing model has been accurate to ~0.004 twice, both times
running ~0.009 above the raw weighted average:

| submission | mix | predicted | observed |
|---|---|---|---|
| 19 | 14 Boltz / 6 AF3 | 0.773 | 0.7823 |
| 20 | 8 Boltz / 12 AF3 | 0.741 | 0.7449 |
| 21 | 17 Boltz / 3 AF3 | **0.789** | expect ~0.798 |

If LDDT-PLI holds near 19's 0.4331 this is a strict improvement of ~+0.016 LDDT-LP. If it
falls toward submission 18's 0.4174, the three reverted AF3 poses were earning their keep
and the margin threshold is the wrong idea.

## Submission 21 scored identically to 18 - and that is informative

21 (17 Boltz / 3 AF3) returned 0.4174 / 3.7523 / 0.8047, **exactly** submission 18's
numbers. A file-by-file diff shows 17 of 20 structures are byte-identical and the three
that differ are OCNT-2312966, OCNT-2315379 and OCNT-2535519 - precisely the AF3
structures 21 kept.

**So none of those three ligands is in the ~10-complex scored subset.** Any future variant
differing only on them will return an identical score and waste a slot.

It also corrects the mixing model: the LDDT-LP prediction of ~0.798 failed because the
model assumed all 20 structures contribute, when only the scored subset does. Submissions
19 and 20 happened to have enough AF3 structures in that subset for the average to look
right.

## Submission 23: ESMFold2 as a third co-folding method

| | coordinating (of 15) | median Fe |
|---|---|---|
| 06 AF3 | 12 | 2.25 A |
| 19 hybrid (best) | 11 | 2.29 A |
| 16 Boltz | 11 | 2.29 A |
| **23 ESMFold2** | **12** | **2.16 A** |

Competitive Fe geometry, and operationally very different: the compute is **remote**
(Biohub API, no SLURM job), ~28 s per ligand against ~3 h for the tuned Boltz run, and
the chain layout needs no fixing (A protein / B HEM at 43 atoms / C LIG).

**But it is deterministic** - two identical calls gave byte-identical mmCIF - so there is
no ensemble to select from, and `pair_chains_iptm` returns `None`, so the PL-iPTM rule
that produced submission 16 is unavailable. Its value here is as a third candidate for
cross-model hybrids, which is the regime where Glide selection demonstrably works
(submission 19 beat both parents).


## Submission 24: cross-model Glide selection is not a real capability

| | 19 AF3+Boltz | 24 +ESMFold2 | 23 pure ESMFold2 |
|---|---|---|---|
| LDDT-PLI | 0.4331 | 0.4348 | **0.4527** |
| BiSyRMSD | **3.6754** | 3.7261 | 3.8551 |
| LDDT-LP | 0.7823 | 0.7798 | **0.7985** |

Adding a third co-folding model and letting Glide choose moved LDDT-PLI by **+0.0017**
over submission 19 - noise - while costing **0.0179** against simply shipping ESMFold2
unaltered. BiSyRMSD and LDDT-LP both got slightly worse.

This was the test submission 24 was built for, and it fails cleanly. The argument stated
when it was built was that a genuinely capable cross-model ranker should have preferred
ESMFold2 far more often than 5 of 20, given ESMFold2 outscores both other models on its
own. It preferred it 5 times, and the result is a hybrid that is 0.018 worse than its
best parent.

The measurement is real rather than invisible: 24 differs from 19 on exactly five
ligands, and the score did move, so at least one of the five is inside the ~10-structure
scored subset. The five were not marginal either - pocket-superposed RMSD between 19 and
24 is 3.6-4.2 A for four of them (analysis 22), squarely in the range where submission
19 appeared to show Glide discriminating.

**Conclusion.** Submission 19 beating both its parents by 0.004-0.016 is best read as
noise on a ~10-structure subset, not as evidence of a selection rule. No scoring
function in this project - Glide gscore, ICM RTCNN, or a co-folding model's confidence -
has been shown to choose structures across co-folding programs. What has repeatedly
worked is trying a better model: AF3 to Boltz-2 to ESMFold2.


## Submissions 22 and 25: post-processing a good submission has never worked

| | parent | result | delta |
|---|---|---|---|
| 22 = 19 + 5 Glide-ranked ICM poses | 0.4331 | 0.4160 | **-0.0171** |
| 25 = 23 + 3 structures from 19 | 0.4527 | 0.4465 | **-0.0062** |

Submission 22 is the Fe-distance heuristic's third defeat, after submissions 11 (-0.005)
and 20. It also imported docked poses, and analysis 17 had already measured ICM redocking
at 4 of 15 ligands coordinating against co-folding's 11-12, so the structures were worse
on the one chemical criterion available before they were ever shipped.

Submission 25 is the more informative loss. Two of its three swaps corrected values that
are not physically possible - OCNT-2315379 at 6.89 A with a nitrogen donor, OCNT-2395460
at 1.79 A with none - and the third, OCNT-2395328, was included on a 2.7 kcal/mol Glide
gap over the objection that it repeats the Fe heuristic. The submission as a whole lost
0.0062, and which of the three is responsible cannot be read off this number.

Submission 26's score has since undercut the reading given here. Its max-accuracy
structures corrected BOTH impossible distances independently (6.89 -> 4.02 A,
1.79 -> 3.82 A) and moved OCNT-2395328 the other way, which looked like confirmation
that the first two swaps were right. But 26 scored 0.4365, losing 0.0162 to submission
23. Two unrelated routes to fixing those two distances - swapping structures in, and
re-folding them - both lost. What 26 corroborated was that two methods AGREE those
values are unusual, not that correcting them helps the score.

**Running count of attempts to improve a good submission by rearranging or replacing its
structures: 11, 20, 21, 22, 24, 25 - six attempts, six losses.** The only thing that has
ever gained is running a better model. Submission 23 at 0.4527, pure ESMFold2 with no
selection step of any kind, remains the best result of the project.


## Submission 26: the best pocket of the project, and a worse interface

| | 23 default endpoint | 26 max-accuracy | delta |
|---|---|---|---|
| LDDT-PLI | **0.4527** | 0.4365 | **-0.0162** |
| BiSyRMSD | **3.8551** | 3.9024 | +0.0473 |
| LDDT-LP | 0.7985 | **0.8063** | **+0.0078** |

0.8063 is the best LDDT-LP of the project, just past submission 16's 0.8060. LDDT-LP
scores the protein pocket itself and is largely independent of where the ligand sits, so
max-accuracy builds the better pocket and then places the ligand in it worse. That is a
coherent result for an endpoint that generates its own MSA and tunes for structural
accuracy, and it is the clearest separation of the two metrics we have seen.

It also partly answers whether our ColabFold MSA was earning its keep. Submission 23 used
it and beat 26, which built its own, by 0.0162 on the primary metric - though the two
differ in endpoint and model weights as well, so the MSA is not isolated.

Seven attempts now, seven losses, and this one changed the structures themselves rather
than rearranging them. Submission 23 stands at 0.4527.


## Submission 27: what the MSA is worth, measured cleanly

This is the only properly controlled comparison in the project. Submission 27 differs
from submission 23 in exactly one input - `msa=None` - with the same `fold_all_atom`
endpoint, the same esmfold2-2026-05 weights, the same sequence and chain layout, and the
same single non-ESMFold2 structure at OCNT-2313038.

| | 23, ColabFold MSA | 27, no MSA | MSA is worth |
|---|---|---|---|
| LDDT-PLI | **0.4527** | 0.4384 | **+0.0143** |
| BiSyRMSD | 3.8551 | **3.8017** | **-0.0534** (the MSA makes this WORSE) |
| LDDT-LP | **0.7985** | 0.7908 | +0.0077 |

So the 14,301-sequence ColabFold MSA earns its keep on the primary metric and on the
pocket, and costs a little on pose RMSD. Worth keeping, but it is a 0.014 effect, not a
0.05 one.

Two corrections to things recorded earlier fall out of this.

First, the claim after submission 27 was folded that the MSA "does essentially nothing"
rested on Fe distances, which are nearly identical between the two sets. Analysis 27
measured pocket-superposed ligand RMSD at a median of 1.53 A, and the score now confirms
those pose differences are real and cost 0.0143.

Second, submission 26's 0.0162 deficit against 23 was attributed to an unresolved mix of
endpoint, weights and self-built MSA. Most of it is the MSA: no-MSA default scores 0.4384
and max-accuracy scores 0.4365, so the endpoint and weight change is worth only about
0.002 beyond losing our MSA.

All three ESMFold2 runs land between 0.4365 and 0.4527, a span of 0.016, while the choice
between AF3, Boltz-2 and ESMFold2 spans 0.4174 to 0.4527. Which model is used still
dominates how it is configured.


## Submission 28: RTCNN across runs, and do PLI and RMSD trade off?

0.4411 / 3.7722 / 0.7965. Against its three parents:

| | 23 | 26 | 27 | **28 RTCNN pick** |
|---|---|---|---|---|
| LDDT-PLI | **0.4527** | 0.4365 | 0.4384 | 0.4411 |
| BiSyRMSD | 3.8551 | 3.9024 | 3.8017 | **3.7722** |
| LDDT-LP | 0.7985 | **0.8063** | 0.7908 | 0.7965 |

RTCNN beats two of three parents on the primary metric and all three on BiSyRMSD - the
best pose RMSD of any ESMFold2 submission, and second best in the project behind
submission 19's 3.6754. It does not beat submission 23, so cross-run RTCNN selection is
not a win, but neither is it the flat failure that within-ensemble selection was
(submission 14 lost 0.034 against its baseline; this loses 0.012 and gains on RMSD).

### Do LDDT-PLI and BiSyRMSD trade off?

Asked because the ESMFold2 submissions look like they do. Spearman between the two,
over every scored submission and within slices:

| scope | n | Spearman PLI vs BiSyRMSD |
|---|---|---|
| all scored submissions | 23 | **-0.596** - they agree |
| PLI >= 0.42 | 12 | -0.545 - agree |
| PLI < 0.42 | 12 | -0.811 - agree strongly |
| ESMFold2 family only | 5 | **+0.100** - mild tradeoff |

(Negative means they agree: higher LDDT-PLI with lower RMSD.)

So globally there is no tradeoff - the metrics track each other. The tradeoff is real but
confined to the ESMFold2 family, where the ranking inverts: 23 has the best PLI with
mediocre RMSD, 28 the best RMSD with middling PLI.

A plausible reading: across the project, variation is dominated by WHICH MODEL produced
the structures, and a better model lifts both metrics together. Within the ESMFold2
family the model is fixed and only the poses are reshuffled or re-folded, leaving the
genuine geometric tension - a pose can sit nearer the true position overall while making
slightly worse interface contacts. At n=5, +0.100 is weak evidence and should not carry
much weight.


## Public standing, 2026-10-26

| rank | user | LDDT-PLI | BiSyRMSD | LDDT-LP |
|---|---|---|---|---|
| 1 | nova | 0.4906 | 3.6471 | 0.7679 |
| 2 | z | 0.4770 | 3.3255 | 0.8051 |
| 3 | CYP-Random | 0.4750 | 3.2924 | 0.8085 |
| 4 | TCB | 0.4578 | 3.7838 | 0.8111 |
| 5 | Fe-Nomenal | 0.4512 | 3.4394 | 0.8169 |
| **9** | **TangerineTrees (us)** | **0.4411** | 3.7722 | 0.7965 |

**THE TABLE RANKS THE LATEST SUBMISSION, NOT THE BEST.** The 0.4411 shown is submission
28, the most recent upload. When submission 23 was the latest the rank was 5th at 0.4527,
and 0.4527 would still place 5th today, ahead of Fe-Nomenal.

So every experiment run from here lowers the public standing while it is the newest
entry. That is acceptable during development, but **the best set must be resubmitted
last, before the deadline** - currently submission 23 at
`submissions/23_esmfold2_default/submission_esmfold2_default.zip`.

The field has also moved: the leader has gone 0.4812 -> 0.4906 since 2026-10-22, and
three entrants now sit above 0.47. Our best is 0.038 off the lead. Note that the three
strongest BiSyRMSD scores in the field (3.29, 3.34, 3.33) are far below anything we have
produced - our best is submission 19 at 3.6754 - which suggests the leaders are doing
something we are not on pose geometry specifically.


## Submission 29: selection DOES work - and the mode was the confound

**0.4638 / 3.7277 / 0.8001. New best LDDT-PLI of the project, and the first submission
ever to beat its best parent on all three metrics simultaneously.**

| | 23 best parent | 26 | 27 | **29 Glide inplace** | 28 RTCNN |
|---|---|---|---|---|---|
| LDDT-PLI | 0.4527 | 0.4365 | 0.4384 | **0.4638** | 0.4411 |
| BiSyRMSD | 3.8551 | 3.9024 | 3.8017 | **3.7277** | 3.7722 |
| LDDT-LP | 0.7985 | **0.8063** | 0.7908 | 0.8001 | 0.7965 |

Against its best parent: **+0.0111 LDDT-PLI, -0.1274 BiSyRMSD, +0.0016 LDDT-LP.**
Against submission 28, the same three sets selected by RTCNN instead: **+0.0227**, from
only 10 differing structures.

### This overturns the project's main negative conclusion

Every Glide selection that failed here - submissions 13 (-0.056), 18, 20, 24 (-0.018
against its best parent) - used `DOCKING_METHOD mininplace`. Every one. The single time
`inplace` has been used to select, it won, and by a clear margin.

That is a confound we introduced ourselves, on one argument: `inplace` returns emodel and
einternal as 0.00, and emodel is Glide's pose-ranking term. The argument was never
tested, and it only matters if emodel is used - the selections all ranked on gscore
alone, so nothing was being lost by using `inplace`.

Analysis 28 had already shown the two modes are not equivalent (r=0.791, 7 of 20 picks
change) and that `mininplace` gives the WORSE score for 38 of 60 structures, median
+0.177. The likely mechanism: restrained minimisation moves the pose off the co-folded
geometry and POSTDOCK rescores from there, so on co-folded input the displacement costs
more than the relaxation gains. Scoring the pose as it stands is the right thing to do
for structures we did not dock.

**So "no scoring function can select structures here" is not supported.** What is
supported is that `mininplace` gscore cannot. `inplace` gscore has one prospective test
and it won. One result is not a rule, but it is the opposite of what was recorded.

### What should be retested

Submissions 13, 18, 20 and 24 are all natural candidates to repeat under `inplace`:
13 and 18 within an ensemble, 24 across three models. Analysis 23's finding that Glide
agrees with PL-iPTM at chance on the Boltz ensemble was also measured with `mininplace`
and should be remeasured.

### Standing

0.4638 would place **4th** on the 2026-10-26 table, ahead of TCB at 0.4578, behind
CYP-Random at 0.4750. Remember the table shows the LATEST submission: 29 is now the one
displayed, which for once is also our best.


## Submission 33: the third candidate was not the problem

**0.4464 / 3.7705 / 0.7952.** Best Glide `inplace` gscore over sets 23 and 27 only,
with set 26 (max-accuracy) dropped from the pool.

| | 23 alone | 29 best-of-3 | **33 best-of-2** |
|---|---|---|---|
| LDDT-PLI | 0.4527 | **0.4638** | 0.4464 |
| BiSyRMSD | 3.8551 | **3.7277** | 3.7705 |
| LDDT-LP | 0.7985 | **0.8001** | 0.7952 |

Dropping set 26 cost **0.0174** against submission 29 and left the result **0.0063**
below plain ESMFold2. The hypothesis was that set 26 was over-picked (9 of 20 against
6.7 expected) because gscore rewards its pocket quality, which gscore cannot judge. The
test says otherwise: those 9 picks were carrying real value, and removing the candidate
cost more than any bias it introduced.

So the two findings stand together. Selection across *different models* fails
(submission 30, -0.0113) because the receptors are not comparable. Selection across
*runs of one model* works, and works better with more candidates in the pool, even when
one of those runs scores worst on its own. The pool size is the lever, not the pruning.

### Where that leaves the project

| rank | submission | LDDT-PLI |
|---|---|---|
| 1 | 29 best Glide `inplace` over 3 ESMFold2 runs | **0.4638** |
| 2 | 23 plain ESMFold2, default + MSA | 0.4527 |
| 3 | 25 | 0.4465 |
| 4 | **33 best-of-2** | 0.4464 |
| 5 | 32 z-score consensus | 0.4438 |

Submission 29 remains the only rule that has beaten plain ESMFold2.


## Submissions 34 and 35: the pool is a local optimum, and the metric is coarse

| | LDDT-PLI | BiSyRMSD | LDDT-LP | vs submission 29 |
|---|---|---|---|---|
| 29 best of 3 runs | **0.4638** | **3.7277** | 0.8001 | — |
| 34 same picks, minimised ligands | 0.4637 | 3.7483 | 0.8001 | **-0.0001** |
| 35 best of 4 runs | 0.4512 | 3.8222 | 0.7975 | **-0.0126** |

### Submission 34: a 0.497 A ligand move is worth 0.0001

Every ligand was replaced by its Glide-minimised pose, a median 0.497 A displacement
and up to 2.503 A, inside receptors that were not touched. LDDT-LP came back identical
to four decimals, which confirms the receptors really were untouched and makes this a
clean single-variable measurement of what ligand placement is worth.

**The answer is: at this scale, nothing.** LDDT-PLI moved by 0.0001 and BiSyRMSD by
0.0206. Moving every ligand half an angstrom is invisible to the primary metric.

That is the most useful calibration this project has produced, and it reframes a lot of
earlier reasoning. Analysis 31 treated a 1.47 A pose shift as a large effect and it cost
0.0411; analysis 27 measured 1.53 A between sets 23 and 27, which differ by 0.0143.
Against those, 0.5 A costing nothing says the metric responds to which contacts exist,
not to their precise lengths. Sub-angstrom refinement cannot win here. Neither can any
rescoring scheme whose effect is to nudge poses.

It also ends a different question. Post-processing a good submission was 0 for 6; this
is the seventh attempt and the first that did not lose. It did not win either - it
reproduced its parent. Changing coordinates without changing which structure is shipped
is, on this evidence, a no-op.

### Submission 35: adding a bad candidate costs, and so does removing a good one

Set 31 (msa_max_depth 2048, 0.4116 alone, the worst run of the project) won 7 of 20
picks against the 5 chance would give it. Those 7 cost 0.0126.

Put beside submission 33, which removed set 26 from the same pool and cost 0.0174, the
pair rules out the simple reading of either result:

| | change to the pool | result |
|---|---|---|
| 33 | remove set 26 (3 -> 2 runs) | -0.0174 |
| 29 | the pool as built (3 runs) | **0.4638** |
| 35 | add set 31 (3 -> 4 runs) | -0.0126 |

**Pool size is not the lever. Submission 29's three-run pool is a local optimum, and
moving in either direction costs about the same.**

The mechanism both results support: Glide gscore carries weak signal for structure
quality. Among candidates of comparable quality, picking the best-scoring one is close
to harmless and slightly helpful - submission 29 beat all three of its parents. Offer
the same selector a systematically worse run and it cannot tell: it took set 31's
structures at 1.4x the chance rate and imported their deficits. Remove a comparable
candidate and it simply has fewer good structures to find.

So the rule is narrower than "select within one model". It is **select among runs of
comparable quality within one model**. Set 31 was 0.0411 below the best of the pool and
that was far enough outside to break it; sets 23, 26 and 27 span 0.0163 and work.

### Where this leaves the project

| rank | submission | LDDT-PLI |
|---|---|---|
| 1 | 29 best Glide `inplace` over 3 ESMFold2 runs | **0.4638** |
| 2 | 34 same, with minimised ligands | 0.4637 |
| 3 | 23 plain ESMFold2, default + MSA | 0.4527 |
| 4 | 35 best of 4 runs | 0.4512 |
| 5 | 33 best of 2 runs | 0.4464 |

Submission 29 still stands, and submission 34 ties it from a different direction. Three
separate attempts to improve on it - prune the pool, extend the pool, refine the poses -
have now each failed, and failed in ways that point at the same thing: the selector has
little signal to give, and the metric is too coarse for pose-level work to matter.
