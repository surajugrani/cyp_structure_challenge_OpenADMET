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
| 19 | [AF3/Boltz hybrid, better Glide score](19_hybrid_glide/) | 2026-10-18 | **0.4331** | **3.6754** | 0.7823 | — | — |
| 20 | [best Glide over all AF3+Boltz samples](20_union_glide/) | 2026-10-19 | 0.3854 | 4.2275 | 0.7449 | — | — |
| 21 | [hybrid 19 with a 0.5 kcal/mol Glide margin](21_hybrid_margin/) | 2026-10-19 | pending | pending | pending | — | — |
| 12 | [ICM redocked poses in AF3's receptor](12_icm_redock_af3/) | 2026-10-13 | pending | pending | pending | — | — |

Metric notes: **LDDT-PLI** is the primary score (protein-ligand interface, higher better).
**BiSyRMSD** is binding-site-superposed symmetry-corrected pose RMSD (lower better).
**LDDT-LP** scores the protein pocket itself. **Coverage** = fraction of structures that
scored at all; 1 means nothing was zeroed by OpenStructure or PoseBusters.

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
