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
