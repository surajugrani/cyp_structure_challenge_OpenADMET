# CYP3A4 Structure Prediction Challenge (OpenADMET)

Submissions to the structure prediction track of the OpenADMET CYP3A4 challenge: predict
protein–ligand complex structures of human cytochrome P450 3A4 (UniProt
[P08684](https://www.uniprot.org/uniprotkb/P08684)) for 20 blinded ligands.

Scoring uses three metrics. **LDDT-PLI** (higher better) is the primary one and measures
the protein–ligand interface. **BiSyRMSD** (lower better) is binding-site-superposed,
symmetry-corrected pose RMSD. **LDDT-LP** (higher better) scores the protein pocket
itself and is largely independent of where the ligand sits. While the challenge is open
the organizers score a subset of roughly 10 of the 20 structures, so small differences
are noisy and three of our ligands appear never to be scored at all (see finding 5).

**Where things stand:** best submission is 23 at **0.4527 LDDT-PLI**, 5th on the public
leaderboard as of 2026-10-22, 0.024 behind the leader.

## Results

Not every attempt is listed. These are the ones that either scored well or settled a
question — several of the most useful were the ones that failed. Full history, including
the attempts omitted here, is in
[`submissions/LEADERBOARD.md`](submissions/LEADERBOARD.md).

| # | Submission | LDDT-PLI | BiSyRMSD | LDDT-LP | What it showed |
|---|---|---|---|---|---|
| 03 | [Boltz-2 defaults](submissions/03_boltz2_prepwizard/) | 0.3921 | 4.1179 | 0.7962 | First baseline. PrepWizard minimisation is what makes structures pass PoseBusters |
| 06 | [AF3, PL-iPTM selection](submissions/06_af3_pliptm/) | 0.4291 | 4.0647 | 0.6985 | Best AF3 variant. Nothing beat AF3's own confidence on AF3's ensemble |
| 07 | [AF3 + cryo-EM templates](submissions/07_af3_templates/) | 0.3853 | 4.2613 | 0.6762 | **Templates hurt.** Holo templates pulled predictions toward the wrong ligand |
| 13 | [AF3 samples picked by Glide](submissions/13_af3_glide_select/) | 0.3730 | 4.5822 | 0.6968 | **Glide cannot rank within one model's ensemble** — 0.056 below confidence-based picking |
| 14 | [AF3 samples picked by ICM RTCNN](submissions/14_af3_icm_rtcnn/) | 0.3949 | 4.4057 | 0.7004 | Same failure with a second scoring function, so it is the regime, not the program |
| 16 | [Boltz-2 tuned, PL-iPTM](submissions/16_boltz2_pliptm/) | 0.4229 | 3.9283 | **0.8060** | First selection rule that ever gained anything (+0.0166). Best LDDT-LP |
| 18 | [Boltz-2, Glide selection](submissions/18_boltz2_glide/) | 0.4174 | 3.7523 | 0.8047 | Confirms 13 on a different ensemble |
| 19 | [AF3/Boltz hybrid by Glide](submissions/19_hybrid_glide/) | 0.4331 | **3.6754** | 0.7823 | **Glide can choose between whole models**, beating both parents. Best BiSyRMSD |
| 20 | [Best Glide over all samples](submissions/20_union_glide/) | 0.3854 | 4.2275 | 0.7449 | Collapsing 19's choice to a free-for-all destroys the gain — it is a scale effect |
| 21 | [19 with a Glide margin](submissions/21_hybrid_margin/) | 0.4174 | 3.7523 | 0.8047 | Scored byte-identically to 18 despite differing structures, exposing the unscored subset |
| 22 | [19 with 5 Glide-ranked ICM poses](submissions/22_icm_fe_swap/) | 0.4160 | 3.7855 | 0.7823 | Fe-distance heuristic's third defeat; imported docked poses already known to be worse |
| 23 | [ESMFold2 + ColabFold MSA](submissions/23_esmfold2_default/) | **0.4527** | 3.8551 | 0.7985 | **Best result.** A third co-folding model beat both others with no selection step at all |
| 24 | [Glide over three models](submissions/24_hybrid_three_way/) | 0.4348 | 3.7261 | 0.7798 | **Kills the cross-model rule.** +0.0017 over 19, but 0.0179 below just using ESMFold2 |
| 25 | [23 with 3 structures from 19](submissions/25_esm_outlier_fix/) | 0.4465 | 3.9500 | 0.7969 | Even correcting physically impossible Fe distances lost 0.0062 |
| 26 | [ESMFold2 max-accuracy endpoint](submissions/26_esmfold2_maxacc/) | 0.4365 | 3.9024 | **0.8063** | Best pocket of the project, worse interface — the two metrics come apart |
| 27 | [ESMFold2, no MSA](submissions/27_esmfold2_nomsa/) | 0.4384 | 3.8017 | 0.7908 | One-variable MSA control: the ColabFold MSA is worth +0.0143 LDDT-PLI |
| 28 | [best RTCNN over 3 ESMFold2 runs](submissions/28_rtcnn_best_of_three/) | 0.4411 | 3.7722 | 0.7965 | Best BiSyRMSD of the ESM family; beats 2 of 3 parents on LDDT-PLI |
| 29 | [best Glide `inplace` over 3 ESMFold2 runs](submissions/29_glide_inplace_best_of_three/) | **0.4638** | 3.7277 | 0.8001 | **Best result.** Selection works — with `inplace`, not `mininplace` |
| 30 | [Glide `inplace` over three models](submissions/30_inplace_three_way/) | 0.4414 | 3.7178 | 0.7808 | **Selection works within a model, not across one.** `inplace` is worth +0.0066, still 0.0113 under ESMFold2 alone |
| 31 | [ESMFold2, MSA depth 1024→2048](submissions/31_esmfold2_depth2048/) | 0.4116 | 4.1680 | 0.7894 | **Biggest single-knob loss, −0.0411.** Pocket intact, ligands worse: the dial moves poses the wrong way |
| 32 | [z-score consensus of Glide + RTCNN](submissions/32_zscore_consensus/) | 0.4438 | 3.7662 | 0.8014 | Landed between its two inputs — consensus averaged them, the decorrelation bought nothing |
| 33 | [best Glide `inplace` over 2 ESMFold2 runs](submissions/33_glide_best_of_two/) | 0.4464 | 3.7705 | 0.7952 | Dropping set 26 cost 0.0174 — the pool mattered more than the suspected bias |
| 34 | [submission 29 with Glide-minimised ligand poses](submissions/34_glide_minimised/) | 0.4637 | 3.7483 | 0.8001 | **A 0.497 Å median ligand move cost 0.0001.** LDDT-PLI is blind at this scale |
| 35 | [best Glide `inplace` over 4 ESMFold2 runs](submissions/35_glide_best_of_four/) | 0.4512 | 3.8222 | 0.7975 | Adding a weak candidate cost 0.0126 — with 33, the three-run pool is a local optimum |
| 36 | [best Glide `emodel` over 3 ESMFold2 runs](submissions/36_emodel_best_of_three/) | pending | pending | pending | Same pool as 29, ranked by Glide's pose term instead of its affinity term; 8/20 change |

## What the project has learned

1. **Co-folding beats docking for this target.** Of the 15 blinded ligands that can
   coordinate the haem iron, co-folded structures place 11–12 in contact with it; ICM
   redocking into the same receptors manages 4 and Glide 0. Glide systematically drifts
   ligands away from the haem. Analyses
   [17](analyses/17_icm_redock_predicted/), [18](analyses/18_glide_redock_predicted/).

2. **Glide selection works — but only with `DOCKING_METHOD inplace`, and we spent most
   of the project using `mininplace`.** Submission 29 picked per ligand among three
   ESMFold2 runs by `inplace` gscore and beat its best parent on all three metrics
   (0.4638 vs 0.4527), the first submission ever to do so. Every Glide selection that
   *failed* — 13, 18, 20, 24 — used `mininplace`. That mode gives the worse score for 38
   of 60 structures and changes 7 of 20 picks
   ([analysis 28](analyses/28_inplace_vs_mininplace/)), apparently because its restrained
   minimisation moves the pose off the co-folded geometry before `POSTDOCK` rescores it.
   The only argument for `mininplace` was that `inplace` zeroes emodel — irrelevant when
   selecting on gscore alone. Submissions 13, 18, 20 and 24 are all worth repeating.

3. **The model matters more than almost anything done afterwards.** The single largest gain in the project came from trying a third
   co-folding method, not from any selection, docking or refinement step: ESMFold2
   scored 0.4527 against 0.4331 for the best hybrid and 0.4291 for the best AF3 variant.
   Against that, seven attempts to improve a good submission by rearranging or replacing
   its structures — submissions 11, 20, 21, 22, 24, 25 and 26 — all lost ground,
   including one (25) that only corrected Fe distances which are not physically possible.
   Submission 29 is the single exception, and finding 2 explains why it differs.

4. **"Closer to the haem" is not a proxy for "more correct."** Both AF3 and Boltz
   already place ligands in the right pocket ([analysis 16](analyses/16_fe_distance_swap/));
   the error is orientation within it. Submissions that swapped structures for ones
   nearer the iron lost ground. Among crystal structures, ligands with no nitrogen donor
   never coordinate the iron ([analysis 19](analyses/19_crystal_fe_basterate/)), so a
   short Fe contact on such a ligand is a red flag rather than a success.

5. **Three ligands are never scored.** Submissions 18 and 21 returned byte-identical
   scores while differing on exactly OCNT-2312966, OCNT-2315379 and OCNT-2535519, which
   places those three outside the ~10-structure scored subset. Any measured difference
   that rests only on them is invisible.

6. **The MSA is worth +0.0143 LDDT-PLI, measured cleanly.** Submission 27 repeats
   submission 23 with `msa=None` and nothing else changed — same endpoint, weights,
   sequence and chain layout — scoring 0.4384 against 0.4527. It is the project's only
   one-variable comparison. The MSA also gains +0.0077 LDDT-LP but *costs* 0.053 on
   BiSyRMSD, so it is worth keeping without being decisive.

7. **Fe distance hides orientation changes; use ligand RMSD instead.** Submissions 23
   and 27 agree within 0.3 Å on Fe–ligand distance for 18 of 20 ligands, yet their
   pocket-superposed ligand RMSD has a median of 1.53 Å and a maximum of 3.23 Å
   ([analysis 27](analyses/27_sub23_vs_nomsa/)). Fe distance is one scalar per complex
   and is blind to rotation about the coordinating atom — the very error mode finding 4
   describes. Keep it for spotting impossible values, not for asking whether two
   predictions differ.

8. **Templates hurt.** Both holo cryo-EM templates (submission 07) and a single holo
   template (08, 09) scored below the template-free prediction, apparently by biasing
   the pocket toward the template's own ligand.

## Layout

```
input/        challenge inputs: sequence, blinded ligands, RCSB report
analyses/     one self-contained folder per question, each with its own README
submissions/  one folder per real submission attempt, each ending in a .zip
tools/        shared helpers reused across analyses and submissions
figures/      generated SVGs, regenerated from committed data by tools/
```

Every analysis folder holds the scripts, the SLURM job that ran them, the results CSV
and a README stating the question, the decisions taken and how to rerun. A chronological
record is in [`PROJECT_LOG.md`](PROJECT_LOG.md) and all returned scores are in
[`submissions/LEADERBOARD.md`](submissions/LEADERBOARD.md).

Figures regenerate from committed structures and CSVs:

```bash
python3 tools/plot_fe_distances.py          # every configured figure
python3 tools/plot_fe_distances.py --list   # what is configured
```

## Analyses

| # | Question | Outcome |
|---|---|---|
| 01 | [Glide score of crystal poses vs IC50](analyses/01_glide_score_vs_ic50/) | No usable correlation (rho = −0.27) |
| 02 | [Can Glide reproduce the crystal pose?](analyses/02_redocking_pose_recovery/) | 18 % top-1 success at 2 Å |
| 03 | [Haem representation and Glide's metal term](analyses/03_heme_representation/) | Metal term activation follows no chemical logic |
| 04 | [Does the wrong chemistry dock better?](analyses/04_wrong_chemistry_replication/) | No; the metal reward is unreachable |
| 05 | [Redocking from SMILES, all structures](analyses/05_full_pdb_redocking/) | 97/100 docked |
| 06 | [Redocking with ICM-Pro](analyses/06_icm_redocking/) | 100/100 docked |
| 07 | [Blinded ligands vs PDB chemistry](analyses/07_test_vs_pdb_chemistry/) | Sets barely overlap; the redocking benchmark is pessimistic here |
| 08 | [ICM scoring of crystal poses vs IC50](analyses/08_icm_score_only/) | No correlation; an RTCNN subset hit is spurious |
| 09 | [Which score picks the native pose?](analyses/09_pose_selection/) | Glide 2/9, ICM 3/9 top-1 on 9 rigid ligands |
| 10 | [Which AF3 metric picks a valid pose?](analyses/10_af3_selection_vs_posebusters/) | No confidence metric does; PoseBusters itself reaches the 19/20 ceiling |
| 11 | [Cryo-EM preprint insights](analyses/11_preprint_insights/) | Pocket plasticity; motivated the template attempts |
| 12 | [Does a metal constraint fix coordination?](analyses/12_metal_constraint/) | No |
| 13 | [How far apart are AF3's own samples?](analyses/13_pose_divergence/) | 1–2 Å typical — the regime where scoring fails |
| 14 | [Do Glide scores of AF3 poses predict potency?](analyses/14_glide_vs_ic50/) | No, including vdW and electrostatic terms separately |
| 15 | [How big a bet is each selection rule?](analyses/15_candidate_rules/) | Glide is the wildest and worst-motivated bet, 8 structures moving >2 Å |
| 16 | [Is AF3 misplacing the pocket?](analyses/16_fe_distance_swap/) | No — the pocket is right, the orientation is not |
| 17 | [ICM redocking into predicted receptors](analyses/17_icm_redock_predicted/) | 4/15 coordinating, worse than co-folding |
| 18 | [Glide redocking into predicted receptors](analyses/18_glide_redock_predicted/) | 0/15 coordinating; ligands drift off the haem |
| 19 | [Do crystal ligands actually coordinate Fe?](analyses/19_crystal_fe_basterate/) | 0 of 7 ligands without a nitrogen donor do |
| 20 | [ICM scoring of the AF3 samples](analyses/20_icm_score_af3/) | Ends the within-ensemble selector question |
| 21 | [Why does Boltz beat AF3 on LDDT-LP?](analyses/21_pocket_lddt_vs_crystals/) | AF3 reuses one pocket (6DAA for 15/20); Boltz adapts it per ligand |
| 22 | [Pose movement between submissions](analyses/22_sub16_vs_17/) | Pocket-superposed RMSD between any two submissions |
| 23 | [Would Glide pick better Boltz samples?](analyses/23_glide_on_boltz/) | No — agreement with PL-iPTM is at chance |
| 24 | [Preprint pocket residues](analyses/24_preprint_residues/) | Reference residue set for pocket comparisons |
| 25 | [Glide rescoring of ICM poses](analyses/25_glide_on_icm_poses/) | Recovers 4→7 of 15 coordinating, still below co-folding |
| 26 | [Glide scoring of the ESMFold2 structures](analyses/26_glide_on_esm/) | Feeds the three-way hybrid, submission 24 |
| 27 | [ESMFold2 with and without the MSA](analyses/27_sub23_vs_nomsa/) | Median 1.53 Å ligand RMSD apart — the MSA changes poses, not just Fe distances |
| 28 | [`inplace` vs `mininplace` on 60 ESMFold2 structures](analyses/28_inplace_vs_mininplace/) | r=0.791, `mininplace` worse on 38/60, 7/20 picks change — the mode was the confound |
| 29 | [ICM RTCNN on the ESMFold2 structures](analyses/29_icm_on_esm/) | 60/60 scored; Glide and RTCNN agree on only 10/20 winners |
| 30 | [`inplace` scores for the AF3 and Boltz-2 sets](analyses/30_inplace_three_model/) | 40 structures, feeds the rebuilt cross-model hybrid, submission 30 |
| 31 | [What the ESMFold2 config knobs do](analyses/31_msa_depth/) | Compute knobs already maxed; ESMFold2 reproduces to 0.31 Å; MSA depth moves poses 1.47 Å |
| 32 | [Glide `inplace` on the depth-2048 set](analyses/32_inplace_depth2048/) | Fourth candidate for the pool; it wins 7/20 despite being the worst run |
| 33 | [ICM redock into submission 29's receptors](analyses/33_icm_redock_sub29/) | 200 poses, 10 per ligand, into the best receptors we have |
| 34 | [Glide `inplace` on the ICM poses](analyses/34_inplace_on_icm_poses/) | 177/200 rejected — overlapping hydrogens, not pose quality |
| 35 | [Does ICM change the receptor?](analyses/35_icm_receptor_check/) | Yes, and it is not the cause: the clash is H···H at 0.73 Å |

## Methods

| Tool | Role |
|---|---|
| [AlphaFold3](https://github.com/google-deepmind/alphafold3) | Co-folding, 5 samples per ligand |
| [Boltz-2](https://github.com/jwohlwend/boltz) | Co-folding, tuned recycling and diffusion sampling |
| [ESMFold2](https://github.com/evolutionaryscale/esm) | Co-folding via the remote Biohub API |
| [Schrödinger](https://www.schrodinger.com/) | PrepWizard preparation, Glide docking and scoring |
| [Molsoft ICM-Pro](https://www.molsoft.com/icm_pro.html) | Docking and RTCNN scoring |
| [PoseBusters](https://github.com/maabuu/posebusters) | Physical validity of every shipped structure |
| [OpenStructure](https://openstructure.org/) | Local lDDT, matching the organizers' scoring |

Every submission is validated before packaging: chain layout and ligand connectivity
against the input SMILES, then PoseBusters in `dock` configuration. Nothing ships below
20/20 on both.

## Related work

This builds on the methods developed for the
[OpenADMET PXR structure challenge](https://github.com/surajugrani/pxr-structure-challenge).
