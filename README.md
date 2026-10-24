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

## What the project has learned

1. **Co-folding beats docking for this target.** Of the 15 blinded ligands that can
   coordinate the haem iron, co-folded structures place 11–12 in contact with it; ICM
   redocking into the same receptors manages 4 and Glide 0. Glide systematically drifts
   ligands away from the haem. Analyses
   [17](analyses/17_icm_redock_predicted/), [18](analyses/18_glide_redock_predicted/).

2. **No scoring function has been shown to select structures, within a model or across
   models.** A physics score cannot rank five samples of the same ligand (submissions
   13, 14, 18, 20). Submission 19 appeared to be the exception — picking per ligand
   between AF3 and Boltz-2 by Glide score beat both parents — but submission 24 tested
   that rule with a third model and it failed: Glide preferred ESMFold2 for only 5 of
   20 ligands even though ESMFold2 alone scores 0.0196 above the hybrid, and the
   three-way result landed 0.0179 *below* simply shipping ESMFold2 unaltered. Submission
   19's margin is best read as noise on a ~10-structure scored subset.

3. **The model matters more than anything done afterwards, and post-processing has never
   once helped.** The single largest gain in the project came from trying a third
   co-folding method, not from any selection, docking or refinement step: ESMFold2
   scored 0.4527 against 0.4331 for the best hybrid and 0.4291 for the best AF3 variant.
   Against that, six separate attempts to improve a good submission by rearranging or
   replacing its structures — submissions 11, 20, 21, 22, 24 and 25 — have all lost
   ground, including one (25) that only corrected Fe distances which are not physically
   possible.

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

6. **Templates hurt.** Both holo cryo-EM templates (submission 07) and a single holo
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
