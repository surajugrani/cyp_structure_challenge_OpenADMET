# Analysis 23 — would Glide pick better Boltz samples than PL-iPTM?

**No. It picks at chance level.**

Submission 16 showed PL-iPTM selection on the Boltz 5-sample ensemble is worth
+0.0166 LDDT-PLI — the first selection rule in this project that helped. The obvious
follow-up is whether a physics score does better. This answers it without spending a
submission slot.

## Method

All 100 Boltz samples (20 ligands × 5 diffusion samples from submission 15) through
submission 04's proven recipe: PrepWizard `-rmsd 0.3` → Glide grid on that sample's own
receptor conformation → `mininplace` scoring. 5 concurrent, 1 h 44 m, **0 failures**.

`mininplace` rather than `inplace` because `inplace` returns emodel and einternal as 0.00
(submission 02), and emodel is Glide's own pose-ranking term.

## Result

| | agrees with PL-iPTM pick |
|---|---|
| Glide `gscore` | **4/20** |
| Glide `emodel` | **5/20** |
| chance (1 of 5 samples) | ~4/20 |

Glide is not flat — median within-ligand `gscore` spread is **1.83 kcal/mol**, against
PL-iPTM's 0.079 — so it discriminates confidently. It simply discriminates on something
uncorrelated with pose quality.

## Three independent lines now agree

| evidence | finding |
|---|---|
| analysis 09, crystal receptors with known RMSD | within-ligand ρ = **+0.086** (gscore), +0.115 (lipo) |
| submission 13, Glide selection on the AF3 ensemble | **0.3730** vs PL-iPTM's 0.4291, −0.056 |
| analysis 23, Glide selection on the Boltz ensemble | agreement at **chance** |

The failure is a property of the scoring function, not of which model generated the
candidates — which is why changing the ensemble changed nothing. Submission 14 found the
same for ICM RTCNN (0.3949, −0.034), despite its much stronger within-ligand signal of
+0.334.

## CORRECTION — this recommendation was wrong

The original text here read "do not spend a submission slot on Glide or ICM selection".
Submission 18 did exactly that and scored **0.4174 LDDT-PLI / 3.7523 Å BiSyRMSD /
0.8047 LDDT-LP — the best BiSyRMSD of the project**, and only 0.0055 LDDT-PLI below the
best rule rather than the 0.056 Glide cost on AF3.

The inferential mistake: **agreement with a working rule is not a proxy for performance.**
Two selectors can disagree on which sample to take and both land on good poses, because
a 5-sample ensemble can contain several acceptable candidates. Chance-level agreement
measures independence, not quality. Only the leaderboard settles it.

The related mistake: **a score's failure on one ensemble does not transfer to another.**
Glide's -0.056 on AF3's samples did not predict its behaviour on Boltz's.

## Why PL-iPTM works where these do not

It is not merely dynamic range — Glide has more of that than PL-iPTM in absolute terms.
It is that PL-iPTM is the model's *own* estimate of how well it placed the ligand, formed
with access to the generation process, whereas a docking score re-judges the finished
pose with a physics model that was never calibrated on this problem.

| file | contents |
|---|---|
| `results/glide_scores.csv` | all 100 samples with the full GlideScore term set |
| `results/pick_comparison.csv` | per ligand: PL-iPTM pick, gscore pick, emodel pick, gscore spread |
