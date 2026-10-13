# Analysis 15 — how big a bet is each candidate selection rule?

A cheap read before spending a submission slot. For each candidate rule, how many of the
20 structures would change versus submission 06 (our best at 0.4291), and how far does the
ligand actually move? No zip is built.

The motivation is analysis 13: a single 3.5 Å pose swap cost submission 09 0.026 LDDT-PLI.
Knowing the number and size of changes up front tells you the variance you are buying.

## Result

| rule | differs from 06 | median shift | max shift | shifts >2 Å |
|---|---|---|---|---|
| best Glide gscore | 16/20 | 2.23 Å | 6.65 Å | 8 |
| consensus centroid, no-template pool | 14/20 | 1.17 Å | 3.04 Å | 4 |
| consensus centroid, pooled 04+08 | 15/20 | 1.68 Å | 5.43 Å | 5 |

The pooled centroid picks a template-run (08) pose for 9 of 20 ligands.

## Reading it

**Glide is the wildest bet and the worst-motivated.** 8 structures moving more than 2 Å,
one by 6.65 Å. Against it: analysis 10 found Glide selected PoseBusters-valid poses *below*
random (10/20 vs 11.4/20 expected) and mildly preferred the poses PB rejects, and analysis
14 found its scores anti-correlated with potency. High variance with no evidence of skill
is the combination to avoid.

**Consensus is the moderate option.** Changing 14/20 with a 1.17 Å median shift is a
smaller perturbation than Glide while still being a genuinely different rule. It is also
the only candidate built on a signal orthogonal to iPTM — agreement among independent
predictions rather than a score of one pose.

**Pooling the two runs adds variance without obvious benefit.** The pooled centroid moves
further (1.68 Å median, max 5.43) and imports template-run poses for 9 ligands, and
submissions 07/08/09 established that templates do not help on CYP3A4.

None of these was submitted. Analysis 16 then tested a physically motivated rule instead.

| file | contents |
|---|---|
| `results/candidate_rules.csv` | per ligand: 06's pick, each rule's pick, and the shift |
