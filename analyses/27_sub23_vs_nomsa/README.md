# Analysis 27 — what does the ColabFold MSA actually change?

**Question.** Submission 27 is submission 23 with one variable altered: same ESMFold2
`fold_all_atom` endpoint, same `esmfold2-2026-05` weights, same sequence and chain
layout, `msa=None`. How different are the resulting ligand poses?

**Answer: substantially — and the Fe distance completely hides it.**

| measure | 23 vs 27 |
|---|---|
| median Fe–ligand distance | 2.24 vs 2.20 Å — looks identical |
| structures within 0.3 Å on **Fe distance** | 18 of 20 |
| median **ligand RMSD** (pocket-superposed, symmetry-corrected) | **1.53 Å** |
| ligand RMSD > 1 Å | 15 of 19 |
| max ligand RMSD | 3.23 Å (OCNT-2395342) |

The two sets place ligands at nearly the same *radius* from the haem iron while putting
them in materially different *orientations*. Reading only the Fe distance, as the first
pass here did, gives the flatly wrong conclusion that the MSA does almost nothing.

OCNT-2313038 returns exactly 0.000 and is the control: it is the one structure copied
from submission 23 into 27, so anything other than zero would indict the measurement.

## Why this matters beyond the MSA

Fe distance has been the project's main structural proxy, and this is the second time it
has concealed a large change. Submission 24 swapped five structures that moved 3.6–4.2 Å
by ligand RMSD while their Fe distances shifted by at most 0.07 Å. **Fe distance is one
scalar per complex and is insensitive to rotation about the coordinating atom**, which is
precisely the error mode both AF3 and Boltz were already known to have (analysis 16: the
pocket is right, the orientation is not).

Use ligand RMSD when asking whether two predictions differ. Keep Fe distance for what it
is good at: flagging values that are chemically impossible, such as a 1.79 Å contact on a
ligand with no nitrogen donor.

## Reproduce

```bash
sbatch analyses/27_sub23_vs_nomsa/run.sh
```

Uses `tools/ligand_rmsd.py`, which generalises the comparison previously copied three
times into analysis 22. Results in `results/sub23_vs_27.csv`.

## Resolved: the pose differences cost 0.0143

Submission 27 scored **0.4384 / 3.8017 / 0.7908** against submission 23's
0.4527 / 3.8551 / 0.7985. The ColabFold MSA is therefore worth **+0.0143 LDDT-PLI** and
+0.0077 LDDT-LP, while *costing* 0.053 on BiSyRMSD.

The median 1.53 Å of pose movement measured here is real and it matters, which the Fe
distances flatly denied. Had we judged on Fe distance alone we would have dropped the
MSA and lost 0.0143 on the primary metric.
