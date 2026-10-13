# Analysis 17 — ICM redocking into the predicted receptors

All 20 blinded ligands redocked into **both** predicted receptors, to answer: is docking
into a co-folded structure better than the co-folded pose itself?

**No. And the receptor that docks best is not the one with the better pocket score.**

## Setup

| receptor | from | LDDT-LP |
|---|---|---|
| `boltz_03` | submission 03 (Boltz-2 + PrepWizard) | 0.7962 |
| `af3_06` | submission 06 (AF3 + PrepWizard) | 0.6985 |

Prospective by construction: the ligand is built in 3D from SMILES with `obabel --gen3d`,
never from the co-folded coordinates, so no pose information leaks in. The co-folded
ligand is used only to define the binding site for the ICM project.

Settings match analysis 06 for comparability: `confs=10 effort=3. seed=42`. 40 dockings,
one complex at a time in a single serial job (ICM licence limit), 46 minutes total,
0 failures.

There is no ground truth for the blinded set, so "RMSD to co-folded" measures **agreement
with submission 06**, not accuracy.

## Result

| | Boltz (03) receptor | AF3 (06) receptor |
|---|---|---|
| top-1 RMSD to co-folded, median | 6.38 Å | **5.04 Å** |
| best-of-10 RMSD, median | 2.36 Å | **1.69 Å** |
| top-1 within 2 Å of co-folded | 0/20 | **2/20** |
| best-of-10 within 2 Å | 6/20 | **12/20** |
| top-1 coordinating Fe (≤2.5 Å) | 0/20 | **4/20** |

## Three findings

**1. Backbone accuracy and dockability come apart.** Boltz's pocket scores 0.098 higher on
LDDT-LP, yet AF3's receptor is the better docking target on every single measure. A
prediction made before this ran — that Boltz's receptor would be the better choice
*because* of its LDDT-LP — was wrong. LDDT-LP scores the pocket backbone against the
reference; it says nothing about whether side-chain rotamers leave a dockable cavity.

**2. ICM destroys the Fe coordination.** 12 of 20 co-folded poses coordinate the iron at
2.06–2.32 Å. ICM keeps it for 4/20 in the AF3 receptor and **0/20** in the Boltz one.
Analysis 12 already showed Glide metal constraints cannot restore this either (0/10 under
2 Å, and the constraint destroyed one structure outright). Co-folding gets the
coordination essentially for free; docking does not.

**3. Sampling is fine, ranking is not — again.** Best-of-10 median is 1.69 Å and 12/20
pools contain a pose within 2 Å of the co-folded one, but top-1 median is 5.04 Å and only
2/20 land under 2 Å. That is exactly the split analysis 09 measured on crystal receptors
(ICM top-1 27 % against a much higher oracle), and the same wall analyses 09, 10 and 13
hit: the right pose is usually in the pool, and no scoring function finds it.

## Note on the hitlist order

**ICM writes `hits.sdf` worst-first** — pose 10 of 10 carries the best Score, monotonically
improving down the file. Anything that treats the first record as the top pose will be
reading the *worst* one. Scores must be read from the SD tags, never from `_dockScan`
stdout.

Also fixed here: analysis 06 looked for an SD tag named `RTCNNscore`, but ICM writes
`RTCNN`, so that column came back empty there. Corrected in this script.

| file | contents |
|---|---|
| `results/icm_poses.csv` | 400 rows — every pose, both receptors, with RMSD to co-folded, Fe distance, and all ICM score terms |
