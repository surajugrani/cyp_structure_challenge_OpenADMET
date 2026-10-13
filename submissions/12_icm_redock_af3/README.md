# Submission 12 — ICM redocked poses in AF3's receptor

A deliberate floor-finding experiment: how much worse is docking into a predicted receptor
than the co-folded pose it replaces? Submitted knowing it should score below 06.

Protein and haem from submission 06 (our best, 0.4291), ligand replaced by ICM's top-1
redocked pose from [analysis 17](../../analyses/17_icm_redock_predicted/).

**20/20 structural, 20/20 PoseBusters** — `--allow-pb-fail` was available but not needed.
The poses are physically clean; they simply disagree with the co-folded ones.

## Why AF3's receptor and not Boltz's

Analysis 17 docked into both. Despite Boltz's pocket scoring 0.098 higher on LDDT-LP,
AF3's receptor was the better docking target on every measure — top-1 RMSD to co-folded
5.04 Å vs 6.38, best-of-10 12/20 vs 6/20 under 2 Å, Fe coordination retained 4/20 vs 0/20.
Backbone accuracy and dockability are not the same property.

## What is being shipped

| | value |
|---|---|
| median top-1 RMSD from 06's pose | 5.04 Å |
| poses within 2 Å of 06's | 2/20 |
| Fe coordination retained | 4/20 (co-folded: 12/20) |

The large disagreement is the point of the experiment, not a defect in the build.

## Expected outcome

Clearly below 06's 0.4291. Two reasons: the ligand moves a median 5 Å from a pose that
already scored our best, and the Fe coordination — present in 12 of 20 co-folded
structures and the defining feature of a type II CYP binder — survives in only 4.

If it lands close to 06 anyway, that would be informative in its own way: it would mean
LDDT-PLI is far less sensitive to ligand placement than the BiSyRMSD figures suggest.

## Build notes

* ICM writes `hits.sdf` **worst-first**; pose 10 of 10 is the best-scoring one.
* An earlier version of this script omitted the altLoc space in the `HETATM` line, pushing
  the residue name into columns 17–19 so `l[17:20]` read `"IG "`. All 20 structures failed
  with "no LIG residue". PDB columns: 13–16 atom name, 17 altLoc, 18–20 resName, 22 chain.

| file | contents |
|---|---|
| `results/selection.csv` | per ligand: ICM pose rank, Score, RTCNN, RMSD from co-folded, Fe distance, validation |
