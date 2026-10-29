# Submission 34 - submission 29 with Glide-minimised ligand poses

**The first submission that changes coordinates rather than rearranging which structure
goes where.**

Submission 29 (0.4638, the best result of the project) picked, for each of the 20
blinded ligands, whichever of the three ESMFold2 runs scored best under Glide
`DOCKING_METHOD inplace`. That pick is kept exactly here: same receptor, same run, same
ligand, same chain layout. The only thing that differs is the ligand's coordinates,
replaced by the pose Glide produced under `mininplace` - a restrained minimisation of
the co-folded pose inside that receptor's own field.

No new docking, no new scoring. Analysis 28 ran both modes over all 60 ESMFold2
structures to compare their scores, and its `mininplace` output poses were still on
disk. This submission is a read of files we already had.

## The poses genuinely moved

| | ligand displacement from submission 29 |
|---|---|
| median | **0.497 A** |
| minimum | 0.126 A (OCNT-2395645) |
| maximum | 2.503 A (OCNT-2312966) |

This was the check worth doing before building: at 0.2 A the submission would have been
a no-op and not worth a slot. At 0.5 A median it is a real test. For scale, analysis 31
measured ESMFold2's run-to-run reproducibility at 0.31 A, so the minimiser moves the
ligand further than re-running the model does.

All 20 matched on atom name between the Glide pose and the submission PDB, so no
fallback to file order was needed.

## The case for, and the case against

**For.** LDDT-PLI scores protein-ligand interface contact quality. Minimisation inside
the receptor field is the one operation available that relieves clashes and tightens
contacts, and the co-folded poses were never relaxed against a physics force field at
all.

**Against, recorded before the result is known.** Analysis 28 found `mininplace` gives
the *worse* gscore for 38 of 60 of these same structures, median +0.177 - the minimiser
moves the pose somewhere Glide likes less, which is the finding that made `inplace` the
right mode for scoring in the first place. And post-processing a good submission is 0
for 6 in this project (submissions 11, 20, 21, 22, 24, 25).

The two are not in direct conflict: gscore and LDDT-PLI are different objectives, and
submission 30 already showed gscore cannot judge things outside its remit. But the prior
is clearly negative.

## What the result means

- **Above 0.4638** - relaxation helps, and refinement is a direction this project has
  never explored. The obvious follow-ups are a gentler restraint and minimising the
  other two runs' poses too.
- **Below 0.4638** - the co-folded geometry is better than a force field can make it,
  and the running count of failed post-processing attempts becomes 0 for 7. In that case
  the lesson is the same one submissions 30 and 33 taught: leave the model's output
  alone, and spend effort on the model.

## Validation

Full re-validation, not a checksum copy - the coordinates changed. Chain layout, HEM
integrity, ligand connectivity against the input SMILES, then PoseBusters in `dock`
configuration: **20/20 structural, 20/20 PoseBusters.**

## Files

| path | what |
|---|---|
| `scripts/01_build.py` | grafts the minimised ligand coordinates onto submission 29's structures |
| `scripts/02_validate.py` | re-validates and packages |
| `results/displacement.csv` | per-ligand displacement and how the atoms were matched |
| `results/validation.csv` | per-ligand structural and PoseBusters result |
| `submission_glide_minimised.zip` | the 20 structures |
