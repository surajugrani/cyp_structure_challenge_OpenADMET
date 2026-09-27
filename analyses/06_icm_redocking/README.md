# 06 — Pose recovery with ICM-Pro

**Question:** Does ICM-Pro recover CYP3A4 crystal poses better than Glide, given the same
ligands and the same binding site?

**Status:** Running, following behind analysis 05.

## Design

| | Glide (analysis 05) | ICM (this analysis) |
|---|---|---|
| Starting structure | `<PDB>_clean.mae` | **same** `<PDB>_clean.mae` |
| Receptor prep | PrepWizard + `normalize_heme` | **ICM `convertObject`** |
| Ligand | LigPrep states from SMILES | **the same LigPrep states** |
| Site | grid centred on crystal ligand | residues within 5 Å of crystal ligand |
| Poses | 10 | 10 |
| RMSD | spyrmsd | same spyrmsd code path |

**Receptor preparation is deliberately ICM's own.** `_dockBatch` calls `convertObject` with
hydrogen optimisation, His optimisation and **hetero charge assignment** — the last is what
charges the heme, and is precisely what the PXR-era invocation
(`convertObject a_rec. no no no no no no no ""`) switched off. ICM therefore starts from
the *pre-PrepWizard* cleaned structure; feeding it the PrepWizard output would be double
preparation and would inherit the heme problems documented in analyses 03–04.

**Ligand preparation is held constant** (the same LigPrep states Glide received), so a
difference in pose recovery is attributable to the docking pipeline rather than to ligand
setup.

**Caveat:** prep and engine therefore change together. This compares *pipelines*, which is
the right question for deciding what to run in the challenge, but it cannot attribute a
difference to the docking algorithm alone.

## Implementation notes

Three things that are easy to get wrong here, all found the hard way:

1. **V3000 vs V2000.** Schrödinger writes V3000 MOL blocks; ICM rejects them with
   *"no atoms in the object"*. obabel converts them. (The PXR repo's
   `v-3000-to-2000_obabel.sh` existed for exactly this reason.)

2. **Scores must come from `output=hits.sdf`, not from stdout.** `_dockScan` prints
   `SCORES>` lines in a *different order* from both the `.ob` conformation stack and the
   hitlist. Parsing them positionally silently mis-assigns every score — on 4D6Z it gave
   the 2.81 Å pose a score of −13.94 when the true value is −12.85. `hits.sdf` binds each
   score to its own molecule as SD tags. Verified: `hits.sdf` is ordered worst→best by
   `Score` and is the exact reverse of the `.ob` stack.

3. **ICM is stochastic; Glide is not.** `_dockScan` is Monte Carlo, so runs are only
   reproducible with `seed=`. (Glide was verified deterministic: 3/3 identical repeats.)
   `_dockBatch` is split from `_dockScan` purely so `seed=` and `confs=10` can be passed —
   `_dockBatch` alone hardcodes `confs=3`.

## Running

```bash
sbatch scripts/run_icm_follow.sh    # follows analysis 05, exits when all are done
sbatch scripts/run_icm.sh --only 4D6Z,6MA7   # single structures
```

ICM uses a separate licence from Schrödinger, so it runs alongside the Glide job. Because
ICM reuses analysis 05's LigPrep output, it can only process structures Glide has already
prepared; `run_icm_follow.sh` loops, docking whatever is ready and sleeping when it catches
up.

## Preliminary result (2 structures — not conclusive)

| | Glide top-1 | Glide best-of-10 | ICM top-1 | ICM best-of-10 |
|---|---|---|---|---|
| 4D6Z (PK9) | 5.84 | 5.00 | 8.38 | **2.81** |
| 6MA7 (TPF) | 4.68 | 3.40 | 4.05 | **2.21** |

On both, ICM's *sampling* found a better pose than Glide generated at all, while its
*scoring* failed to select it — for 4D6Z the 2.81 Å pose ranked 6th of 10. That is the same
sampling-versus-scoring split analysis 02 found in Glide, but starker: here a near-native
pose was available and simply not chosen.

Two structures prove nothing; this is recorded as something to check against the full set.
