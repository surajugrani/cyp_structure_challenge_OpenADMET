# 01 — Boltz-2, default settings

**Status:** complete. 20/20 predicted (job 47806615, 1:50:42) and passing all validation checks. PoseBusters not yet run.

## Inputs

One YAML per blinded ligand: full-length CYP3A4 (P08684, 503 aa) chain A, haem as CCD
`HEM` chain B, ligand SMILES chain C. MSAs from Boltz's own server. No templates, no
pocket restraints, no potentials — everything else default.

The organizers' tutorial YAML omits the haem; that is not followed here. Their own sample
predictions (`inputs/boltz-2_cyp3a4_predictions/*.cif`) do contain HEM with 43 atoms, and
the validator explicitly allows a third chain for it.

## Running it

```bash
python3 scripts/01_make_yamls.py
sbatch scripts/02_run_boltz.sh
$SCHRODINGER/run python3 scripts/03_package.py
```

## SLURM notes — four failures worth not repeating

Reusing the PXR `run_boltz.sh` verbatim would have avoided the first three.

1. `micromamba activate` fails under sbatch — the shell function lives in `~/.bashrc`,
   which batch scripts do not source. Call `~/micromamba/envs/boltz/bin/boltz` directly.
2. `module load cuda/11.8` is required; torch here is `2.7.1+cu118`.
3. `--gres=gpu:rtxpro6000:1` never schedules — both such nodes are drained/reserved. Use
   `gtx1080ti`. Also `--cpus-per-task=4`, not 8: the 1080Ti nodes have only 2–6 idle CPUs.
4. **`srun` is mandatory.** Without it the script body runs outside the job-step cgroup,
   `nvidia-smi` reports both GPUs on the node while one is allocated, and torch fails with
   "CUDA driver initialization failed" — while the identical torch call succeeds under
   `srun` on the same node.

Runtime is about 23 min per complex on a 1080 Ti, so roughly 7 h for 20.

## Packaging — Boltz names the haem `LIG` too

Boltz labels *every* non-polymer chain `LIG`, so raw output has two LIG residues (haem
chain B, 43 atoms; ligand chain C). The validator demands exactly one and would reject it.
`03_package.py` renames chain B to `HEM`, then re-runs every check the organizers do
(.pdb not .cif, 20 files, one LIG, ≤3 chains, unique residue numbers / atom names / atom
serials, ligand connectivity vs the given SMILES) and refuses to write the zip unless all
20 pass.

## Result

20/20 structures produced and all pass the replicated validation checks.
`submission_boltz2_default.zip`, 1.5 MB.

Boltz confidence is uniformly high — ligand_iptm median 0.971 (range 0.931–0.987),
complex pLDDT 0.920–0.948, no outliers. Treat that cautiously: analysis 07 found these
ligands share no Murcko scaffold with the CYP3A4 PDB set (max Tanimoto 0.34), so
self-reported confidence is not independent evidence of accuracy.

**Not yet done: PoseBusters.** It is a hard gate — failing it scores LDDT-PLI 0 — and the
organizers run it separately (`evaluation/posebusters_utils.py`) from the submission
validator replicated in `03_package.py`. Run it before uploading.
