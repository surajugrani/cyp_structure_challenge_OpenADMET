# Submission 23 — ESMFold2, a third co-folding method

First use of **ESMFold2** (EvolutionaryScale), with the same ColabFold MSA that AF3 used.
19 structures are pure ESMFold2; OCNT-2313038 comes from submission 19 (see below).
20/20 structural, 20/20 PoseBusters.

## How this one differs operationally

**The folding is not run here.** ESMFold2 is a remote API — compute happens on Biohub's
GPUs and this machine is only an HTTPS client, so there is no SLURM job for the
prediction and it must run where there is outbound internet (login node: yes; compute
nodes: no).

**It is fast and deterministic.** ~28 s per ligand with the MSA (21 s without), so all 20
in under 10 minutes against ~3 h for the tuned Boltz run. Two identical calls returned
**byte-identical** mmCIF, so there is no ensemble: `num_sampling_steps` is the number of
denoising steps within one trajectory, not a count of samples, and `FoldingConfig` has no
`diffusion_samples`. The selection strategy that produced our best results has nothing to
select from.

**`pair_chains_iptm` comes back `None`** despite `include_pair_chains_iptm=True`, so the
PL-iPTM analogue that drove submission 16 is unavailable. Only the saturated global
`iptm` (0.875–0.966) is returned.

## Setup notes

| | |
|---|---|
| env | `esmfold2` (Python 3.12) — `esm` >= 3.4 **requires 3.12**; on 3.11 pip silently keeps 3.2.1, which has no `esmfold2_client` |
| torch | CUDA wheels were installed by accident; **CPU-only is correct** for a client, ~2 GB smaller |
| MSA | `cyp3a4_colabfold.a3m`, 14,301 sequences, shared by all 20; MMseqs2 NUL bytes stripped first (they broke AF3) |
| config | defaults plus `include_pae` and `include_pair_chains_iptm` |
| token | `~/.cyp_biohub_token`, mode 600, git-ignored |

Chain layout needs **no fixing**: ESMFold2 emits A = protein, B = HEM (exactly 43 atoms),
C = LIG, unlike Boltz which labels every het chain `LIG`.

## Fe geometry — competitive

| | coordinating (of 15 capable) | median Fe |
|---|---|---|
| 06 AF3 | 12 | 2.25 Å |
| 19 hybrid (best) | 11 | 2.29 Å |
| 16 Boltz | 11 | 2.29 Å |
| **23 ESMFold2** | **12** | **2.16 Å** |

Matches the best coordination count with the tightest median. Two specific defects:
**OCNT-2315379 ejected to 7.11 Å** (it has an N donor — a clear miss, and the same ligand
ICM threw to 8.62 Å), and **OCNT-2395460 at 1.64 Å**, far shorter than a real Fe–N bond
(2.0–2.3 Å) on a ligand with **no** nitrogen donor. PoseBusters flagged neither — it does
not check metal geometry.

## The one swapped structure

Raw ESMFold2 scored 17/20 PoseBusters: two `minimum_distance_to_protein` clashes and one
`non-aromatic_ring_non-flatness`. PrepWizard at `-rmsd 0.3` fixed both clashes → 19/20.

The ring would not flatten. **Dropping `-rmsd` does nothing — 0.3 is PrepWizard's own
default** (`-r RMSD, -rmsd RMSD   Minimization: RMSD cutoff in angstrom (default 0.3)`),
a mistake caught before it wasted a run. Re-minimising at an explicit `-rmsd 1.0` moved
the ligand only **0.19 Å mean / 0.32 Å max** and the ring stayed flagged, so the restraint
was never the limiting factor. OCNT-2313038 was therefore taken from submission 19.

| file | contents |
|---|---|
| `raw/` | ESMFold2 mmCIF output, one per ligand |
| `results/fold_metrics.csv` | per ligand: runtime, pTM, ipTM, mean pLDDT |
| `results/validation.csv` | per ligand: atom counts, grafted atoms, source, validation |
