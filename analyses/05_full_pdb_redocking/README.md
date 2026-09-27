# 05 — Pose recovery across all CYP3A4 structures, docking from SMILES

**Question:** How often can Glide reproduce the crystal pose when the ligand is built from
a SMILES string, across every usable CYP3A4 crystal structure?

**Status:** Running.

## Why this differs from analysis 02

| | Analysis 02 | Analysis 05 |
|---|---|---|
| Structures | 22 (those with BindingDB affinities) | **102** |
| Ligand source | Crystal geometry | **SMILES → LigPrep (Epik, pH 7.4)** |
| RMSD | Glide `rmsd_to_input` | **spyrmsd** (graph-matched, symmetry-corrected) |

Analysis 02's `confgen` regenerated torsions but inherited crystal bond lengths, angles
and stereochemistry. Here the only crystal-derived inputs are the **receptor** and the
**grid centre** — the same situation as the 20 blinded challenge compounds, where we have
a SMILES string and nothing else. Expect worse numbers than 18 % top-1; the gap measures
how much the crystal-geometry shortcut was flattering the earlier estimate.

Pose prediction does not need affinity data, so widening to 102 structures also tightens
the confidence interval on the recovery rate considerably.

## Set selection

102 of the 122 entries in the RCSB report. Dropped:

- **6 apo structures** (heme only, or heme + cryoprotectants)
- **12 iridium complexes** — most docking programs have no Ir parameters
- **4 ligands under 150 Da** (imidazole in 9YK4 and 4D6Z, metformin in 5G5J, PG0 in 9BV6).
  Both imidazole and metformin are genuine ligands, but with so few heavy atoms a 2 Å RMSD
  cutoff is nearly free and would inflate the success rate.

92 unique ligands across the 102 structures; all have stereo-defined SMILES from the RCSB
chemical component dictionary, so no stereoisomer enumeration is needed (`ligprep -g`).

## RMSD

Building the ligand from SMILES changes the atom ordering, so Glide's
`r_i_glide_rmsd_to_input` becomes meaningless — it measures distance from the LigPrep
starting conformer, not from the crystal pose. RMSD is computed with **spyrmsd**
(graph-matched, symmetry-corrected, no superposition).

Validated against Glide's own values on crystal-geometry poses: exact agreement on 6 of 8
test poses; the two that differ (0.16 and 0.08 Å) are both on TPF, a ligand with genuine
symmetry, where spyrmsd is the stricter of the two.

spyrmsd is installed locally rather than into the read-only Schrödinger tree. Install with
`--no-deps` — a plain install pulls NumPy 2.x, which shadows Schrödinger\'s NumPy 1.24 and
breaks its compiled extensions:

```bash
module load schrodinger/2025-3
$SCHRODINGER/run python3 -m pip install --no-deps --target=$PWD/.pylibs spyrmsd
```

`scripts/run_all.sh` puts `.pylibs` on `PYTHONPATH`.

## Running

```bash
sbatch scripts/run_all.sh --limit 2   # smoke test
sbatch scripts/run_all.sh             # all 102
```

**Serial by design** — only one Schrödinger licence token is available, so parallel array
tasks would fail on licence checkout rather than finish sooner. The script is resumable
and appends results after every structure, so an interrupted run continues where it
stopped; just resubmit.

## Open issue: an unexplained irreproducibility

Analyses 02 and 05 produced different pose sets for **4D6Z** from inputs that differ only
by a title string (receptor coordinates, ligand coordinates, atom order and `gridgen.in`
all verified identical). Both Glide docking and grid generation were separately verified
deterministic — 3/3 and 2/2 identical repeats. The cause is **unexplained**; it is recorded
here rather than rationalised. The 22 structures shared with analysis 02 give a direct
empirical measure of how often this occurs.
