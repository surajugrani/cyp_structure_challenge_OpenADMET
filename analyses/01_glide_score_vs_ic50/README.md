# 01 — Does the Glide score of crystal poses correlate with IC50?

**Question:** For CYP3A4 crystal structures with a measured IC50 for the bound ligand, how well does the Glide score of the *experimental* pose track binding affinity? This sets a baseline for how much to trust Glide scores when ranking poses for the challenge.

**Status:** Steps 1–4 implemented and smoke-tested on one structure (8SPD). Full run and step 5 pending.

## Pipeline

| Step | Script | Runs on | Output |
|---|---|---|---|
| 1 | `scripts/01_fetch_affinities.py` | login node | `data/rcsb_raw.json`, `data/affinities_all.csv`, `data/entries.csv` |
| 2 | `scripts/02_curate.py` | login node | `data/ic50_values_used.csv`, `data/curated_set.csv` |
| 3 | `scripts/03_prep_structures.py` | sbatch | `structures/<PDB>/` (git-ignored) |
| 4 | `scripts/04_grid_and_score.py` | sbatch | Glide `inplace` + `mininplace` scores |
| 5 | `scripts/05_correlate.py` | login node | `results/` table + plots |

Steps 3–4 are CPU-heavy and **must not be run on the login node**. Smoke-test a single
structure first, then submit the full set:

```bash
sbatch scripts/test_one.sh 8SPD      # one structure, ~5 min
sbatch scripts/run_prep_and_score.sh # all 22
squeue -u $USER                      # check progress
```

Run from this folder with Schrödinger's Python (it already has pandas, scipy and requests):

```bash
module load schrodinger/2025-3
$SCHRODINGER/run python3 scripts/01_fetch_affinities.py
$SCHRODINGER/run python3 scripts/02_curate.py
```

## Data source

The starting point is an RCSB custom report (`input/rcsb_pdb_custom_report_20260926222041.csv`, 122 CYP3A4 entries). Only the **entry IDs** are taken from it. Affinities are fetched from the RCSB Data API (GraphQL), which ties each value to a specific ligand (`comp_id`) and its source database. All 70 values came from BindingDB.

## Curation decisions

| Decision | Reason |
|---|---|
| IC50 only; drop Kd/Ki | Different quantities; 58 of 70 values are IC50 |
| Drop qualified values (`<`, `>`) | They are bounds, not measurements (2 dropped) |
| Drop Ir-containing complexes | Most docking programs have no Ir parameters (12 entries in the report, none had IC50) |
| Reference sequence identity ≥ 95 % | All passed (98–100 %) |
| pIC50 = −log10(IC50 [M]), **median** per ligand | Log scale, so a single 80 µM value can't dominate; median resists outlier assays |
| One structure per ligand: **best resolution** | Repeat structures of one ligand aren't independent points |
| Ligand copy: chain A (lowest chain ID), lowest resnum | Consistent and deterministic |

**Flags kept for review (not removed automatically):**
- `spread_flag`: pIC50 spread > 1.5 log units. **KLN/2V0M** (ketoconazole, 18 values, spread 3.8) and **MWV/6OO9** (mibefradil, 5 values, spread 3.0). The medians look reasonable, but the spread shows how assay-dependent CYP3A4 IC50s are.
- `shared_value`: **G0V/6DAB, FZV/6DAJ, G0Y/6DAG** are different stereoisomers that all have exactly 310 nM in BindingDB, most likely one measurement attached to all of them. Only the best-resolution one (G0V/6DAB) is in the main set.
- `lig_copies_in_chain` = 2: **2V0M** (two stacked ketoconazole molecules in the pocket, a known feature of this structure) and **4K9T**. Step 3 must decide how to handle the second copy.

## Structure preparation (step 3)

Per structure: keep chain A (protein + HEM + every copy of the ligand of interest);
drop waters, cryoprotectants and other chains. The **scoring copy is the ligand copy
closest to the heme Fe**, chosen by measured distance rather than residue number.

Measured Fe–ligand distances confirm the choice: most ligands coordinate the iron at
1.9–2.4 Å (classic type-II binding). Four sit further out: 6OO9 (3.76 Å), 7LXL (3.73 Å),
3UA1 (3.65 Å) and 4D6Z (6.09 Å) — the last being the weakest binder in the set.

PrepWizard: `-rehtreat -captermini -fillsidechains -propka_pH 7.4 -epik_pH 7.4
-f OPLS_2005 -rmsd 0.3` (falls back without `-fillsidechains` if Prime fails).

**Post-prep fix:** PrepWizard leaves the two heme propionates as neutral COOH. At pH 7.4
they are carboxylates, so the script removes the acidic H and sets each oxygen to −1
(heme formal charge goes +1 → −1). This is systematic across all structures so it barely
affects the *correlation*, but it makes the electrostatics defensible.

**Checks recorded per structure:** ligand atom count, ligand and heme formal charge, and
the atoms bonded to Fe (expected: 4 × HEM:N + CYS:S, plus the ligand N for type-II binders).

## Glide scoring (step 4)

Grid centred on the crystal ligand (inner box 10 Å, outer box derived from ligand size),
then two scoring modes on the same pose:

- `inplace` — scores the crystal coordinates as-is
- `mininplace` — local minimization within the site before scoring

Scores and all component terms (lipo, hbond, **metal**, evdw, ecoul, esite, emodel) are
read from the Glide results CSV.

### Version notes (schrodinger/2025-3)

Three things differ from the older scripts used in the PXR challenge:

| Issue | Resolution |
|---|---|
| `prepwizard -f 2005` is rejected | The flag now takes the force field *name*: `-f OPLS_2005` |
| `generate_glide_grids -j NAME` | Writes `NAME-gridgen.zip`, not `NAME.zip` |
| No `.rept` file produced | Glide writes a results **CSV**; parsing that is simpler and more robust |

## Result of curation

22 unique ligands; **20 in the main set**, pIC50 range 4.1–7.4. See `data/curated_set.csv`.
