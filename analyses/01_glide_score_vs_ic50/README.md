# 01 — Does the Glide score of crystal poses correlate with IC50?

**Question:** For CYP3A4 crystal structures with a measured IC50 for the bound ligand, how well does the Glide score of the *experimental* pose track binding affinity? This sets a baseline for how much to trust Glide scores when ranking poses for the challenge.

**Status:** Steps 1–2 (data curation) done. Steps 3–5 (Glide) pending.

## Pipeline

| Step | Script | Runs on | Output |
|---|---|---|---|
| 1 | `scripts/01_fetch_affinities.py` | login node | `data/rcsb_raw.json`, `data/affinities_all.csv`, `data/entries.csv` |
| 2 | `scripts/02_curate.py` | login node | `data/ic50_values_used.csv`, `data/curated_set.csv` |
| 3 | `scripts/03_prep_structures.py` | sbatch | `structures/<PDB>/` (git-ignored) |
| 4 | `scripts/04_grid_and_score.py` | sbatch | Glide `inplace` + `mininplace` scores |
| 5 | `scripts/05_correlate.py` | login node | `results/` table + plots |

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

## Result of curation

22 unique ligands; **20 in the main set**, pIC50 range 4.1–7.4. See `data/curated_set.csv`.
