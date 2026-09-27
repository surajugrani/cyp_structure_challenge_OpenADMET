# 01 — Does the Glide score of crystal poses correlate with IC50?

**Question:** For CYP3A4 crystal structures with a measured IC50 for the bound ligand, how well does the Glide score of the *experimental* pose track binding affinity? This sets a baseline for how much to trust Glide scores when ranking poses for the challenge.

**Status:** Complete.

**Answer: no usable correlation.** Even with the experimentally correct pose, the Glide
score does not rank these 20 CYP3A4 ligands by potency (Spearman rho = −0.27, 95 % CI
[−0.72, +0.28], p = 0.26).

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


## Results

Main set, n = 20 (`results/correlations.csv`, `results/correlation_plot.png`):

| Score | Spearman rho [95 % CI] | p | Pearson r [95 % CI] | p |
|---|---|---|---|---|
| `inplace` docking score | **−0.265** [−0.72, +0.28] | 0.26 | **−0.363** [−0.74, +0.19] | 0.12 |
| `inplace` evdw | −0.134 [−0.63, +0.37] | 0.57 | −0.102 [−0.56, +0.35] | 0.67 |
| `inplace` ecoul | +0.153 [−0.32, +0.53] | 0.52 | +0.107 [−0.23, +0.47] | 0.65 |
| `mininplace` docking score | −0.093 [−0.59, +0.46] | 0.70 | −0.159 [−0.65, +0.37] | 0.50 |
| `mininplace` ecoul | −0.081 [−0.54, +0.42] | 0.73 | −0.129 [−0.47, +0.33] | 0.59 |
| `mininplace` evdw | −0.077 [−0.61, +0.53] | 0.75 | −0.113 [−0.61, +0.47] | 0.64 |

Bootstrap: 1000 resamples. Terms were examined in stages (`docking_score` + `evdw`, then
`ecoul`) rather than scanning all 13 at once — over 20 points, testing every term would
almost guarantee something crossed p < 0.05 by chance.

A more negative Glide score means better predicted binding, so a *negative* rho is the
expected direction. The sign is right but the interval spans zero.

**Sensitivity checks** (on `inplace_docking_score`) all agree — nothing rescues the
correlation, and nothing was inflating it:

| Subset | n | rho | 95 % CI |
|---|---|---|---|
| all (main set) | 20 | −0.265 | [−0.70, +0.25] |
| drop pIC50-spread-flagged | 18 | −0.166 | [−0.63, +0.39] |
| drop multi-copy (2V0M, 4K9T) | 18 | −0.298 | [−0.71, +0.25] |
| Fe-coordinating only (< 3 Å) | 16 | −0.256 | [−0.74, +0.35] |

### Why it fails: the score tracks molecular size

| Relationship | Spearman rho |
|---|---|
| MW vs Glide score | **−0.47** (bigger ligands score better) |
| MW vs pIC50 | 0.12 (size barely relates to potency here) |

Glide rewards size; size does not predict potency in this set. Two cases make it concrete:

- **08Y** (bromoergocryptine, MW 655) has the *best* score (−11.5) but only mid potency (pIC50 5.52).
- **CL6** (clotrimazole, MW 345) is the second most potent (pIC50 7.24) but scores −4.78.

### Glide's metal term is inactive, and that likely explains the failure

With the chemically correct heme (dianionic porphyrin, see *Structure preparation*), the
Glide **metal term is identically 0.00 for all 20 structures** — so it is excluded from
the table above; there is nothing to correlate.

This matters: **16 of the 20 ligands coordinate the heme iron at 1.9–2.4 Å**, and that
coordination is the defining interaction of a type-II CYP inhibitor. The total Glide score
gives it no credit at all.

The term is sensitive to the **iron's first coordination shell**, and to nothing else we
changed. 8SPD was scored under three heme states during development:

| 8SPD heme state | Propionates | Pyrrole N | Heme charge | Metal term |
|---|---|---|---|---|
| **A** — PrepWizard raw | COOH (neutral) | one N⁻, one N–H | +1 | −2.30 |
| **B** — after propionate fix | both −1 | one N⁻, one N–H | −1 | −2.30 |
| **C** — after `normalize_heme` | both −1 | both N⁻ | −2 | **0.00** |

- **A → B** changes only the propionates, 10 Å from the iron. Net heme charge moves by
  2 units and the metal term does not budge. Remote carboxylates are irrelevant to it,
  as expected for a ligand–metal interaction term.
- **B → C** changes only the pyrrole N–H. Net charge moves by 1 unit and the term goes
  to zero.

So this is not a net-charge effect (a 2-unit change did nothing, a 1-unit change zeroed
it) — it is specifically the iron's first coordination shell. The Fe bond count is
identical in all three states (4 × HEM:N + CYS:S), so it is the perceived valence or
charge state of the metal, not its connectivity.

A metalated porphyrin has all four nitrogens coordinating the iron and none protonated,
so **C is the correct chemistry** and is what this analysis uses. A plausible reading is
that the chemically wrong N–H makes Glide perceive the iron as coordinatively
unsatisfied, leaving a site to reward; completing the dianion makes it look satisfied.
That is an inference from three measurements on one structure, not something verified
inside Glide — but the practical consequence is solid either way: **as correctly
prepared, Glide scores these complexes with no metal-coordination reward at all.**

### Other observations

- **Electrostatics carry no signal either**, and `inplace_ecoul` even has the wrong sign
  (rho = +0.15).
- **`mininplace` is worse than `inplace`** (−0.09 vs −0.27). Minimization moves poses off
  their crystal geometry and loses information — consistent with the PXR finding that
  Glide post-processing degraded good AF3 poses.
- **Two crystal poses score as unfavourable** (positive docking score): 6AW (+0.64) and
  PK9 (+0.62, the 6.1 Å non-coordinating ligand). An experimentally observed pose scoring
  positive is a scoring-function failure, not a pose problem.

### Caveats

n = 20, with IC50s pooled from heterogeneous BindingDB assays (different substrates and
labs). This cannot distinguish "weak correlation" from "no correlation" — the CI is too
wide for that. What it does support is the practical conclusion: **Glide score is not a
usable affinity proxy for CYP3A4 here, even given a perfect pose.** Since docking a
predicted pose can only add error, this is an upper bound for score-based ranking in the
challenge.
