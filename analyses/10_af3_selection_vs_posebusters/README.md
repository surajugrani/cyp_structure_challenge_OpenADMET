# 10 — Which AF3 metric picks a PoseBusters-valid pose?

**Question:** AF3 emits 5 samples per ligand. Submission 04 selected on protein-ligand
iPTM and got only 13/20 through PoseBusters. Does any reported confidence metric predict
physical validity, so selection alone could raise that?

**Answer: no metric does — but PoseBusters itself can, and it reaches the ceiling.**

## Method

All 100 samples (20 ligands x 5) converted to submission-style PDB and run through
`PoseBusters(config="dock")` exactly as the organizers do, then joined to every metric in
AF3's `summary_confidences.json`.

## Results

**57/100 samples pass PoseBusters.** 19/20 ligands have at least one passing sample, so
**19/20 is the ceiling for any selection rule.**

| selection rule | PB passes |
|---|---|
| random sample | ~11.4/20 (57 %) |
| best `pl_iptm` (what submission 04 used) | 13/20 |
| best `iptm` | 12/20 |
| best `ptm` | 13/20 |
| best `ranking_score` (AF3's own) | 13/20 |
| best `chain_iptm_lig` | 14/20 |
| lowest `has_clash` | 13/20 |
| lowest `fraction_disordered` | 13/20 |
| **PoseBusters-pass, then best `pl_iptm`** | **19/20** |

### AF3's confidences carry no signal about physical validity

The best metric beats random by 2.6 structures out of 20. Three candidates are
**constant across all 100 samples** and cannot discriminate at all:

| metric | distinct values in 100 samples |
|---|---|
| `has_clash` | **1** (0.0 everywhere) |
| `fraction_disordered` | **1** (0.03) |
| `ptm` | 2 (0.91–0.92) |
| `iptm` | 9 (0.89–0.97) |
| `pl_iptm` | 20 (0.76–0.95) |

`has_clash` being identically zero while 43 % of samples fail a PoseBusters distance check
is the clearest statement of the problem: AF3 does not know when its own poses are
sterically invalid.

### Use PoseBusters as the selector

It needs no ground truth and is the organizers' own gate, so filtering on it is legitimate
rather than leakage. Selecting the highest-`pl_iptm` sample *among those that pass* gives
**19/20**, changing the chosen sample for 6 ligands.

**OCNT-2535519 is unreachable this way** — all 5 of its samples fail
`minimum_distance_to_protein` (it coordinates Fe at 2.03 Å). That one needs PrepWizard.

### Failure modes over all 100 samples

| count | check |
|---|---|
| 41 | `minimum_distance_to_protein` |
| 3 | `volume_overlap_with_protein` |
| 2 | `minimum_distance_to_organic_cofactors` |

Essentially a single problem: contacts that are too close.

## Caveat

This optimises **physical validity, not accuracy.** PoseBusters says a pose is not absurd,
not that it is right. Filtering on it protects the Coverage term and avoids zeroed scores,
but a PB-valid pose in the wrong orientation still scores badly on LDDT-PLI. Whether
PB-filtered selection actually improves the leaderboard is an empirical question the next
submission answers.

## Rerun

```bash
sbatch scripts/run.sh
```
