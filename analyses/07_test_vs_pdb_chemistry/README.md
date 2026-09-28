# 07 — How similar are the blinded challenge ligands to the CYP3A4 PDB ligands?

**Question:** Analyses 05/06 measure docking performance on ~97 CYP3A4 crystal structures.
Does that benchmark actually tell us what to expect on the 20 blinded compounds?

**Status:** Complete.

**Answer: the two sets barely overlap — but the difference probably works in our favour.**

## 1. Chemically distinct

| | Blinded (20) | PDB (92) | Mann-Whitney p |
|---|---|---|---|
| Max Tanimoto to any PDB ligand | **median 0.20**, max 0.34 | — | — |
| Shared Murcko scaffolds | **0 of 18** | — | — |
| MW | 329 | 531 | 3 × 10⁻⁸ |
| cLogP | 2.8 | 5.2 | 1 × 10⁻⁸ |
| Rotatable bonds | **4.5** | **12** | 7 × 10⁻⁷ |
| Heavy atoms | 23.5 | 37 | 1 × 10⁻⁸ |
| TPSA | 58 | 80 | 0.002 |

Morgan r=2, 2048-bit fingerprints. **Not one blinded ligand reaches Tanimoto 0.4** to
anything in the benchmark, and there are no shared scaffolds. The blinded compounds look
like a fragment-to-lead series; the PDB set is dominated by large drugs (ritonavir,
ketoconazole, ergot alkaloids).

## 2. Why that is probably good news

Docking accuracy degrades sharply with ligand flexibility, and our own data show it
clearly (n = 97, Glide):

| Rotatable bonds | n | top-1 ≤ 2 Å | median RMSD |
|---|---|---|---|
| 0–4 | 15 | **27 %** | 3.5 Å |
| 5–8 | 15 | **40 %** | 3.5 Å |
| 9–12 | 30 | 6.7 % | 6.0 Å |
| 13+ | 37 | 5.4 % | 4.8 Å |

Flexibility vs best-of-N RMSD: rho = +0.44 (p = 7 × 10⁻⁶); vs top-1 rho = +0.29 (p = 0.004).
Ligands with ≤ 8 rotatable bonds succeed at **33 %** (10/30) versus **6 %** (4/67) above that.

**Every blinded ligand has ≤ 7 rotatable bonds** (10 in the 0–4 bin, 10 in 5–8, none in the
hard bins). Re-weighting the benchmark's per-bin success rates by the blinded set's
flexibility profile gives a **flexibility-matched expectation of ~33 % top-1** (≈7 of 20),
against the 13 % measured on the benchmark as a whole.

So the headline number from analyses 05/06 is likely **pessimistic for the actual
challenge, by roughly 2.5×**.

## 3. Fe-coordinating motifs

Analyses 01–04 showed the benchmark is dominated by type-II inhibitors coordinating the
heme iron, and that this is exactly what Glide's scoring mishandles. Motif prevalence:

| Motif | Blinded | PDB |
|---|---|---|
| any Fe-coordinating motif | **60 %** | 90 % |
| pyridine | 30 % | 64 % |
| imidazole | 25 % | 13 % |
| thiazole | 0 % | 12 % |
| nitrile | 15 % | 0 % |

**40 % of the blinded ligands carry no obvious coordinating motif**, versus 10 % of the
benchmark. Those compounds must bind in the lipophilic pocket without an iron anchor —
the regime where our benchmark performed worst (4D6Z, the one clearly non-coordinating
case in the affinity subset, had the worst pose recovery of all).

## Caveats

- **This is extrapolation.** With zero scaffold overlap and max Tanimoto 0.34, the blinded
  set is outside the benchmark's applicability domain. The flexibility-matched estimate
  corrects for one axis only.
- Bin sizes are small (15 and 15 in the two relevant bins), so ±10 % on that 33 % would not
  be surprising.
- The correction does **not** account for the lower Fe-motif prevalence, which pushes the
  other way.

## Implication for the challenge

Do not read 13 % as the expected outcome. The blinded compounds are small and rigid, which
is the regime where docking works best in our own measurements. But 40 % of them lack an
iron anchor, and our heme-related findings (analyses 01–04) apply most directly to the
type-II compounds. A pose-selection method that does not rely on the Glide score — see
[[pose-selector-idea]] — remains the highest-value direction.
