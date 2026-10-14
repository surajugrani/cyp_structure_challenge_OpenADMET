# Analysis 21 — why does Boltz beat AF3 on LDDT-LP?

**Because AF3 predicts essentially the same pocket every time, while Boltz adapts it to
the ligand.**

The leaderboard reports Boltz-2 (submission 03) at LDDT-LP 0.7962 against AF3
(submission 06) at 0.6985 — a 0.098 gap on the protein pocket — with no per-structure
breakdown to explain it.

## Method

The organizers' LDDT-LP cannot be reproduced without the blinded references. But CYP3A4
is the same protein in every PDB entry, so a crystal pocket is a legitimate yardstick for
"does this look like a real CYP3A4 conformation".

* **References:** ligand-bound crystals at ≤2.5 Å, 38 of the 100 in analysis 05. Apo
  pockets are a different shape; poor maps add noise.
* **Every model against every reference**, not one chosen reference — the CYP3A4 pocket is
  famously plastic, so a single reference would bake in its own conformation. Reported as
  `max` (best-matching crystal: does this resemble *any* real conformation?) and `median`
  (how typical is it across the ensemble).
* **Fixed pocket:** residues with a heavy atom within 10 Å of the haem iron in the
  highest-resolution reference. Fixed rather than per-model, so both methods are scored on
  identical atoms and the comparison is not confounded by each model's ligand placement.
* **lDDT from OpenStructure** — the organizers' own tool, in the `pxr-eval` environment.

## Result

| | max lDDT (mean) | median-vs-ensemble (mean) |
|---|---|---|
| Boltz-2 (03) | **0.9726** | **0.9513** |
| AF3 (06) | 0.9650 | 0.8955 |

Boltz wins on **19 of 20 ligands**; the one exception, OCNT-2395789, is +0.002.

### The gap is in the median, not the max

0.056 on the median against only 0.008 on the max. Both models match *some* crystal
closely — but Boltz's pocket resembles the whole ensemble, while AF3's fits one
conformation and drifts from the others.

### The best-match pattern says the same thing outright

| model | best-matching crystals |
|---|---|
| AF3 (06) | **6DAA ×15**, 6UNG ×3, 6DA5 ×1, 6DA3 ×1 |
| Boltz (03) | 7UFB ×4, 4D75 ×4, 6BCZ ×4, 6DAJ ×2, 6OOB ×2 |

**AF3 lands on 6DAA for 15 of 20 ligands** — the same pocket conformation almost
regardless of what it is binding. Boltz spreads across at least five distinct crystals.
AF3 is producing a near ligand-independent pocket; Boltz is varying conformation with the
ligand.

## Why that would cost LDDT-LP

CYP3A4's binding site is one of the most plastic in the proteome, and the blinded
references will span several conformational states. A model locked to one state scores
well on the few references that happen to match it and poorly on the rest; a model that
adapts tracks them individually. That is exactly the asymmetry measured here, and it is a
plausible mechanism for the 0.098 leaderboard gap.

It also fits an earlier observation that looked contradictory: analysis 17 found AF3's
receptor is the *better docking target* despite the worse LDDT-LP. A single
well-formed, induced-fit pocket can be easy to dock into while still being the wrong
conformation for most of the reference set.

## Caveat on magnitude

These values (0.96–0.97) sit far above the leaderboard's 0.70–0.80 because the pocket
here is a 10 Å heme shell on a fixed residue set, not the organizers' ligand-defined
binding site, and the references are different structures entirely. **The direction and
the consistency are the result — not the absolute numbers.**

| file | contents |
|---|---|
| `results/pocket_lddt.csv` | per model and ligand: max / median / min lDDT, number of references, best-matching PDB and its resolution |
| `results/pocket_residues.json` | the fixed pocket residue numbers |
