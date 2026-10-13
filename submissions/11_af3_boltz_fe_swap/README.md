# Submission 11 — submission 06 with 3 structures taken from Boltz-2

Submission 06 (AF3, PL-iPTM selection, our best at 0.4291) with three structures replaced
by submission 03's Boltz-2 predictions:

| ligand | AF3 (06) Fe dist | Boltz (03) Fe dist | Boltz closer by |
|---|---|---|---|
| OCNT-2395645 | 4.66 Å | 4.02 Å | 0.63 Å |
| OCNT-2395460 | 4.43 Å | 3.52 Å | 0.91 Å |
| OCNT-2313429 | 4.33 Å | 3.86 Å | 0.47 Å |

17/20 from AF3, 3/20 from Boltz. 20/20 structural, 20/20 PoseBusters.

## Rationale and its limits

All three sit in the non-coordinating (type I) cluster, and the working assumption is that
the pose nearer the catalytic iron is the likelier one. This is a deliberately crude test,
and the caveats are real:

* **"Closer to Fe is better" is an assumption.** For a type I binder there is no reason
  the true distance should be minimal; these ligands legitimately sit 3.5–5 Å out and the
  correct value for any one of them is unknown.
* **[Analysis 16](../../analyses/16_fe_distance_swap/) found no pocket failures to fix.**
  Nothing in either submission lies beyond 4.7 Å of the iron, and AF3 coordinates *more*
  ligands than Boltz (12 vs 11). These swaps are refinements within the right pocket, not
  rescues.
* **The likeliest outcome is no change at all.** With ~10 of 20 complexes scored, three
  swapped structures may include zero that are actually scored — in which case the result
  is identical to 06 and uninformative, rather than a negative result.

## One part is not a gamble

Boltz's protein is better: submission 03 scored **LDDT-LP 0.7962** against 06's 0.6985.
Taking 3 of 20 structures from Boltz should raise LDDT-LP by roughly 3/20 × 0.098 =
**+0.015** fairly reliably, because that is a backbone-quality difference rather than a
per-pose bet. LDDT-PLI is the part that could go either way.

## Deliberately NOT swapped: OCNT-2535422

The largest AF3/Boltz disagreement (2.28 vs 3.80 Å, a 1.52 Å gap) is excluded on purpose.
Its swap is justified by **chemistry** — no imidazole, pyridine or triazole nitrogen, so
it should not be a type II binder and Boltz's non-coordinating pose looks right — whereas
the three above are justified by the **distance heuristic**. Mixing them would confound
the two: if 11 improves, we could not tell which reasoning worked. See analysis 16.

| file | contents |
|---|---|
| `results/sources.csv` | per structure: which submission it came from, heavy-atom count |
