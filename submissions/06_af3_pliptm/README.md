# 06 — AF3 + PrepWizard, selected on protein–ligand iPTM alone

**Status:** ready. 20/20 structural, 20/20 PoseBusters.
`submission_af3_pliptm.zip`

A controlled A/B against [05](../05_af3_prepwizard/). Identical in every respect except
the selection key: `chain_pair_iptm[0][2]` alone, instead of averaged with the
haem–ligand term. **Exactly 2 of 20 structures differ**, both reverting to
`seed-1_sample-0`:

| ligand | 05 chose (mean) | 06 chooses (PL only) | PL | HL |
|---|---|---|---|---|
| OCNT-2312554 | sample-2 | sample-0 | 0.88 | 0.68 |
| OCNT-2535422 | sample-1 | sample-0 | 0.83 | 0.53 |

In both, `pl_iptm` was tied at 2 dp, so the haem–ligand term was breaking a tie rather
than overriding a preference.

**Caveat on interpreting the result:** the leaderboard scores ~10 of the 20 complexes. If
neither differing structure is in that subset, the metrics will be identical — which would
show nothing either way, not that the rules are equivalent. Roughly a 1-in-4 chance both
are excluded.
