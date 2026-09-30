# 03 — Boltz-2 + PrepWizard restrained minimisation

**Status:** submitted. 20/20 structural checks, 20/20 PoseBusters.
**Leaderboard: LDDT-PLI 0.3921, BiSyRMSD 4.1179 Å, LDDT-LP 0.7962, Coverage 1 — rank 5/9.**
See [`../LEADERBOARD.md`](../LEADERBOARD.md).

Fixes the three `minimum_distance_to_protein` failures in
[submission 01](../01_boltz2_default/) without disturbing the poses.

## Method

PrepWizard on each of the 20 Boltz-2 predictions, restrained to 0.3 Å global heavy-atom
RMSD (`-rmsd 0.3`, analysis 01's proven option set — `-f OPLS_2005`, not `-f 2005`, which
schrodinger/2025-3 rejects). 0.3 Å is a *global* threshold over ~4100 atoms, so clashing
side-chain tips can move the ~0.5 Å needed while the average stays tight; 0.5 Å was not
needed.

Minimised coordinates are **grafted** onto the Boltz PDBs — only the x/y/z columns of
matching heavy atoms are replaced, keyed on (chain, resnum, atom name). PrepWizard adds
hydrogens and rewrites the file, and the validator is sensitive to atom naming and serial
uniqueness; the Boltz files already pass all of that.

## Results

| | before | after |
|---|---|---|
| OCNT-2395328 | 2.09 Å | **2.81 Å** |
| OCNT-2395342 | 2.45 Å | **2.89 Å** |
| OCNT-2395712 | 2.36 Å | **2.78 Å** |

Displacement 0.16–0.18 Å median across all 20. Fe coordination preserved: the 11
coordinating ligands went from 2.01–2.27 Å to 2.09–2.33 Å, so the concern that PrepWizard
would pull ligands off the iron did not materialise. For those 11 the closest protein
contact *is* the Fe — type-II binding — and PoseBusters passes them, confirming it does
not penalise metal coordination.

`submission_boltz2_prepwizard.zip` — this is the submission candidate.

## Rejected alternative: ICM self-tethered minimisation

Attempted first (see `../02_boltz2_icm_tethered/`) and abandoned. ICM's own recipe is
`convertObject` (which stores each atom's input coordinates as a **selftether**) followed
by `minimize cartesian ... selftether=as_` with the **`ts`** term — not `tz`, which is for
tethers between two objects. Six configurations failed before the documentation
(`man/icm/icm-commands.html`, `man/icm/terms.html`) was read properly. Even once correct it
moved heavy atoms 1.0–7.6 Å median and corrupted 6 ligands, versus PrepWizard's 0.17 Å.
