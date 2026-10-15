# Submission 15 — Boltz-2 with PXR's tuned sampling

**Result: 0.4063 LDDT-PLI / 4.0017 Å BiSyRMSD / 0.7983 LDDT-LP.**

**Best of any submission on BiSyRMSD and LDDT-LP.** Submission 06 still leads the primary
metric at 0.4291.

Identical to submission 03 except three Boltz flags, so the comparison isolates the
sampling settings:

| | 03 default | 15 tuned | change |
|---|---|---|---|
| LDDT-PLI | 0.3921 | 0.4063 | **+0.0142** |
| BiSyRMSD | 4.1179 | **4.0017** | **−0.116 Å** |
| LDDT-LP | 0.7962 | **0.7983** | +0.0021 |

## What changed

```
--recycling_steps 5     (default 3)
--sampling_steps  300   (default 200)
--diffusion_samples 5   (default 1)
```

Same YAMLs, same MSA server, same GPU recipe, same PrepWizard options.

## Why these values

From the PXR challenge, where the same knobs were tuned:

| PXR run | LDDT-PLI |
|---|---|
| 8 ColabFold MSA, defaults | 0.4439 |
| **9 recycling 5, sampling 300** | **0.4621** (+0.0182) |
| 10 recycling 6, sampling 400 | 0.4548 — more is *not* better |
| 11 recycling 5, sampling 300 + 3 templates | 0.4516 |

5/300 was the sweet spot there; 6/400 gave back a third of the gain. **The transfer held:
+0.0142 here against +0.0182 on PXR.**

`--diffusion_samples 5` goes beyond what PXR varied. Analysis 21 found Boltz's advantage
over AF3 is that it *adapts* the pocket per ligand (AF3 matched crystal 6DAA for 15 of 20
ligands; Boltz spread across five crystals), so a diffusion ensemble should give real
conformational spread. We cannot yet pick the best member — analyses 09, 10, 20 and
submissions 13/14 all failed at exactly that — so this submission simply takes **Boltz's
own top-ranked sample** (`model_0`).

## Boltz's own ranking

Boltz ranks by `confidence_score`, a blend of complex pLDDT and ipTM, and writes the
models already sorted, so `model_0` is its best pick (verified on 8 of 8 spot-checked
ligands).

**The spread is tiny** — across 5 samples `confidence_score` varies by only 0.001–0.007
(e.g. 0.949 / 0.949 / 0.948 / 0.948 / 0.948). That is the same flatness that made AF3's
iPTM near-useless as a selector, so `model_0` is close to arbitrary within the ensemble.

Worth a future attempt: the confidence JSON also carries `pair_chains_iptm[0][2]`, Boltz's
protein–ligand iPTM and the direct analogue of the metric that won for AF3 (submission
06). It reads ~0.70 where `confidence_score` reads ~0.95, so it has far more dynamic
range. The other four samples per ligand are kept on disk for when a working selector
exists.

## Build note — the mistake that cost a rerun

**Boltz names every het chain `LIG`, including the haem.** Submission 01 fixed this in
`03_package.py` by renaming chain B to `HEM`; feeding raw Boltz output straight into
PrepWizard skipped that step, and the whole set failed validation 0/20 with
`all_atoms_connected`, `minimum_distance_to_protein` and `Fe nan->nan` — the validator saw
two `LIG` residues and could not find the iron. Renaming the 43 chain-B atoms per
structure before PrepWizard fixed it.

| file | contents |
|---|---|
| `raw_model0/` | Boltz `model_0` per ligand, haem renamed to HEM |
| `results/` | full Boltz output — all 5 diffusion samples, confidence JSONs, PAE/pLDDT |
| `work/` | PrepWizard output (git-ignored) |
