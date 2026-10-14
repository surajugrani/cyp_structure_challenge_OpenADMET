# figures/

Regenerate everything with:

```bash
python3 tools/plot_fe_distances.py        # --list to see what is configured
```

All plots are rebuilt from committed CSVs, so rerunning an analysis and then this script
keeps the figures in step. They are SVG: open in a browser, or drop into a slide.

| file | shows |
|---|---|
| `fe_af3_vs_boltz.svg` | submission 06 (AF3) against submission 03 (Boltz-2) |
| `fe_pliptm_vs_glide_pick.svg` | PL-iPTM selection (06, 0.4291) against Glide-score selection (13, 0.3730) |
| `fe_af3_vs_icm_redock.svg` | submission 06 against ICM redocked top-1 (submission 12) |
| `fe_three_way.svg` | 06, Glide-picked and Glide-redocked together |

## Reading them

One row per ligand, sorted by the first series. The bar joins the series, so a long bar
means the methods disagree. The shaded band is the <=2.5 A coordination cutoff.

A degree sign after the ligand name means **no nitrogen donor** - that ligand should not
coordinate the iron. The classification comes from analysis 16 and is broader than the
classic type II motif: it counts aromatic N, aliphatic amines, aniline/indoline N and
nitriles, not just pyridine/imidazole/triazole. An earlier aromatic-only version wrongly
flagged three ligands, OCNT-2535422 among them.

Base rate from analysis 19, for judging a long distance: in 81 CYP3A4 crystal structures,
60 of 74 Fe-capable ligands (81 %) coordinate, and **0 of 7** without a donor do. So
coordinating *without* a donor is a strong error signal, while failing to coordinate *with*
one is only weak evidence - and part of that 19 % is structures where a second small
ligand occupies the axial site and blocks access.

## Adding a plot

Append an entry to `PLOTS` in `tools/plot_fe_distances.py`, naming each series and the
function that supplies its per-ligand distances.
