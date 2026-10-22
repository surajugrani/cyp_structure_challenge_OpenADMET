"""Fe-ligand distance dumbbell plots, written to figures/ as SVG.

    python3 tools/plot_fe_distances.py            # every configured plot
    python3 tools/plot_fe_distances.py --list     # show what is configured

Each plot is one row per ligand, sorted by the first series, with a bar joining the
series so disagreement is visible at a glance. The shaded band is the <=2.5 A
coordination cutoff; a degree sign marks ligands with no nitrogen donor, which should not
coordinate (analysis 19: 0 of 7 such crystal ligands do).

Adding a plot: append a PLOTS entry naming the series and where each one's distances come
from. Everything is read from committed CSVs and structure files, so the figures
regenerate from scratch after any analysis is rerun, and row order is deterministic.

For a ONE-OFF comparison, do not add a PLOTS entry and do not write a throwaway script -
point --series at the directories:

    python3 tools/plot_fe_distances.py \
        --series "23 ESMFold2 - 0.4527=submissions/23_esmfold2_default/submission" \
        --series "max accuracy=submissions/23_esmfold2_default/max_accuracy/raw" \
        --out fe_maxacc_vs_sub23 --title "max-accuracy against submission 23"

Each directory may hold .pdb or .cif; only ligands present in every series are drawn.
Distances come from tools/fe_distances.py, which caches on file mtime - a full
regeneration of all figures costs ~0.2 s once warm, against ~8 s cold.
"""
import argparse, csv, glob, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fe_distances as fd                                               # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
A16 = os.path.join(ROOT, "analyses", "16_fe_distance_swap", "results")
A18 = os.path.join(ROOT, "analyses", "18_glide_redock_predicted", "results")
S12 = os.path.join(ROOT, "submissions", "12_icm_redock_af3", "results")
COL = {"blue": "#2a78d6", "orange": "#eb6834", "teal": "#1baf7a", "yellow": "#eda100"}


def rows(path):
    return list(csv.DictReader(open(path)))


def donors():
    return {r["ligand"]: int(r["n_donors"]) > 0
            for r in rows(os.path.join(A16, "donor_recount.csv"))}


def src_fe_distance(col):
    return {r["ligand"]: float(r[col])
            for r in rows(os.path.join(A16, "fe_distance.csv")) if r[col]}


def src_donor_recount(col):
    return {r["ligand"]: float(r[col])
            for r in rows(os.path.join(A16, "donor_recount.csv")) if r[col]}


def src_glide_sel(col):
    return {r["ligand"]: float(r[col])
            for r in rows(os.path.join(A16, "glide_selection_fe.csv")) if r[col]}


def src_icm_sub12():
    return {r["ligand"]: float(r["fe_dist"])
            for r in rows(os.path.join(S12, "selection.csv")) if r["fe_dist"]}


def src_submission_fe(subdir):
    """Shipped submission PDBs."""
    return fd.series(os.path.join(ROOT, "submissions", subdir, "submission"))


def src_raw_sel(subdir):
    """A submission's SELECTED but not-yet-PrepWizard-ed structures.

    Only for a submission still in the queue - PrepWizard moves heavy atoms ~0.15 A
    median, so these are representative but not final.
    """
    return fd.series(os.path.join(ROOT, "submissions", subdir, "raw_sel"))


def src_esmfold2(sub="23_esmfold2_default", sub_dir="raw"):
    """ESMFold2 mmCIF output (remote Biohub API)."""
    return fd.series(os.path.join(ROOT, "submissions", sub, sub_dir))


def src_icm_glide_picked():
    """Fe distance of the ICM redocked pose that GLIDE ranks best (analysis 25).

    ICM's own Score ships a badly drifted pose for two ligands; this asks what the pool
    looks like when Glide does the ranking instead.
    """
    import csv as _csv
    f = os.path.join(ROOT, "analyses", "25_glide_on_icm_poses", "results",
                     "pick_comparison.csv")
    return {r["ligand"]: float(r["glide_fe"]) for r in _csv.DictReader(open(f))}


def src_icm_own_picked():
    import csv as _csv
    f = os.path.join(ROOT, "analyses", "25_glide_on_icm_poses", "results",
                     "pick_comparison.csv")
    return {r["ligand"]: float(r["icm_fe"]) for r in _csv.DictReader(open(f))}


def src_glide_redock():
    """Top-1 per ligand; poseviewer order is best-first, so take the first row seen."""
    out = {}
    for r in rows(os.path.join(A18, "glide_poses.csv")):
        if r.get("status") != "ok" or not r.get("fe_dist"):
            continue
        out.setdefault(r["ligand"], float(r["fe_dist"]))
    return out


def _sub17():
    """Final PDBs once submission 17 is packaged; selected raw structures until then."""
    if os.path.isdir(os.path.join(ROOT, "submissions", "17_boltz2_ipde", "submission")) \
            and glob.glob(os.path.join(ROOT, "submissions", "17_boltz2_ipde",
                                       "submission", "*.pdb")):
        return src_submission_fe("17_boltz2_ipde")
    return src_raw_sel("17_boltz2_ipde")


def _sub18():
    """Final PDBs once packaged; the selected raw structures until then."""
    d = os.path.join(ROOT, "submissions", "18_boltz2_glide", "submission")
    if os.path.isdir(d) and glob.glob(os.path.join(d, "*.pdb")):
        return src_submission_fe("18_boltz2_glide")
    return src_raw_sel("18_boltz2_glide")


PLOTS = [
    {"name": "fe_esmfold2_vs_best",
     "title": "ESMFold2 against the leading submissions",
     "series": [("06 AF3, PL-iPTM \u2014 0.4291", "blue",
                 lambda: src_submission_fe("06_af3_pliptm")),
                ("19 hybrid \u2014 0.4331  BEST", "teal",
                 lambda: src_submission_fe("19_hybrid_glide")),
                ("16 Boltz, PL-iPTM \u2014 0.4229", "yellow",
                 lambda: src_submission_fe("16_boltz2_pliptm")),
                ("ESMFold2 + ColabFold MSA \u2014 not submitted", "orange",
                 src_esmfold2)]},
    {"name": "fe_icm_glidepicked_vs_subs",
     "title": "ICM poses ranked by Glide, against submissions 06 and 19",
     "series": [("06 AF3, PL-iPTM \u2014 0.4291", "blue",
                 lambda: src_submission_fe("06_af3_pliptm")),
                ("19 hybrid \u2014 0.4331  BEST", "teal",
                 lambda: src_submission_fe("19_hybrid_glide")),
                ("ICM poses, Glide-ranked (analysis 25)", "orange",
                 src_icm_glide_picked)]},
    {"name": "fe_af3_icm_best",
     "title": "AF3 (06) vs ICM redock of it (12) vs best submission (19)",
     "series": [("06 AF3, PL-iPTM \u2014 0.4291 / 4.0647", "blue",
                 lambda: src_submission_fe("06_af3_pliptm")),
                ("12 ICM redocked into 06 \u2014 not submitted", "orange",
                 lambda: src_submission_fe("12_icm_redock_af3")),
                ("19 hybrid \u2014 0.4331 / 3.6754  BEST", "teal",
                 lambda: src_submission_fe("19_hybrid_glide"))]},
    {"name": "fe_hybrid_19_vs_20",
     "title": "Hybrids: pairwise Glide (19) vs union Glide (20)",
     "series": [("19 better-of-two \u2014 0.4331 / 3.6754", "blue",
                 lambda: src_submission_fe("19_hybrid_glide")),
                ("20 best of all \u2014 pending", "orange",
                 lambda: src_submission_fe("20_union_glide"))]},
    {"name": "fe_ipde_vs_glide",
     "title": "Boltz ensemble: complex_iPDE (17) vs Glide gscore (18)",
     "series": [("17 complex_iPDE \u2014 0.4157 / 3.8021", "teal",
                 lambda: src_submission_fe("17_boltz2_ipde")),
                ("18 Glide gscore \u2014 pending", "orange", _sub18)]},
    {"name": "fe_boltz_all_four",
     "title": "Boltz-2: defaults (03), tuned (15), PL-iPTM (16), iPDE (17)",
     "series": [("03 defaults \u2014 0.3921", "blue",
                 lambda: src_submission_fe("03_boltz2_prepwizard")),
                ("15 tuned, Boltz's own pick \u2014 0.4063", "orange",
                 lambda: src_submission_fe("15_boltz2_tuned")),
                ("16 tuned + PL-iPTM \u2014 0.4229", "teal",
                 lambda: src_submission_fe("16_boltz2_pliptm")),
                ("17 tuned + complex_iPDE \u2014 pending", "yellow", _sub17)]},
    {"name": "fe_boltz_progression",
     "title": "Boltz-2: defaults (03) vs tuned (15) vs tuned+PL-iPTM (16)",
     "series": [("03 Boltz defaults \u2014 0.3921", "blue",
                 lambda: src_submission_fe("03_boltz2_prepwizard")),
                ("15 tuned sampling \u2014 0.4063", "orange",
                 lambda: src_submission_fe("15_boltz2_tuned")),
                ("16 tuned + PL-iPTM \u2014 0.4229", "teal",
                 lambda: src_submission_fe("16_boltz2_pliptm"))]},
    {"name": "fe_af3_vs_boltz",
     "title": "AF3 (sub 06) vs Boltz-2 (sub 03)",
     "series": [("06 AF3, best PL-iPTM", "blue", lambda: src_fe_distance("af3_06_fe")),
                ("03 Boltz-2", "orange", lambda: src_fe_distance("boltz_03_fe"))]},
    {"name": "fe_pliptm_vs_glide_pick",
     "title": "PL-iPTM selection (sub 06) vs Glide-score selection (sub 13)",
     "series": [("06 AF3, best PL-iPTM — 0.4291", "blue",
                 lambda: src_glide_sel("pliptm_06_fe")),
                ("13 AF3, best Glide score — 0.3730", "orange",
                 lambda: src_glide_sel("glide_fe"))]},
    {"name": "fe_af3_vs_icm_redock",
     "title": "AF3 (sub 06) vs ICM redocked (sub 12)",
     "series": [("06 AF3, best PL-iPTM", "blue",
                 lambda: src_donor_recount("af3_fe_dist")),
                ("12 ICM redocked top-1", "orange", src_icm_sub12)]},
    {"name": "fe_three_way",
     "title": "AF3 vs Glide-picked vs Glide-redocked",
     "series": [("06 AF3, best PL-iPTM — 0.4291", "blue",
                 lambda: src_glide_sel("pliptm_06_fe")),
                ("13 AF3, best Glide score — 0.3730", "orange",
                 lambda: src_glide_sel("glide_fe")),
                ("Glide redocked top-1 (analysis 18)", "teal", src_glide_redock)]},
]


def draw(cfg, can):
    series = [(lab, COL[c], fn()) for lab, c, fn in cfg["series"]]
    # Ligand name breaks ties, so a rerun cannot reorder rows that share a distance.
    ligs = sorted(set.intersection(*[set(s[2]) for s in series]),
                  key=lambda l: (-series[0][2][l], l))
    if not ligs:
        return None
    vals = [v for _, _, m in series for l, v in m.items() if l in ligs]
    xmin, xmax = 1.6, max(6.3, max(vals) + 0.3)
    nlab = len(series)
    L, R = 140, 636
    TOP, ROW = 50 + 20 * ((nlab + 1) // 2) + 16, 21
    H = TOP + ROW * len(ligs) + 40

    def X(v):
        return L + (v - xmin) / (xmax - xmin) * (R - L)

    p = [f'<rect x="{X(xmin):.1f}" y="{TOP-12}" width="{X(2.5)-X(xmin):.1f}" '
         f'height="{ROW*len(ligs)+6}" fill="#1baf7a" opacity="0.08"/>',
         f'<text x="{X(2.5)-4:.1f}" y="{TOP-20}" text-anchor="end" class="ts">'
         f'coordinating, ≤2.5 Å</text>']
    v = 2
    while v <= xmax:
        p.append(f'<line x1="{X(v):.1f}" y1="{TOP-12}" x2="{X(v):.1f}" '
                 f'y2="{TOP+ROW*len(ligs)-6}" stroke="#9a9a94" stroke-opacity="0.3"/>')
        p.append(f'<text x="{X(v):.1f}" y="{TOP+ROW*len(ligs)+12}" '
                 f'text-anchor="middle" class="ts">{v}</text>')
        v += 1
    for i, lg in enumerate(ligs):
        y = TOP + i * ROW
        short = lg[5:] if lg.startswith("OCNT-") else lg
        p.append(f'<text x="{L-12}" y="{y+4}" text-anchor="end" class="ts">'
                 f'{short}{"" if can.get(lg, True) else " °"}</text>')
        vs = [m[lg] for _, _, m in series]
        if max(vs) - min(vs) > 0.02:
            p.append(f'<line x1="{X(min(vs)):.1f}" y1="{y}" x2="{X(max(vs)):.1f}" '
                     f'y2="{y}" stroke="#888780" stroke-width="1.3" stroke-opacity="0.6"/>')
        for (_, col, m) in reversed(series):
            p.append(f'<circle cx="{X(m[lg]):.1f}" cy="{y}" r="4.5" fill="{col}"/>')
    p.append(f'<line x1="{L}" y1="{TOP+ROW*len(ligs)-6}" x2="{R}" '
             f'y2="{TOP+ROW*len(ligs)-6}" stroke="#9a9a94"/>')
    p.append(f'<text x="{(L+R)//2}" y="{TOP+ROW*len(ligs)+32}" text-anchor="middle" '
             f'class="ts">minimum Fe–ligand heavy-atom distance (Å)</text>')
    for k, (lab, col, _) in enumerate(series):
        lx, ly = 146 + (k % 2) * 290, 20 + (k // 2) * 20
        p.append(f'<circle cx="{lx}" cy="{ly}" r="4.5" fill="{col}"/>')
        p.append(f'<text x="{lx+11}" y="{ly+4}" class="ts">{lab}</text>')
    nco = sum(1 for lg in ligs if not can.get(lg, True))
    p.append(f'<text x="146" y="{20+20*((nlab+1)//2)+4}" class="ts">'
             f'° = no nitrogen donor ({nco} of {len(ligs)}) · '
             + " · ".join(
                 f"{lab.split(',')[0]}: {sum(1 for lg in ligs if m[lg] <= 2.5)} coordinating"
                 for lab, _, m in series) + '</text>')
    desc = (f"Minimum haem iron to ligand distance for {len(ligs)} CYP3A4 ligands, "
            + " versus ".join(lab for lab, _, _ in series) + ".")
    return (f'<svg viewBox="0 0 700 {H}" xmlns="http://www.w3.org/2000/svg" role="img">'
            f'<title>{cfg["title"]}</title><desc>{desc}</desc>'
            + "".join(p) + '</svg>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--series", action="append", metavar="LABEL=DIR",
                    help="ad-hoc comparison: a directory of .pdb or .cif complexes, "
                         "repeatable, up to 4. Skips PLOTS entirely.")
    ap.add_argument("--out", metavar="NAME",
                    help="figure basename for --series (default: fe_adhoc)")
    ap.add_argument("--title", default="Ad-hoc Fe-distance comparison")
    a = ap.parse_args()
    if a.series:
        cols = ["blue", "orange", "teal", "yellow"]
        cfg = {"name": a.out or "fe_adhoc", "title": a.title, "series": []}
        for k, spec in enumerate(a.series):
            lab, _, d = spec.partition("=")
            if not d:
                raise SystemExit(f"--series needs LABEL=DIR, got {spec!r}")
            cfg["series"].append((lab, cols[k], lambda d=d: fd.series(d)))
        svg = draw(cfg, donors())
        if svg is None:
            raise SystemExit("no ligands common to every series")
        path = os.path.join(FIG, cfg["name"] + ".svg")
        open(path, "w").write(svg)
        print(f"{cfg['name']:28} -> {os.path.relpath(path, ROOT)}")
        return 0
    if a.list:
        for c in PLOTS:
            print(f"{c['name']:28} {c['title']}")
        return 0
    os.makedirs(FIG, exist_ok=True)
    can = donors()
    for c in PLOTS:
        try:
            svg = draw(c, can)
        except FileNotFoundError as e:
            print(f"{c['name']:28} skipped, missing input: {e.filename}")
            continue
        if svg is None:
            print(f"{c['name']:28} skipped, no overlapping ligands")
            continue
        path = os.path.join(FIG, c["name"] + ".svg")
        open(path, "w").write(svg)
        print(f"{c['name']:28} -> {os.path.relpath(path, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
