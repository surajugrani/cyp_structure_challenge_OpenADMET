"""
Analysis 04 step 2: compare the metal-term-activating heme (fe2_nh) against the
chemically correct baseline (fe2_dianion, analyses 01/02).

Three questions, in order of how much they matter:

  1. Pose recovery  - does activating the metal term help Glide find the crystal pose?
  2. Affinity rank  - does it help the score track pIC50?
  3. Attribution    - is any change due to the metal term or to the residual ecoul shift?

Plus a characterisation the baseline could not provide: **metal term vs Fe distance**
across every redocked pose. In the baseline the term is identically 0, so it has no
variance; here it varies, and analysis 03's smoke test hinted the geometry window may
favour longer contacts than a real dative bond (8SPD: the near-native 1.73 A pose scored
metal = 0.00 while a displaced 3.96 A pose scored -0.79). If that holds across the set,
the term would be actively mis-ranking type-II poses rather than merely missing them.

Usage (light - no docking):
    $SCHRODINGER/run python3 scripts/02_compare.py

Outputs: results/comparison.csv, results/metal_vs_distance.csv
"""

import os

import numpy as np
import pandas as pd
from scipy import stats
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
A02 = os.path.join(ANALYSES, "02_redocking_pose_recovery")
RESULTS = os.path.join(HERE, "results")
STRUCT = os.path.join(HERE, "structures")

N_BOOT = 1000
RNG = np.random.default_rng(0)


def boot_ci(x, y, n_boot=N_BOOT):
    out = []
    for _ in range(n_boot):
        i = RNG.integers(0, len(x), len(x))
        if len(np.unique(x[i])) < 3 or len(np.unique(y[i])) < 3:
            continue
        out.append(stats.spearmanr(x[i], y[i]).statistic)
    return (np.nan, np.nan) if not out else tuple(np.percentile(out, [2.5, 97.5]))


def recovery(df, col="rmsd_top1", n=None):
    n = n or len(df)
    return {f"<=2A": int((df[col] <= 2).sum()), f"<=1A": int((df[col] <= 1).sum()),
            "median": float(df[col].median()), "n": n}


def main():
    cur = pd.read_csv(os.path.join(A01, "data", "curated_set.csv"))

    # ── 1 & 2: baseline vs wrong chemistry ──────────────────────────────────
    base_score = pd.read_csv(os.path.join(A01, "results", "glide_scores.csv"))
    base_pose = pd.read_csv(os.path.join(A02, "results", "pose_rmsd.csv"))
    new_score = pd.read_csv(os.path.join(RESULTS, "inplace_scores.csv"))
    new_pose = pd.read_csv(os.path.join(RESULTS, "redock_summary.csv"))

    rows = []
    print("=== Pose recovery ===")
    for label, d in (("baseline (fe2_dianion)", base_pose), ("wrong chem (fe2_nh)", new_pose)):
        d = d[d["rmsd_top1"].notna()]
        for metric, col in (("top-1", "rmsd_top1"), ("best-of-N", "rmsd_best")):
            r = recovery(d, col)
            print(f"  {label:<24} {metric:<10} <=2A {r['<=2A']:2d}/{r['n']} "
                  f"({100*r['<=2A']/r['n']:3.0f}%)  <=1A {r['<=1A']:2d}  median {r['median']:.2f}")
            rows.append({"comparison": "pose_recovery", "set": label, "metric": metric, **r})

    print("\n=== Affinity correlation (main set) ===")
    main_ids = set(cur.loc[cur["in_main_set"], "pdb_id"])
    for label, d, col in (("baseline (fe2_dianion)", base_score, "inplace_docking_score"),
                          ("wrong chem (fe2_nh)", new_score, "inplace_docking_score")):
        m = cur.merge(d, on=["pdb_id", "comp_id"]).query("in_main_set")
        sub = m[[col, "pIC50_median"]].dropna()
        x, y = sub[col].to_numpy(), sub["pIC50_median"].to_numpy()
        rho = stats.spearmanr(x, y)
        lo, hi = boot_ci(x, y)
        pr = stats.pearsonr(x, y)
        print(f"  {label:<24} n={len(sub):2d}  rho={rho.statistic:+.3f} [{lo:+.2f},{hi:+.2f}] "
              f"p={rho.pvalue:.3f}   r={pr.statistic:+.3f} p={pr.pvalue:.3f}")
        rows.append({"comparison": "affinity", "set": label, "n": len(sub),
                     "spearman_rho": rho.statistic, "spearman_lo": lo, "spearman_hi": hi,
                     "spearman_p": rho.pvalue, "pearson_r": pr.statistic})

    # ── 3: attribution - metal vs ecoul ─────────────────────────────────────
    print("\n=== Attribution: what changed the score? ===")
    j = base_score[["pdb_id", "inplace_docking_score", "inplace_ecoul", "inplace_evdw"]].merge(
        new_score[["pdb_id", "inplace_docking_score", "inplace_ecoul", "inplace_evdw",
                   "inplace_metal"]], on="pdb_id", suffixes=("_base", "_new"))
    j["d_score"] = j["inplace_docking_score_new"] - j["inplace_docking_score_base"]
    j["metal_contrib"] = j["inplace_metal"]                       # weight ~1.0
    j["ecoul_contrib"] = (j["inplace_ecoul_new"] - j["inplace_ecoul_base"]) * 0.19
    j["evdw_contrib"] = (j["inplace_evdw_new"] - j["inplace_evdw_base"]) * 0.046
    print(j[["pdb_id", "d_score", "metal_contrib", "ecoul_contrib", "evdw_contrib"]]
          .to_string(index=False, float_format="%.2f"))
    print(f"\n  mean total change {j.d_score.mean():+.2f} | metal {j.metal_contrib.mean():+.2f} "
          f"| ecoul {j.ecoul_contrib.mean():+.2f} | evdw {j.evdw_contrib.mean():+.2f}")

    # ── metal term vs Fe distance, per pose ─────────────────────────────────
    print("\n=== Metal term vs Fe distance (all redocked poses) ===")
    poses = pd.read_csv(os.path.join(RESULTS, "redock_poses.csv"))
    dists = []
    for pdb, grp in poses.groupby("pdb_id"):
        pv = os.path.join(STRUCT, pdb, f"{pdb}_fe2_nh_redock_pv.mae")
        rec = os.path.join(STRUCT, pdb, f"{pdb}_fe2_nh_rec.mae")
        if not (os.path.exists(pv) and os.path.exists(rec)):
            continue
        fe = [a for a in structure.StructureReader.read(rec).atom if a.element == "Fe"]
        if not fe:
            continue
        fexyz = np.array(fe[0].xyz)
        for i, st in enumerate(structure.StructureReader(pv)):
            if i == 0:
                continue
            d = min(np.linalg.norm(np.array(a.xyz) - fexyz) for a in st.atom if a.element != "H")
            dists.append({"pdb_id": pdb, "rank": i, "fe_dist": d})
    dd = poses.merge(pd.DataFrame(dists), on=["pdb_id", "rank"], how="left")
    dd.to_csv(os.path.join(RESULTS, "metal_vs_distance.csv"), index=False, float_format="%.3f")

    v = dd.dropna(subset=["fe_dist", "metal"])
    if len(v):
        print(f"  n = {len(v)} poses; metal non-zero in {int((v.metal < -0.01).sum())}")
        print(f"  metal vs Fe distance: rho = {stats.spearmanr(v.fe_dist, v.metal).statistic:+.3f} "
              f"(p = {stats.spearmanr(v.fe_dist, v.metal).pvalue:.2g})")
        bins = pd.cut(v.fe_dist, [0, 2.2, 2.6, 3.0, 3.5, 4.5, 99])
        print(v.groupby(bins, observed=True).agg(n=("metal", "size"),
                                                 mean_metal=("metal", "mean"),
                                                 mean_rmsd=("rmsd", "mean")).to_string())
        # The key question: do near-native poses get MORE metal reward than displaced ones?
        near, far = v[v.rmsd <= 2], v[v.rmsd > 2]
        if len(near) and len(far):
            print(f"\n  near-native poses (RMSD<=2A, n={len(near)}): mean metal {near.metal.mean():+.3f}")
            print(f"  displaced poses  (RMSD> 2A, n={len(far)}): mean metal {far.metal.mean():+.3f}")
            print(f"  Mann-Whitney p = {stats.mannwhitneyu(near.metal, far.metal).pvalue:.3f}")

    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "comparison.csv"), index=False,
                              float_format="%.4f")
    print("\nWrote results/comparison.csv, results/metal_vs_distance.csv")


if __name__ == "__main__":
    main()
