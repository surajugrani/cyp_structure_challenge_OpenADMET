"""
Step 5: Correlate Glide scores of crystal poses with measured pIC50.

Joins data/curated_set.csv (pIC50) with results/glide_scores.csv and reports, for both
scoring modes and every Glide component term:

  - Spearman rho (primary: monotonic, robust to the non-linear score/affinity relation)
  - Pearson r
  - bootstrap 95 % CIs for both

With only ~20 points the CIs are wide, so they are reported alongside every rho: a point
estimate alone would overstate what this dataset can support.

Sensitivity checks re-run the primary correlation while dropping, in turn:
  - flagged ligands (pIC50 spread > 1.5 log units)
  - the multi-copy structures (2V0M, 4K9T)
  - ligands not coordinating the heme Fe (> 3 A)

Usage (cheap, fine on the login node):
    $SCHRODINGER/run python3 scripts/05_correlate.py

Outputs:
    results/correlations.csv     rho/r + CIs for every score column
    results/merged_data.csv      per-ligand pIC50 + scores (the analysis table)
    results/correlation_plot.png scatter plots for both scoring modes
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(HERE, "data")
RESULTS_DIR = os.path.join(HERE, "results")

N_BOOT = 1000
RNG = np.random.default_rng(0)
FE_COORD_CUTOFF = 3.0  # A; beyond this the ligand is not coordinating the iron

# Glide terms to analyse, examined in stages rather than all at once: looking at every
# component term over ~20 points invites reading noise as signal. Start with the total
# score and the van der Waals term (shape complementarity), then metal and ecoul.
# "metal" is deliberately excluded: with the chemically correct (dianionic) porphyrin it
# is identically 0.00 for every structure, so there is nothing to correlate. See README.
DEFAULT_TERMS = ["docking_score", "evdw", "ecoul"]
MODES = ["inplace", "mininplace"]


def boot_ci(x, y, func, n_boot=N_BOOT, alpha=0.05):
    """Percentile bootstrap CI for a correlation statistic."""
    n = len(x)
    if n < 4:
        return (np.nan, np.nan)
    stats_ = []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        if len(np.unique(x[idx])) < 3 or len(np.unique(y[idx])) < 3:
            continue
        stats_.append(func(x[idx], y[idx]).statistic)
    if not stats_:
        return (np.nan, np.nan)
    return tuple(np.percentile(stats_, [100 * alpha / 2, 100 * (1 - alpha / 2)]))


def correlate(df, score_col, y_col="pIC50_median"):
    sub = df[[score_col, y_col]].dropna()
    if len(sub) < 4:
        return None
    x, y = sub[score_col].to_numpy(), sub[y_col].to_numpy()
    sp = stats.spearmanr(x, y)
    pe = stats.pearsonr(x, y)
    sp_lo, sp_hi = boot_ci(x, y, stats.spearmanr)
    pe_lo, pe_hi = boot_ci(x, y, lambda a, b: stats.pearsonr(a, b))
    return {
        "score": score_col, "n": len(sub),
        "spearman_rho": sp.statistic, "spearman_p": sp.pvalue,
        "spearman_lo": sp_lo, "spearman_hi": sp_hi,
        "pearson_r": pe.statistic, "pearson_p": pe.pvalue,
        "pearson_lo": pe_lo, "pearson_hi": pe_hi,
    }


def make_plot(df, out_path):
    modes = [("inplace_docking_score", "Glide inplace"),
             ("mininplace_docking_score", "Glide mininplace")]
    fig, axes = plt.subplots(1, len(modes), figsize=(11, 5), sharey=True)
    for ax, (col, label) in zip(np.atleast_1d(axes), modes):
        sub = df.dropna(subset=[col, "pIC50_median"])
        if sub.empty:
            continue
        flagged = sub["spread_flag"].fillna(False).astype(bool)
        multi = sub["pdb_id"].isin(["2V0M", "4K9T"])
        ax.scatter(sub.loc[~flagged & ~multi, col], sub.loc[~flagged & ~multi, "pIC50_median"],
                   s=45, label="single value / narrow spread")
        ax.scatter(sub.loc[flagged, col], sub.loc[flagged, "pIC50_median"],
                   s=45, marker="s", facecolors="none", edgecolors="tab:red",
                   label="pIC50 spread > 1.5 log")
        ax.scatter(sub.loc[multi, col], sub.loc[multi, "pIC50_median"],
                   s=70, marker="^", facecolors="none", edgecolors="tab:green",
                   label="2 ligand copies")
        for _, r in sub.iterrows():
            ax.annotate(r["comp_id"], (r[col], r["pIC50_median"]),
                        fontsize=6, alpha=0.6, xytext=(3, 3), textcoords="offset points")

        rho = stats.spearmanr(sub[col], sub["pIC50_median"])
        ax.set_title(f"{label}\nSpearman rho = {rho.statistic:.2f} (p = {rho.pvalue:.3f}, n = {len(sub)})")
        ax.set_xlabel("Glide docking score (kcal/mol)")
        ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("pIC50 (median, BindingDB)")
    np.atleast_1d(axes)[0].legend(fontsize=7, loc="best")
    fig.suptitle("CYP3A4: Glide score of the crystal pose vs measured IC50", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--terms", default=",".join(DEFAULT_TERMS),
                    help=f"comma-separated Glide terms (default: {','.join(DEFAULT_TERMS)}; "
                         "next stage: metal,ecoul)")
    args = ap.parse_args()
    terms = args.terms.split(",")

    cur = pd.read_csv(os.path.join(DATA_DIR, "curated_set.csv"))
    scores = pd.read_csv(os.path.join(RESULTS_DIR, "glide_scores.csv"))
    prep = pd.read_csv(os.path.join(RESULTS_DIR, "prep_summary.csv"))

    df = cur.merge(scores, on=["pdb_id", "comp_id"], how="left")
    df = df.merge(prep[["pdb_id", "score_min_dist_to_fe"]], on="pdb_id", how="left")
    main_set = df[df["in_main_set"]].copy()
    main_set.to_csv(os.path.join(RESULTS_DIR, "merged_data.csv"), index=False,
                    float_format="%.4f")

    score_cols = [f"{mode}_{term}" for term in terms for mode in MODES
                  if f"{mode}_{term}" in df.columns]
    missing = [f"{m}_{t}" for t in terms for m in MODES if f"{m}_{t}" not in df.columns]
    if missing:
        print(f"WARNING: not in glide_scores.csv: {missing}")
    rows = [r for c in score_cols if (r := correlate(main_set, c))]
    corr = pd.DataFrame(rows).sort_values("spearman_rho")
    suffix = "" if sorted(terms) == sorted(DEFAULT_TERMS) else "_" + "-".join(terms)
    corr.to_csv(os.path.join(RESULTS_DIR, f"correlations{suffix}.csv"), index=False,
                float_format="%.4f")

    print(f"\n=== Main set: n = {main_set['inplace_docking_score'].notna().sum()} ===")
    print(corr.to_string(index=False, float_format="%.3f"))

    # ── Sensitivity checks on the primary score ─────────────────────────────
    print("\n=== Sensitivity (inplace_docking_score) ===")
    checks = {
        "all (main set)": main_set,
        "drop spread-flagged": main_set[~main_set["spread_flag"].fillna(False)],
        "drop multi-copy (2V0M, 4K9T)": main_set[~main_set["pdb_id"].isin(["2V0M", "4K9T"])],
        f"Fe-coordinating only (<{FE_COORD_CUTOFF} A)":
            main_set[main_set["score_min_dist_to_fe"] < FE_COORD_CUTOFF],
    }
    for label, sub in checks.items():
        r = correlate(sub, "inplace_docking_score")
        if r:
            print(f"  {label:<38} n={r['n']:2d}  rho={r['spearman_rho']:+.3f} "
                  f"[{r['spearman_lo']:+.2f}, {r['spearman_hi']:+.2f}]  p={r['spearman_p']:.3f}")

    make_plot(main_set, os.path.join(RESULTS_DIR, "correlation_plot.png"))


if __name__ == "__main__":
    main()
