"""
Analysis 09, step 7: the within-ligand linear selector.

Task: learn to identify native poses from ICM's score terms. Features are the score terms,
target is the pose's RMSD to the crystal pose.

Why within-ligand. Step 3 and step 6 fitted the same features pooled across ligands and
both came out WORSE than ICM's Score alone, with a negative weight on Score - the model
concluding that a better score means a worse pose. The cause is a between-ligand confound:
across structures, mean Score vs oracle RMSD gives rho = -0.305 (ligands that score less
favourably are the smaller, more rigid ones, and they dock better), which swamps the
within-ligand signal (Score rho = +0.017, RTCNN +0.334).

The fix is a feature transform, NOT 30 separate models: z-score each ligand's rows against
each other, then fit ONE global set of weights on the standardised values. The ligand-level
offset is removed; the pose-to-pose differences that ranking actually depends on survive.

Set: 11 rows per ligand - the best 10 docked poses by Score plus the `-n` near-native row
(329 rows; 8SO2 has only 9 docked poses so contributes 10). Capping at 10 docked rows
matters: on the uncapped 499-row pool, ligands with 20-30 poses dominated the fit and
performance halved (4/30 vs 8/30).

The `-n` row is included in training and as a test candidate, so the oracle is 30/30.

Evaluation is leave-one-structure-out: fit on 29 structures, standardise and rank the
held-out structure's rows, take the top-1. Reported metric is top-1 < 2 A.

Two target variants are fitted, since it is not obvious which is right:
  raw   - predict RMSD in Angstroms. Keeps the notion of "how wrong", but a ligand whose
          poses are all 5-11 A contributes large residuals that dominate least squares.
  z     - predict within-ligand z-scored RMSD. Every ligand contributes equally and the
          model only has to learn relative ordering, which is all ranking needs.

Usage (light - pandas/numpy only, no licence, no docking):
    $SCHRODINGER/run python3 scripts/07_within_ligand.py

Outputs:
    results/within_ligand_loso.csv
    results/within_ligand_weights.csv
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")

CUT = 2.0
MAX_DOCKED = 10
FSETS = {"Score+RTCNN": ["Score", "RTCNN"],
         "Score+RTCNN+mfScore": ["Score", "RTCNN", "mfScore_sdf"]}


def zwithin(df, cols):
    """z-score each column inside each pdb_id; constant column -> 0, not NaN."""
    out = df.copy()
    for c in cols:
        g = out.groupby("pdb_id")[c]
        sd = g.transform("std")
        out[c + "_z"] = np.where(sd > 1e-9, (out[c] - g.transform("mean")) / sd, 0.0)
    return out


def fit_ols(X, y):
    A = np.column_stack([X, np.ones(len(X))])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c[:-1], c[-1]


def build_set():
    p = pd.read_csv(os.path.join(RESULTS, "icm_pool11_poses.csv"))
    parts = []
    for q in sorted(p.pdb_id.unique()):
        s = p[p.pdb_id == q]
        parts.append(pd.concat([s[~s.is_native_row].nsmallest(MAX_DOCKED, "Score"),
                                s[s.is_native_row]]))
    return pd.concat(parts, ignore_index=True)


def main():
    d = build_set()
    allf = sorted({c for F in FSETS.values() for c in F})
    d = zwithin(d, allf)
    g = d.groupby("pdb_id").rmsd
    sd = g.transform("std")
    d["rmsd_z"] = np.where(sd > 1e-9, (d.rmsd - g.transform("mean")) / sd, 0.0)

    pdbs = sorted(d.pdb_id.unique())
    n = len(pdbs)
    print(f"{len(d)} rows over {n} structures "
          f"({d.groupby('pdb_id').size().min()}-{d.groupby('pdb_id').size().max()} each)")
    print(f"oracle: {sum(d[d.pdb_id == q].rmsd.min() < CUT for q in pdbs)}/{n}")

    base = sum(d[(d.pdb_id == q) & (~d.is_native_row)]
               .sort_values("Score").iloc[0].rmsd < CUT for q in pdbs)
    base11 = sum(d[d.pdb_id == q].sort_values("Score").iloc[0].rmsd < CUT for q in pdbs)
    print(f"baseline, ICM Score over the 10 docked rows: {base}/{n}")
    print(f"baseline, ICM Score over all 11 rows:        {base11}/{n}")

    wrows, rows = [], []
    print(f"\n=== leave-one-structure-out, top-1 < {CUT} A, n = {n} ===")
    for fname, F in FSETS.items():
        Z = [c + "_z" for c in F]
        for tname, tgt in (("raw", "rmsd"), ("z", "rmsd_z")):
            w, b = fit_ols(d[Z].values, d[tgt].values)
            wrows.append({"features": fname, "target": tname,
                          **{c: wi for c, wi in zip(F, w)}, "intercept": b})
            hits = nat = 0
            per = []
            for h in pdbs:
                tr, te = d[d.pdb_id != h], d[d.pdb_id == h]
                w, b = fit_ols(tr[Z].values, tr[tgt].values)
                pred = te[Z].values @ w + b
                pick = te.iloc[int(np.argmin(pred))]
                hits += pick.rmsd < CUT
                nat += bool(pick.is_native_row)
                per.append({"pdb_id": h, "features": fname, "target": tname,
                            "top1": pick.rmsd, "picked_native": bool(pick.is_native_row),
                            "oracle": te.rmsd.min()})
            rows += per
            print(f"  {fname:22} target={tname:3}  {hits:2d}/{n}   "
                  f"(picked -n row {nat}/{n})")

    pd.DataFrame(wrows).to_csv(os.path.join(RESULTS, "within_ligand_weights.csv"),
                               index=False, float_format="%.4f")
    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "within_ligand_loso.csv"),
                              index=False, float_format="%.4f")

    print("\nfull-set weights on z-scored features "
          "(positive = higher value means WORSE pose):")
    for r in wrows:
        terms = "  ".join(f"{k}={r[k]:+.4f}" for k in r
                          if k not in ("features", "target", "intercept"))
        print(f"  {r['features']:22} target={r['target']:3}  {terms}")

    print("\nwithin-ligand Spearman of each z-scored term vs RMSD "
          "(median over structures):")
    for c in allf:
        v = [stats.spearmanr(s[c], s.rmsd).statistic
             for _, s in d.groupby("pdb_id") if s[c].nunique() > 2]
        v = [x for x in v if not np.isnan(x)]
        print(f"  {c.replace('_sdf',''):10} {np.median(v):+.3f}  "
              f"(positive in {sum(x > 0 for x in v)}/{len(v)})")
    print(f"\nWrote 2 CSVs to {RESULTS}")


if __name__ == "__main__":
    main()
