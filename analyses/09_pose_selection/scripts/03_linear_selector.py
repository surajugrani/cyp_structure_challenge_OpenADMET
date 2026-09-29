"""
Analysis 09, step 3: fit a linear combination of ICM's Score, RTCNN and mfScore to
predict pose RMSD, and test whether it selects poses better than ICM's Score alone.

Rationale. Step 2a showed the three ICM scores succeed on DIFFERENT structures - Score
8/30, RTCNN 6/30, mfScore 3/30 top-1 under 2 A, with mfScore uniquely rescuing 7LXL and
9BV6 and RTCNN uniquely rescuing 6MA8. A weighted combination should be able to use all
three. This is the pooled version: one set of weights, fitted across all poses of all
ligands at once.

Two design decisions:

  1. TRAINING SET. In 10 of the 30 structures no pose is within 2 A, so all their label
     variation is among poses that are all wrong - noise for this purpose, and ~1/3 of the
     training poses pulling the fit toward the bad-pose regime. Both training sets are
     tried: all 30, and only the 20 where a near-native pose exists.

     The TEST set is always all 30. On the blind set you do not know in advance which
     ligands have a good pose available, so restricting evaluation would inflate the
     success rate (oracle would be 20/20 by construction).

  2. SPLIT BY STRUCTURE, not by pose. Poses of one ligand are highly correlated, so a
     pose-level split leaks badly. Effective n is 30 structures, not ~470 poses.
     Evaluation is leave-one-structure-out: fit on the other structures, rank the held-out
     structure's poses, take the top-1.

Reported metric is top-1 success under 2 A, not R^2 - a model can improve R^2 while
selecting worse poses, because most pooled RMSD variance is between ligands rather than
between poses of one ligand.

Usage (light - pandas/numpy only, no licence, no docking):
    $SCHRODINGER/run python3 scripts/03_linear_selector.py

Inputs:  results/rtcnn_rescore_poses.csv  (written by 02_rtcnn_rescore.py)
Outputs: results/linear_selector_weights.csv
         results/linear_selector_loso.csv
"""

import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")

FEATS = ["Score", "RTCNN", "mfScore_sdf"]
LABELS = {"Score": "Score", "RTCNN": "RTCNN", "mfScore_sdf": "mfScore"}
CUT = 2.0


def fit_ols(X, y):
    """Least squares with intercept. Returns (weights, intercept)."""
    A = np.column_stack([X, np.ones(len(X))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef[:-1], coef[-1]


def top1_rmsd(sub, pred):
    """Actual RMSD of the pose the model ranks first (lowest predicted RMSD)."""
    return sub.rmsd.values[int(np.argmin(pred))]


def main():
    p = pd.read_csv(os.path.join(RESULTS, "rtcnn_rescore_poses.csv"))
    pdbs = sorted(p.pdb_id.unique())
    n = len(pdbs)

    solvable = {q for q in pdbs if (p[p.pdb_id == q].rmsd < CUT).any()}
    print(f"{len(p)} poses over {n} structures; "
          f"{len(solvable)} have at least one pose < {CUT} A")
    print(f"poses < {CUT} A: {(p.rmsd < CUT).sum()}/{len(p)} "
          f"({100 * (p.rmsd < CUT).mean():.1f}%)")

    TRAIN_SETS = {"train_all30": set(pdbs), "train_solvable20": solvable}

    # ── full-data fits, for the weights themselves ──────────────────────────
    wrows = []
    for name, keep in TRAIN_SETS.items():
        tr = p[p.pdb_id.isin(keep)]
        w, b = fit_ols(tr[FEATS].values, tr.rmsd.values)
        r = np.corrcoef(tr[FEATS].values @ w + b, tr.rmsd.values)[0, 1]
        print(f"\n[{name}] {len(tr)} poses. "
              f"Weights (positive = higher value means WORSE pose):")
        for f, wi in zip(FEATS, w):
            print(f"    {LABELS[f]:9} {wi:+.4f}")
        print(f"    {'intercept':9} {b:+.4f}   in-sample r = {r:+.3f}")
        wrows.append({"train_set": name, "n_poses": len(tr),
                      **{LABELS[f]: wi for f, wi in zip(FEATS, w)},
                      "intercept": b, "in_sample_r": r})
    pd.DataFrame(wrows).to_csv(os.path.join(RESULTS, "linear_selector_weights.csv"),
                               index=False, float_format="%.4f")

    # ── leave-one-structure-out, always tested on all 30 ────────────────────
    rows = []
    for held in pdbs:
        te = p[p.pdb_id == held]
        rec = {"pdb_id": held, "comp_id": te.comp_id.iloc[0], "n_poses": len(te),
               "solvable": held in solvable, "oracle": te.rmsd.min()}
        for f in FEATS:
            rec[LABELS[f] + "_top1"] = top1_rmsd(te, te[f].values)
        for name, keep in TRAIN_SETS.items():
            tr = p[p.pdb_id.isin(keep - {held})]
            w, b = fit_ols(tr[FEATS].values, tr.rmsd.values)
            rec[name + "_top1"] = top1_rmsd(te, te[FEATS].values @ w + b)
        rows.append(rec)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(RESULTS, "linear_selector_loso.csv"),
             index=False, float_format="%.4f")

    print(f"\n=== leave-one-structure-out, tested on all {n}, top-1 < {CUT} A ===")
    for c in [LABELS[f] + "_top1" for f in FEATS] + \
             [t + "_top1" for t in TRAIN_SETS] + ["oracle"]:
        k = int((d[c] < CUT).sum())
        print(f"  {c:22} {k:2d}/{n} ({100 * k / n:4.1f}%)")
    exp = np.mean([(p[p.pdb_id == q].rmsd < CUT).mean() for q in pdbs])
    print(f"  {'random pose':22} {exp * n:4.1f}/{n} ({100 * exp:4.1f}%)  [expectation]")

    for name in TRAIN_SETS:
        c = name + "_top1"
        sv = d[(d.Score_top1 >= CUT) & (d[c] < CUT)]
        br = d[(d.Score_top1 < CUT) & (d[c] >= CUT)]
        print(f"\n{name} vs Score: saves {len(sv)}, breaks {len(br)}")
        for lab, sub in (("SAVE", sv), ("BREAK", br)):
            for r in sub.itertuples():
                print(f"  {lab:5} {r.pdb_id} ({r.comp_id}) "
                      f"Score {r.Score_top1:.2f} -> {getattr(r, c):.2f} "
                      f"(oracle {r.oracle:.2f})")
    print(f"\nWrote 2 CSVs to {RESULTS}")


if __name__ == "__main__":
    main()
