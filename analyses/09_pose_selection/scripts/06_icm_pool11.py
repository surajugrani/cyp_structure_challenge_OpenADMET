"""
Analysis 09, step 6: build the 11-pose ICM pool (docked + `-n` near-native) and fit the
pooled linear selector on it, measuring how much of any gain is the `-n` route artefact.

Step 4 recorded the problem: `-n` scores the same geometry 0.25-13.4 Score units better
than `_dockScan`, in 7 cases out of 7. Since every ligand's `-n` row is also its only
guaranteed near-native row, a model that simply learns "pick the row with the suspiciously
good score" scores 30/30 and is worthless blind.

So three numbers are reported, not one:

  1. top-1 < 2 A on the full 11-pose pool          - the headline, but inflatable
  2. how often the model picks the `-n` row        - the artefact meter. 30/30 means the
                                                     model learned the route, not the pose
  3. top-1 < 2 A with the `-n` row REMOVED AT TEST TIME but kept in training - the honest
     number, directly comparable to step 3's 8/30 baseline, because the candidate set is
     then the same docked poses a blind run would produce

Evaluation is leave-one-structure-out throughout.

Usage (light - pandas/numpy only, no licence, no docking):
    $SCHRODINGER/run python3 scripts/06_icm_pool11.py

Outputs:
    results/icm_pool11_poses.csv
    results/icm_pool11_loso.csv
"""

import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A06 = os.path.join(ANALYSES, "06_icm_redocking")
RESULTS = os.path.join(HERE, "results")

FEATS = ["Score", "RTCNN", "mfScore_sdf"]
CUT = 2.0


def sdf_tags(path):
    out = []
    for b in open(path).read().split("$$$$"):
        if "<Score>" not in b:
            continue
        out.append({m.group(1): m.group(2).strip()
                    for m in re.finditer(r"> *<([^>]+)>\s*\n([^\n]*)", b)})
    return out


def fit_ols(X, y):
    A = np.column_stack([X, np.ones(len(X))])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c[:-1], c[-1]


def build_pool():
    docked = pd.read_csv(os.path.join(RESULTS, "rtcnn_rescore_poses.csv"))
    docked = docked[["pdb_id", "comp_id", "rmsd"] + FEATS].copy()
    docked["is_native_row"] = False

    nr = pd.read_csv(os.path.join(RESULTS, "native_pose_rmsd.csv")).set_index("pdb_id")
    rows = []
    for pdb in sorted(docked.pdb_id.unique()):
        f = os.path.join(A06, "structures", pdb, f"{pdb}_native_hits.sdf")
        t = sdf_tags(f)[0]
        rows.append({"pdb_id": pdb,
                     "comp_id": docked[docked.pdb_id == pdb].comp_id.iloc[0],
                     "rmsd": float(nr.loc[pdb, "n_rmsd"]),
                     "Score": float(t["Score"]), "RTCNN": float(t["RTCNN"]),
                     "mfScore_sdf": float(t["mfScore"]), "is_native_row": True})
    pool = pd.concat([docked, pd.DataFrame(rows)], ignore_index=True)
    pool.to_csv(os.path.join(RESULTS, "icm_pool11_poses.csv"),
                index=False, float_format="%.4f")
    return pool


def main():
    p = build_pool()
    pdbs = sorted(p.pdb_id.unique())
    n = len(pdbs)
    print(f"pool: {len(p)} rows over {n} structures "
          f"({int((~p.is_native_row).sum())} docked + {int(p.is_native_row.sum())} -n)")
    print(f"oracle on full pool: {sum(p[p.pdb_id == q].rmsd.min() < CUT for q in pdbs)}/{n}")

    # how separable is the -n row by score alone, per ligand?
    ranks = []
    for q in pdbs:
        s = p[p.pdb_id == q].sort_values("Score").reset_index(drop=True)
        ranks.append(int(s.index[s.is_native_row][0]) + 1)
    ranks = np.array(ranks)
    print(f"\n-n row's rank by raw Score within its ligand: "
          f"1st in {(ranks == 1).sum()}/{n}, median rank {np.median(ranks):.0f}")

    rows = []
    for held in pdbs:
        tr, te = p[p.pdb_id != held], p[p.pdb_id == held]
        w, b = fit_ols(tr[FEATS].values, tr.rmsd.values)
        pred = te[FEATS].values @ w + b
        pick = te.iloc[int(np.argmin(pred))]

        # same model, but the -n row is not a candidate at test time
        td = te[~te.is_native_row]
        pd_ = td[FEATS].values @ w + b
        pick_d = td.iloc[int(np.argmin(pd_))]

        rows.append({"pdb_id": held, "comp_id": te.comp_id.iloc[0],
                     "n_rows": len(te),
                     "picked_native_row": bool(pick.is_native_row),
                     "top1_full_pool": pick.rmsd,
                     "top1_docked_only": pick_d.rmsd,
                     "score_baseline": td.loc[td.Score.idxmin(), "rmsd"],
                     "oracle_docked": td.rmsd.min()})
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(RESULTS, "icm_pool11_loso.csv"),
             index=False, float_format="%.4f")

    print(f"\n=== leave-one-structure-out, n = {n} ===")
    print(f"  1. top-1 < {CUT} A on the full 11-row pool   "
          f"{int((d.top1_full_pool < CUT).sum())}/{n}")
    print(f"  2. picked the -n row                        "
          f"{int(d.picked_native_row.sum())}/{n}   <- artefact meter")
    print(f"  3. top-1 < {CUT} A, -n row removed at test   "
          f"{int((d.top1_docked_only < CUT).sum())}/{n}   <- honest number")
    print(f"     ICM Score baseline (docked only)         "
          f"{int((d.score_baseline < CUT).sum())}/{n}")
    print(f"     oracle (docked only)                     "
          f"{int((d.oracle_docked < CUT).sum())}/{n}")

    # of the times it picked the -n row, how often was that the right call?
    pn = d[d.picked_native_row]
    if len(pn):
        rescued = pn[pn.oracle_docked >= CUT]
        print(f"\n  picked -n in {len(pn)}: of those, {len(rescued)} were ligands with NO "
              f"docked pose < {CUT} A\n  (i.e. cases where only the -n row could succeed)")


if __name__ == "__main__":
    main()
