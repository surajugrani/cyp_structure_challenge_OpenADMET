"""
Analysis 09, step 1: how well do Glide and ICM's OWN scoring functions identify the
native pose, on the 9 most rigid ligands?

This is the baseline for the pose-selector work. The competition task is pose prediction,
so the metric that matters is not affinity correlation (analyses 01/08: no signal) but
top-1 success: is the pose the scoring function ranks first within 2 A of the crystal pose?

Subset: the 9 lowest-rotatable-bond ligands of the 22 affinity-curated structures. RotB
runs 3,4,4,5,5,7,8,10,10 - there is a natural gap between 8 and 10, so "9 most rigid"
needs no tie-breaking. Rigid ligands are used because analysis 07 showed docking success
falls off steeply with flexibility, and a selector has no chance where sampling never
produces a near-native pose at all.

Poses are reused as-is:
  - Glide: analysis 05 (docked from SMILES, the realistic blind-set protocol)
  - ICM:   analysis 06 (ICM's own convertObject prep, confs=10)

Three numbers per program:
  top-1      success of the deployed scoring function  <- the headline
  best-of-N  success if an oracle picked the best pose <- the ceiling a selector could reach
  best-of-k  same with k = min(n_Glide, n_ICM), since the two return different pose counts

Usage (light - pandas only, no licence, no docking):
    $SCHRODINGER/run python3 scripts/01_baseline_success.py

Outputs:
    results/rigid9_baseline.csv   per-structure RMSDs and score-ranks
    results/rigid9_summary.csv    success rates
"""

import os
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
RESULTS = os.path.join(HERE, "results")

# comp_id, pdb_id, rotatable bonds (RDKit, from the deposited SMILES)
RIGID9 = [("MYT", "6MA6", 3), ("CL6", "8SPD", 4), ("YNV", "7LXL", 4),
          ("08Y", "3UA1", 5), ("TPF", "6MA7", 5), ("KLN", "2V0M", 7),
          ("PK9", "4D6Z", 8), ("J9K", "4D78", 10), ("PKT", "4D7D", 10)]

CUTOFFS = (2.0, 2.5)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    r9 = pd.DataFrame(RIGID9, columns=["comp_id", "pdb_id", "RotB"])

    gp = pd.read_csv(os.path.join(ANALYSES, "05_full_pdb_redocking", "results",
                                  "redock_poses.csv"))
    ip = pd.read_csv(os.path.join(ANALYSES, "06_icm_redocking", "results",
                                  "icm_poses.csv"))

    rows = []
    for r in r9.itertuples():
        # Glide poses are emitted in score order; ICM's "pose" index is NOT score order,
        # so ICM must be sorted by Score explicitly (lower = better in both programs).
        g = gp[gp.pdb_id == r.pdb_id].sort_values("rank").reset_index(drop=True)
        i = ip[ip.pdb_id == r.pdb_id].sort_values("Score").reset_index(drop=True)
        k = min(len(g), len(i))
        rows.append({
            "comp_id": r.comp_id, "pdb_id": r.pdb_id, "RotB": r.RotB, "k": k,
            "G_n": len(g), "G_top1": g.rmsd.iloc[0], "G_best": g.rmsd.min(),
            "G_rank_of_best": int(g.rmsd.idxmin()) + 1,
            "G_n_under2": int((g.rmsd < 2).sum()), "G_best_k": g.rmsd.head(k).min(),
            "I_n": len(i), "I_top1": i.rmsd.iloc[0], "I_best": i.rmsd.min(),
            "I_rank_of_best": int(i.rmsd.idxmin()) + 1,
            "I_n_under2": int((i.rmsd < 2).sum()), "I_best_k": i.rmsd.head(k).min(),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RESULTS, "rigid9_baseline.csv"), index=False,
              float_format="%.4f")
    print(df.to_string(index=False, float_format="%.2f"))

    srows = []
    for prog, pre in (("Glide", "G"), ("ICM", "I")):
        for cut in CUTOFFS:
            srows.append({
                "program": prog, "cutoff_A": cut, "n": len(df),
                "top1": int((df[f"{pre}_top1"] < cut).sum()),
                "best_of_N": int((df[f"{pre}_best"] < cut).sum()),
                "best_of_k": int((df[f"{pre}_best_k"] < cut).sum()),
            })
    s = pd.DataFrame(srows)
    s["top1_pct"] = 100 * s.top1 / s.n
    s["best_of_N_pct"] = 100 * s.best_of_N / s.n
    s.to_csv(os.path.join(RESULTS, "rigid9_summary.csv"), index=False,
             float_format="%.1f")
    print("\n" + s.to_string(index=False, float_format="%.1f"))
    print(f"\nWrote 2 CSVs to {RESULTS}")


if __name__ == "__main__":
    main()
