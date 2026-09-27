"""
Step 2: Symmetry-corrected RMSD of every redocked pose to the crystal pose.

Glide reports this as the per-pose property `r_i_glide_rmsd_to_input` (enabled by
CALC_INPUT_RMS), and since the ligand file IS the crystal pose that is exactly what we
want. Glide's value is used as primary.

Cross-check: Glide preserves the input atom order, so a naive index-wise heavy-atom RMSD
uses the true atom correspondence. Symmetry correction can only *lower* the RMSD, so a
valid Glide value must satisfy `glide <= naive` (within rounding). Any violation is
flagged in the summary.

(An earlier version cross-checked with schrodinger.structutils.rmsd.ConformerRmsd. That
was dropped: its SMARTS-based renumbering picked degenerate atom mappings for some
ligands and *inflated* the RMSD — e.g. 4K9T pose 1 gave 5.85 A where both the naive
index-wise value and Glide agree on 1.708 A.)

RMSD is measured **in place** (no superposition): the receptor is identical in both
cases, so the poses are already in the same frame. Superimposing first would hide exactly
the displacement we are trying to measure. Both methods are symmetry-corrected, so a
phenyl ring flipped 180 degrees is not counted as an error.

For each ligand we report:
    rmsd_top1     : RMSD of the pose Glide ranked first  -> can the score PICK the pose?
    rmsd_best     : lowest RMSD among all poses          -> can sampling FIND the pose?
    rank_of_best  : where that best pose was ranked

Usage (via scripts/run_redock.sh):
    $SCHRODINGER/run python3 scripts/02_rmsd.py

Outputs (committed): results/pose_rmsd.csv, results/pose_rmsd_all_poses.csv
"""

import os

import numpy as np
import pandas as pd
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
STRUCT_DIR = os.path.join(A01, "structures")
RESULTS_DIR = os.path.join(HERE, "results")

SUCCESS_A = 2.0   # field-standard "correct pose" cutoff
STRICT_A = 1.0


def heavy_xyz(st):
    return np.array([a.xyz for a in st.atom if a.element != "H"])


def pose_rmsds(crystal_sdf, pv_mae):
    """(glide_rmsds, naive_rmsds) for each docked pose vs the crystal ligand.

    The poseviewer file is [receptor, pose1, pose2, ...]; the first entry is skipped.
    """
    ref = structure.StructureReader.read(crystal_sdf)
    ref_xyz = heavy_xyz(ref)
    glide, naive = [], []
    for i, st in enumerate(structure.StructureReader(pv_mae)):
        if i == 0:
            continue  # receptor entry
        glide.append(st.property.get("r_i_glide_rmsd_to_input"))
        pose_xyz = heavy_xyz(st)
        if pose_xyz.shape != ref_xyz.shape:
            raise RuntimeError(f"atom count mismatch: {ref_xyz.shape} vs {pose_xyz.shape}")
        naive.append(float(np.sqrt(((ref_xyz - pose_xyz) ** 2).sum(axis=1).mean())))
    return glide, naive


def main():
    cur = pd.read_csv(os.path.join(A01, "data", "curated_set.csv"))
    os.makedirs(RESULTS_DIR, exist_ok=True)

    summary, per_pose = [], []
    for _, row in cur.iterrows():
        pdb, comp = row["pdb_id"], row["comp_id"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        crystal = os.path.join(wdir, f"{pdb}_lig.sdf")
        pv = os.path.join(wdir, f"{pdb}_redock_pv.mae")
        glide_csv = os.path.join(wdir, f"{pdb}_redock.csv")
        rec = {"pdb_id": pdb, "comp_id": comp}
        try:
            if not os.path.exists(pv):
                raise RuntimeError("no redock poseviewer file (run step 1 first)")
            rmsds, check = pose_rmsds(crystal, pv)
            if not rmsds:
                raise RuntimeError("no poses in poseviewer file")
            if any(r is None for r in rmsds):   # CALC_INPUT_RMS missing -> use our own
                rmsds = check

            scores = pd.read_csv(glide_csv)["r_i_docking_score"].tolist() \
                if os.path.exists(glide_csv) else [None] * len(rmsds)

            for rank, r in enumerate(rmsds, 1):
                per_pose.append({"pdb_id": pdb, "comp_id": comp, "rank": rank,
                                 "rmsd": r, "rmsd_naive": check[rank - 1],
                                 "docking_score": scores[rank - 1] if rank <= len(scores) else None})

            best_i = min(range(len(rmsds)), key=lambda i: rmsds[i])
            rec.update({
                "n_poses": len(rmsds),
                "rmsd_top1": rmsds[0],
                "rmsd_best": rmsds[best_i],
                "rank_of_best": best_i + 1,
                "score_top1": scores[0] if scores else None,
                # symmetry correction can only lower RMSD; glide > naive would be a bug
                "max_symmetry_gain": max(b - a for a, b in zip(rmsds, check)),
                "check_violation": max(a - b for a, b in zip(rmsds, check)) > 0.01,
                "status": "ok",
            })
            print(f"{pdb} {comp:>5}: top1={rmsds[0]:6.2f}  best={rmsds[best_i]:6.2f} "
                  f"(rank {best_i + 1}/{len(rmsds)})")
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"{pdb} {comp:>5}: ERROR {e}")
        summary.append(rec)

    df = pd.DataFrame(summary)
    df.to_csv(os.path.join(RESULTS_DIR, "pose_rmsd.csv"), index=False, float_format="%.3f")
    pd.DataFrame(per_pose).to_csv(os.path.join(RESULTS_DIR, "pose_rmsd_all_poses.csv"),
                                  index=False, float_format="%.3f")

    ok = df[df["status"] == "ok"]
    if len(ok):
        print(f"\n=== Pose recovery, n = {len(ok)} ===")
        for label, col in (("top-1", "rmsd_top1"), ("best-of-N", "rmsd_best")):
            for cut in (SUCCESS_A, STRICT_A):
                n = (ok[col] <= cut).sum()
                print(f"  {label:<10} RMSD <= {cut:.1f} A: {n:2d}/{len(ok)} ({100*n/len(ok):.0f} %)")
            print(f"  {label:<10} median RMSD: {ok[col].median():.2f} A")
        nbad = int(ok["check_violation"].sum())
        print(f"  cross-check: max symmetry gain {ok['max_symmetry_gain'].max():.3f} A; "
              f"{nbad} structure(s) with glide > naive"
              + ("  <-- INVESTIGATE" if nbad else " (OK)"))
    print("\nWrote results/pose_rmsd.csv, results/pose_rmsd_all_poses.csv")


if __name__ == "__main__":
    main()
