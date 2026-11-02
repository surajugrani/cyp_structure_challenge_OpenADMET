"""Pick the best-scoring prepped ICM pose per ligand, then measure how far it sits from
submission 29's ligand.

Analysis 34 scored the raw ICM poses against the co-folded receptor and 177 of 200 came
back at the gscore 10000 rejection ceiling, because the receptor's hydrogen network had
been optimised around the co-folded ligand and overlapped the ICM pose at 0.73-1.6 A.
Analysis 36 rebuilt each complex, re-ran PrepWizard with `-rehtreat` at `-rmsd 0.6` so
the H-network relaxes around the ICM pose, regenerated the grid, and scored `inplace`.

This step does the selection: lowest gscore per ligand across its 10 poses.

The RMSD column is a DISTANCE, not a verdict. Submission 29 scored 0.4638, the best this
project has managed, so the separation says how far a submission built from these poses
would move from the one that works - nothing about which of the two is closer to the
truth. Only the leaderboard answers that.

rmsd_to_cofolded is taken from analysis 33, which computed it symmetry-corrected with
spyrmsd against submission 29's own ligand, in the same coordinate frame (the ICM
docking used submission 29's receptors), so no superposition is involved.
"""
import collections
import csv
import glob
import os
import statistics

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A33 = os.path.join(ROOT, "analyses", "33_icm_redock_sub29", "results", "icm_poses.csv")
OUT = os.path.join(HERE, "results", "best_pose_per_ligand.csv")
REJECT = 100.0          # Glide's rejection ceiling is 10000.0; nothing real is near it


def main():
    scored = []
    for f in sorted(glob.glob(os.path.join(HERE, "results", "scores_shard*.csv"))):
        scored += list(csv.DictReader(open(f)))
    print(f"poses scored: {len(scored)}")

    a33 = {(r["ligand"], int(r["pose_rank"])): r for r in csv.DictReader(open(A33))}

    by = collections.defaultdict(list)
    for r in scored:
        if r["status"] == "ok" and r.get("gscore") and float(r["gscore"]) < REJECT:
            by[r["ligand"]].append(r)
    n_rej = sum(1 for r in scored if r["status"] != "ok"
                or not r.get("gscore") or float(r["gscore"]) >= REJECT)
    print(f"usable: {len(scored) - n_rej}/{len(scored)}   "
          f"ligands with >=1 usable pose: {len(by)}/20")

    rows = []
    for lig in sorted(by):
        best = min(by[lig], key=lambda r: float(r["gscore"]))
        rank = int(best["pose_rank"])
        src = a33.get((lig, rank), {})
        rows.append({
            "ligand": lig,
            "picked_set": best["picked_set"],
            "pose_rank": rank,
            "n_usable_poses": len(by[lig]),
            "gscore": round(float(best["gscore"]), 3),
            "lig_shift_during_prep": best.get("lig_shift_during_prep", ""),
            "rmsd_to_sub29_ligand": src.get("rmsd_to_cofolded", ""),
            "icm_score": src.get("Score", ""),
            "fe_dist": src.get("fe_dist", ""),
            "rmsd_best_available": min(
                (float(a33[(lig, k)]["rmsd_to_cofolded"])
                 for k in range(1, 11) if (lig, k) in a33), default=""),
        })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    d = sorted(float(r["rmsd_to_sub29_ligand"]) for r in rows
               if r["rmsd_to_sub29_ligand"])
    print(f"\nwrote {os.path.relpath(OUT, ROOT)}  ({len(rows)} ligands)")
    print(f"RMSD of picked pose to submission 29's ligand:")
    print(f"  median {statistics.median(d):.2f} A   min {d[0]:.2f}   max {d[-1]:.2f}")
    print(f"  within 1 A: {sum(1 for x in d if x < 1):2d}/{len(d)}")
    print(f"  within 2 A: {sum(1 for x in d if x < 2):2d}/{len(d)}")
    print(f"  beyond 4 A: {sum(1 for x in d if x >= 4):2d}/{len(d)}")
    print(f"\n{'ligand':<16}{'pose':>5}{'n':>4}{'gscore':>9}{'rmsd_to_29':>12}"
          f"{'best_avail':>12}")
    for r in rows:
        print(f"{r['ligand']:<16}{r['pose_rank']:>5}{r['n_usable_poses']:>4}"
              f"{r['gscore']:>9.2f}{float(r['rmsd_to_sub29_ligand']):>12.2f}"
              f"{float(r['rmsd_best_available']):>12.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
