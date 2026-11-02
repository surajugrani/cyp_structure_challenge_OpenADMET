"""Submission 38: the same 200 ICM redocked poses, picked by ICM's OWN score.

The companion attempt - ICM redock of submission 29 with the Glide-picked pose and the
PrepWizard-relaxed receptor - scored 0.3293 / 5.0863 / 0.7959 against submission 29's
0.4638 / 3.7277 / 0.8001. LDDT-LP barely moved, so the relaxed receptor was not the
cost; the ligand placements were.

This changes ONLY the selector: lowest ICM `Score` per ligand instead of lowest Glide
`inplace` gscore, over the identical pool of 200 poses from the ICM redock of submission
29's receptors. It changes 18 of 20 structures against the Glide-picked version - the two
selectors agree almost nowhere, and neither ever picks ICM's own rank-1 pose.

EXPECT THIS TO LOSE. The pool it draws from is the one that just produced the project's
worst result by a wide margin. A different ranking of poses that are collectively 5.09 A
from the target has no obvious route back to 0.4638. The question it answers is narrow
and worth one attempt: was the 0.3293 the fault of Glide scoring poses against a receptor
it did not see, or of the ICM poses themselves? If ICM's own selector lands near 0.3293
too, the poses are the problem and this whole redocking route closes.

The receptor shipped is the PrepWizard-relaxed one from the per-pose rebuild, same as the
Glide-picked version, because the raw ICM poses fail PoseBusters 14/20 on heavy-atom
clashes against submission 29's untouched receptor.

ONE FILTER, NOT A RANKING. A first pass ranked purely on ICM `Score` and failed
PoseBusters 17/20: OCNT-2395415, OCNT-2395460 and OCNT-2395503 all picked pose 10, and
all three had been rejected by Glide at the 10000.0 ceiling. The Glide-picked version
passed 20/20 for free, because a clashing pose scores 10000.0 and drops out on its own;
ICM's score knows nothing about the Glide receptor and will happily pick one.

So poses without a usable Glide score are excluded first, then ICM `Score` ranks what is
left. Glide decides only whether a pose is physically admissible, never which is best -
that stays ICM's call, which is the whole point of this attempt. It changes exactly those
3 ligands against the pure-ICM ranking.
"""
import collections
import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A33 = os.path.join(ROOT, "analyses", "33_icm_redock_sub29", "results", "icm_poses.csv")
A36 = os.path.join(ROOT, "analyses", "36_prepped_icm_poses", "work")
GLIDE_PICK = os.path.join(ROOT, "submissions", "37_icm_redock_glide_select",
                          "results", "selection.csv")
SHARDS = os.path.join(ROOT, "analyses", "36_prepped_icm_poses", "results",
                      "scores_shard*.csv")
REJECT = 100.0          # Glide's rejection ceiling is 10000.0
OUT = os.path.join(HERE, "results", "selection.csv")


def main():
    poses = collections.defaultdict(dict)
    for r in csv.DictReader(open(A33)):
        poses[r["ligand"]][int(r["pose_rank"])] = r
    if len(poses) != 20:
        raise SystemExit(f"expected 20 ligands, found {len(poses)}")
    prev = {r["ligand"]: int(r["pose_rank"]) for r in csv.DictReader(open(GLIDE_PICK))}

    admissible = set()
    for f in glob.glob(SHARDS):
        for r in csv.DictReader(open(f)):
            if r["status"] == "ok" and r.get("gscore") and float(r["gscore"]) < REJECT:
                admissible.add((r["ligand"], int(r["pose_rank"])))

    rows = []
    for lig in sorted(poses):
        ok = [k for k in poses[lig] if (lig, k) in admissible]
        if not ok:
            raise SystemExit(f"{lig}: no pose with a usable Glide score")
        rank = min(ok, key=lambda k: float(poses[lig][k]["Score"]))
        raw = min(poses[lig], key=lambda k: float(poses[lig][k]["Score"]))
        mae = os.path.join(A36, lig, f"pose{rank:02d}", "complex.mae")
        if not os.path.exists(mae):
            raise SystemExit(f"{lig}: no relaxed complex for pose {rank} ({mae})")
        src = poses[lig][rank]
        rows.append({"ligand": lig, "pose_rank": rank,
                     "icm_score": src["Score"], "rtcnn": src["RTCNN"],
                     "rmsd_to_sub29_ligand": src["rmsd_to_cofolded"],
                     "fe_dist": src["fe_dist"],
                     "n_admissible_poses": len(ok),
                     "pure_icm_rank": raw,
                     "filtered_by_glide": "yes" if raw != rank else "no",
                     "glide_pick_rank": prev.get(lig, ""),
                     "changed_vs_glide_pick":
                         "yes" if rank != prev.get(lig) else "no"})

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    d = sorted(float(r["rmsd_to_sub29_ligand"]) for r in rows)
    chg = sum(1 for r in rows if r["changed_vs_glide_pick"] == "yes")
    nf = sum(1 for r in rows if r["filtered_by_glide"] == "yes")
    print(f"selected 20/20 by ICM Score among Glide-admissible poses")
    print(f"differs from the pure ICM ranking on {nf}/20 (clashing poses excluded)")
    print(f"changed vs the Glide-picked version: {chg}/20")
    print(f"ligand RMSD to submission 29: median {(d[9]+d[10])/2:.2f} A, "
          f"min {d[0]:.2f}, max {d[-1]:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
