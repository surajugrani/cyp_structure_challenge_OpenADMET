"""Does Glide rank ICM's 10 poses better than ICM's own Score?

No ground truth exists for the blinded set, so the yardstick is the co-folded pose from
submission 06 - the structure that scored 0.4291 LDDT-PLI. "Better" here means a pose
closer to that, which is a proxy, not truth.
"""
import csv, os, statistics as st
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
A25 = os.path.join(ROOT, "analyses", "25_glide_on_icm_poses", "results")

by = defaultdict(list)
for r in csv.DictReader(open(os.path.join(A25, "glide_on_icm_poses.csv"))):
    if r["status"] != "ok" or not r.get("gscore") or not r.get("rmsd_to_cofolded"):
        continue
    try:
        by[r["ligand"]].append({
            "rank": int(r["icm_rank"]), "g": float(r["gscore"]),
            "em": float(r["emodel"]) if r.get("emodel") else None,
            "icm": float(r["icm_score"]), "rtcnn": float(r["icm_rtcnn"]),
            "rmsd": float(r["rmsd_to_cofolded"]), "fe": float(r["fe_dist"])})
    except (ValueError, TypeError):
        pass

print(f"{len(by)} ligands, {sum(len(v) for v in by.values())} poses scored\n")
print(f"{'ligand':16} {'ICM pick':>18} {'Glide pick':>18} {'oracle':>8}")
rows, agg = [], defaultdict(list)
for lg in sorted(by):
    p = by[lg]
    icm = min(p, key=lambda x: x["icm"])
    gl = min(p, key=lambda x: x["g"])
    em = min(p, key=lambda x: x["em"] if x["em"] is not None else 1e9)
    orc = min(p, key=lambda x: x["rmsd"])
    for k, v in (("icm_Score", icm), ("glide_gscore", gl), ("glide_emodel", em),
                 ("oracle", orc)):
        agg[k].append(v["rmsd"])
    agg["fe_icm"].append(icm["fe"]); agg["fe_glide"].append(gl["fe"])
    rows.append({"ligand": lg, "icm_rank": icm["rank"], "icm_rmsd": icm["rmsd"],
                 "icm_fe": icm["fe"], "glide_rank": gl["rank"],
                 "glide_rmsd": gl["rmsd"], "glide_fe": gl["fe"],
                 "oracle_rank": orc["rank"], "oracle_rmsd": orc["rmsd"]})
    print(f"{lg:16} r{icm['rank']:<2} {icm['rmsd']:6.2f}A Fe{icm['fe']:5.2f}  "
          f"r{gl['rank']:<2} {gl['rmsd']:6.2f}A Fe{gl['fe']:5.2f}  {orc['rmsd']:8.2f}")

with open(os.path.join(A25, "pick_comparison.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

n = len(rows)
print(f"\n{'ranker':16} {'median RMSD to co-folded':>26} {'<2A':>6} {'median Fe':>10}")
for k in ("icm_Score", "glide_gscore", "glide_emodel", "oracle"):
    v = agg[k]
    fe = (f"{st.median(agg['fe_icm']):10.2f}" if k == "icm_Score" else
          f"{st.median(agg['fe_glide']):10.2f}" if k == "glide_gscore" else " " * 10)
    print(f"{k:16} {st.median(v):26.2f} {sum(1 for x in v if x < 2):4}/{n} {fe}")
same = sum(1 for r in rows if r["icm_rank"] == r["glide_rank"])
print(f"\nGlide and ICM choose the same pose in {same}/{n} (chance ~{n/10:.0f})")
print("the two ligands ICM ejected:")
for lg in ("OCNT-2313429", "OCNT-2315379"):
    r = next((x for x in rows if x["ligand"] == lg), None)
    if r:
        print(f"  {lg}: ICM r{r['icm_rank']} Fe {r['icm_fe']:.2f}A  ->  "
              f"Glide r{r['glide_rank']} Fe {r['glide_fe']:.2f}A "
              f"(RMSD {r['icm_rmsd']:.2f} -> {r['glide_rmsd']:.2f})")
