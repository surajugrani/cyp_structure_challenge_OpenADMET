"""Which sample would each ICM score pick, and how far from submission 06's choice?

Rules compared, all selecting among a ligand's PoseBusters-PASSING prepped AF3 samples:
  Score         ICM's default VLS score (lowest wins)
  RTCNN         the neural-network rescore - best within-ligand signal in analysis 09
                (Spearman vs RMSD +0.334, against Glide gscore's +0.086)
  RTCNN+Hbond   analysis 09's preferred ICM model, as a simple sum of z-scores within
                each ligand's pool (the two terms are on different scales)
  mfScore       included only for completeness; weakest in analysis 09 (3/30 top-1)
"""
import csv, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
A20 = os.path.join(ROOT, "analyses", "20_icm_score_af3", "results")
S04 = os.path.join(ROOT, "submissions", "04_af3_default", "glide_rescoring", "results")
A16 = os.path.join(ROOT, "analyses", "16_fe_distance_swap", "results")

pbpass, fe = set(), {}
for r in csv.DictReader(open(os.path.join(S04, "pb_after_prepwizard.csv"))):
    if r["after_pass"] == "True":
        pbpass.add((r["ligand"], r["sample"]))
    fe[(r["ligand"], r["sample"])] = float(r["fe_dist_after"])
sel06 = {r["ligand"]: r["chosen"] for r in csv.DictReader(
    open(os.path.join(ROOT, "submissions", "06_af3_pliptm", "results", "selection.csv")))}
don = {r["ligand"]: int(r["n_donors"]) > 0
       for r in csv.DictReader(open(os.path.join(A16, "donor_recount.csv")))}

pool = defaultdict(list)
nerr = 0
for r in csv.DictReader(open(os.path.join(A20, "icm_sample_scores.csv"))):
    if r["status"] != "ok":
        nerr += 1
        continue
    if (r["ligand"], r["sample"]) not in pbpass:
        continue
    try:
        pool[r["ligand"]].append({"sample": r["sample"],
                                  "Score": float(r["Score"]),
                                  "RTCNN": float(r["RTCNN"]),
                                  "mfScore": float(r["mfScore"]),
                                  "Hbond": float(r["Hbond"]),
                                  "fe": fe[(r["ligand"], r["sample"])]})
    except (ValueError, KeyError):
        nerr += 1
print(f"{nerr} samples unusable; {sum(len(v) for v in pool.values())} scored and "
      f"PB-passing over {len(pool)} ligands\n")


def z(vals):
    m = sum(vals) / len(vals)
    sd = (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5
    return [0.0] * len(vals) if sd == 0 else [(v - m) / sd for v in vals]


rules = {}
for lig, cs in pool.items():
    zr = z([c["RTCNN"] for c in cs])
    zh = z([c["Hbond"] for c in cs])
    for c, a, b in zip(cs, zr, zh):
        c["RTCNN_Hbond"] = a + b
    rules.setdefault("Score", {})[lig] = min(cs, key=lambda c: c["Score"])
    rules.setdefault("RTCNN", {})[lig] = min(cs, key=lambda c: c["RTCNN"])
    rules.setdefault("RTCNN_Hbond", {})[lig] = min(cs, key=lambda c: c["RTCNN_Hbond"])
    rules.setdefault("mfScore", {})[lig] = min(cs, key=lambda c: c["mfScore"])

cap = [l for l in pool if don.get(l)]
print(f"{'rule':14} {'differs from 06':>16} {'coordinating of capable':>26} {'median Fe':>10}")
base_coord = sum(1 for l in cap if fe[(l, sel06[l])] <= 2.5)
import statistics as st
print(f"{'06 PL-iPTM':14} {'-':>16} {base_coord:>19}/{len(cap)} "
      f"{st.median(fe[(l, sel06[l])] for l in pool):10.2f}")
rows = []
for name, pick in rules.items():
    d = sum(1 for l in pick if pick[l]["sample"] != sel06[l])
    co = sum(1 for l in cap if pick[l]["fe"] <= 2.5)
    print(f"{name:14} {d:13}/{len(pick)} {co:>19}/{len(cap)} "
          f"{st.median(pick[l]['fe'] for l in pick):10.2f}")
    for l in sorted(pick):
        rows.append({"rule": name, "ligand": l, "chosen": pick[l]["sample"],
                     "sub06_chose": sel06[l],
                     "differs": pick[l]["sample"] != sel06[l],
                     "fe_dist": pick[l]["fe"], "can_coordinate": don.get(l),
                     "Score": pick[l]["Score"], "RTCNN": pick[l]["RTCNN"],
                     "Hbond": pick[l]["Hbond"]})
with open(os.path.join(A20, "rule_comparison.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("\nwrote results/rule_comparison.csv")
