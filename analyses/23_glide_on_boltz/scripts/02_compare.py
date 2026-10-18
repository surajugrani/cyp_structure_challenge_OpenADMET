"""Would Glide pick different Boltz samples than PL-iPTM, and does its score vary?

Context: submission 16 showed PL-iPTM selection on this ensemble is worth +0.0166
LDDT-PLI. Submissions 13 and 14 showed Glide and ICM scores LOSE 0.056 and 0.034 on the
AF3 ensemble. This asks whether Glide behaves any differently on Boltz's samples.

Agreement with PL-iPTM near 1 in 5 means Glide is choosing essentially at random relative
to the rule that works.
"""
import csv, glob, json, os, statistics as st
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
A23 = os.path.join(ROOT, "analyses", "23_glide_on_boltz", "results")

gl = defaultdict(dict)
nerr = 0
for r in csv.DictReader(open(os.path.join(A23, "glide_scores.csv"))):
    if r["status"] != "ok" or not r["gscore"] or abs(float(r["gscore"])) > 100:
        nerr += 1
        continue
    gl[r["ligand"]][r["model"]] = {t: (float(r[t]) if r.get(t) else None)
                                   for t in ("gscore", "emodel", "evdw", "lipo")}

conf = defaultdict(dict)
for d in sorted(glob.glob(os.path.join(ROOT, "submissions", "15_boltz2_tuned",
                                       "results", "*", "predictions", "*"))):
    lig = os.path.basename(d)
    for i in range(5):
        f = os.path.join(d, f"confidence_{lig}_model_{i}.json")
        if os.path.exists(f):
            j = json.load(open(f))
            conf[lig][f"model_{i}"] = {"pl": j["pair_chains_iptm"]["0"]["2"],
                                       "ipde": j["complex_ipde"]}

print(f"{nerr} samples unusable\n")
print(f"{'ligand':16} {'PL':>4} {'gscore':>7} {'emodel':>7} {'gscore spread':>14}")
rows, agree_g, agree_e = [], 0, 0
spreads = defaultdict(list)
for lig in sorted(gl):
    common = [m for m in gl[lig] if m in conf.get(lig, {})]
    if len(common) < 2:
        continue
    pl = max(common, key=lambda m: conf[lig][m]["pl"])
    gs = min(common, key=lambda m: gl[lig][m]["gscore"])
    em = min(common, key=lambda m: gl[lig][m]["emodel"])
    sp = (max(gl[lig][m]["gscore"] for m in common)
          - min(gl[lig][m]["gscore"] for m in common))
    spreads["gscore"].append(sp)
    spreads["pl"].append(max(conf[lig][m]["pl"] for m in common)
                         - min(conf[lig][m]["pl"] for m in common))
    agree_g += gs == pl
    agree_e += em == pl
    rows.append({"ligand": lig, "n": len(common), "pl_pick": pl, "gscore_pick": gs,
                 "emodel_pick": em, "gscore_spread": round(sp, 3),
                 "gscore_at_pl_pick": round(gl[lig][pl]["gscore"], 3),
                 "gscore_best": round(gl[lig][gs]["gscore"], 3)})
    print(f"{lig:16} {pl.replace('model_',''):>4} {gs.replace('model_',''):>7} "
          f"{em.replace('model_',''):>7} {sp:14.3f}")

with open(os.path.join(A23, "pick_comparison.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
n = len(rows)
print(f"\nGlide gscore agrees with the PL-iPTM pick: {agree_g}/{n}   "
      f"emodel: {agree_e}/{n}   (chance is about {n/5:.0f}/{n})")
print(f"median within-ligand spread: gscore {st.median(spreads['gscore']):.3f} kcal/mol, "
      f"PL-iPTM {st.median(spreads['pl']):.3f}")
