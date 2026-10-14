"""If Glide's docking score picked the best of 5 AF3 samples, what would change?

Compares against submission 06 (AF3, PL-iPTM selection) - our best leaderboard result and
the natural baseline, since both rules select from the SAME 5 AF3 samples per ligand.
(Comparing against submission 03 would be meaningless: that is Boltz-2, a different
method, so all 20 structures differ by construction.)

Reports, for both selection rules: which sample is chosen, its Fe-ligand distance, and
how many of the ligands chemically ABLE to coordinate iron actually do so.

"Can coordinate" = SMARTS [nX2;!$([nH])] - an aromatic nitrogen with a free lone pair
(pyridine, imidazole N3, triazole), excluding pyrrole-type NH.
"""
import csv

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
PATT = Chem.MolFromSmarts("[nX2;!$([nH])]")
S04 = "submissions/04_af3_default/glide_rescoring/results"

smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
          csv.DictReader(open("input/cyp-challenge-TEST-BLINDED_structures.csv"))
          if r["Molecule_Name"].strip()}
sel06 = {r["ligand"]: r["chosen"] for r in
         csv.DictReader(open("submissions/06_af3_pliptm/results/selection.csv"))}

pb = {}
for r in csv.DictReader(open(f"{S04}/pb_after_prepwizard.csv")):
    pb[(r["ligand"], r["sample"])] = {"pass": r["after_pass"] == "True",
                                      "fe": float(r["fe_dist_after"])}
gl = {}
for r in csv.DictReader(open(f"{S04}/glide_scores.csv")):
    if r["status"] == "ok" and r["gscore"] and abs(float(r["gscore"])) <= 100:
        gl[(r["ligand"], r["sample"])] = float(r["gscore"])

rows = []
for n in sorted(smiles):
    m = Chem.MolFromSmiles(smiles[n])
    can = bool(m and m.HasSubstructMatch(PATT))
    pool = [s for (lg, s) in pb if lg == n and pb[(lg, s)]["pass"]]
    scored = [s for s in pool if (n, s) in gl]
    if not scored:
        continue
    g_pick = min(scored, key=lambda s: gl[(n, s)])
    p_pick = sel06[n]
    rows.append({"ligand": n, "can_coordinate": can,
                 "pliptm_06_sample": p_pick, "pliptm_06_fe": pb[(n, p_pick)]["fe"],
                 "glide_sample": g_pick, "glide_fe": pb[(n, g_pick)]["fe"],
                 "glide_gscore": round(gl[(n, g_pick)], 3),
                 "differs": g_pick != p_pick})

with open("analyses/16_fe_distance_swap/results/glide_selection_fe.csv", "w",
          newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

print(f"{'ligand':16} {'can':>4} {'06 PL-iPTM':>22} {'Glide best':>22}  diff")
for r in rows:
    print(f"{r['ligand']:16} {'yes' if r['can_coordinate'] else 'no':>4} "
          f"{r['pliptm_06_sample'][7:]:>10} {r['pliptm_06_fe']:6.2f} A   "
          f"{r['glide_sample'][7:]:>10} {r['glide_fe']:6.2f} A  "
          f"{'*' if r['differs'] else ''}")

d = sum(r["differs"] for r in rows)
cap = [r for r in rows if r["can_coordinate"]]
print(f"\ndiffers from submission 06: {d}/{len(rows)}")
print(f"\nligands able to coordinate Fe: {len(cap)}/{len(rows)}")
print(f"  coordinate (<=2.5 A) with PL-iPTM selection (06): "
      f"{sum(r['pliptm_06_fe'] <= 2.5 for r in cap)}/{len(cap)}")
print(f"  coordinate (<=2.5 A) with Glide selection:        "
      f"{sum(r['glide_fe'] <= 2.5 for r in cap)}/{len(cap)}")
print(f"\nall {len(rows)}:  06 {sum(r['pliptm_06_fe'] <= 2.5 for r in rows)}   "
      f"Glide {sum(r['glide_fe'] <= 2.5 for r in rows)}")
import statistics as st
print(f"median Fe distance:  06 {st.median(r['pliptm_06_fe'] for r in rows):.2f} A   "
      f"Glide {st.median(r['glide_fe'] for r in rows):.2f} A")
