"""Which ligands CAN coordinate Fe, and which actually do — AF3 (06) vs ICM redock (12).

The type II / type I split only means something against what each ligand is chemically
capable of. A ligand with no free aromatic nitrogen should not coordinate, so failing to
is correct, not an error. This counts coordination only among the ligands that could.

"Can coordinate" = SMARTS [nX2;!$([nH])], an aromatic nitrogen with two connections and
no hydrogen, i.e. a free lone pair pointing out of the ring - pyridine, imidazole N3,
triazole. Excludes pyrrole-type NH and 3-substituted aromatic N, whose lone pairs are in
the ring system.
"""
import csv

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
PATT = Chem.MolFromSmarts("[nX2;!$([nH])]")

smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
          csv.DictReader(open("input/cyp-challenge-TEST-BLINDED_structures.csv"))
          if r["Molecule_Name"].strip()}
icm = {r["ligand"]: float(r["fe_dist"]) for r in
       csv.DictReader(open("submissions/12_icm_redock_af3/results/selection.csv"))}
af3 = {r["ligand"]: float(r["af3_06_fe"]) for r in
       csv.DictReader(open("analyses/16_fe_distance_swap/results/fe_distance.csv"))}

rows = []
print(f"{'ligand':16} {'can':>4} {'AF3':>6} {'ICM':>6}  note")
for n in sorted(smiles):
    m = Chem.MolFromSmiles(smiles[n])
    can = bool(m and m.HasSubstructMatch(PATT))
    a, i = af3[n], icm[n]
    note = []
    if can and a <= 2.5 and i > 2.5:
        note.append("ICM lost coordination")
    if i > 6:
        note.append("ICM far from haem")
    rows.append({"ligand": n, "can_coordinate": can, "af3_06_fe": a, "icm_12_fe": i,
                 "note": "; ".join(note)})
    print(f"{n:16} {'yes' if can else 'no':>4} {a:6.2f} {i:6.2f}  {'; '.join(note)}")

with open("analyses/16_fe_distance_swap/results/coordination.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

cap = [r for r in rows if r["can_coordinate"]]
print(f"\nligands able to coordinate (free aromatic N): {len(cap)}/20")
print(f"  actually coordinate (<=2.5 A) in AF3 (06): "
      f"{sum(r['af3_06_fe'] <= 2.5 for r in cap)}/{len(cap)}")
print(f"  actually coordinate (<=2.5 A) in ICM (12): "
      f"{sum(r['icm_12_fe'] <= 2.5 for r in cap)}/{len(cap)}")
print(f"  lost by ICM: {sum(r['af3_06_fe'] <= 2.5 < r['icm_12_fe'] for r in cap)}")
nc = [r for r in rows if not r["can_coordinate"]]
print(f"\nligands unable to coordinate: {len(nc)}/20  "
      f"(AF3 nonetheless <=2.5 A for {sum(r['af3_06_fe'] <= 2.5 for r in nc)}, "
      f"ICM for {sum(r['icm_12_fe'] <= 2.5 for r in nc)})")
print(f"\nall 20: AF3 {sum(r['af3_06_fe'] <= 2.5 for r in rows)}/20   "
      f"ICM {sum(r['icm_12_fe'] <= 2.5 for r in rows)}/20")
far = [(r["ligand"], r["icm_12_fe"]) for r in rows if r["icm_12_fe"] > 6]
print(f"ICM >6 A from haem: {far or 'none'}")
for t in ("OCNT-2313429", "OCNT-2315379"):
    r = next(x for x in rows if x["ligand"] == t)
    print(f"{t}: can={r['can_coordinate']}  AF3 {r['af3_06_fe']} A  ICM {r['icm_12_fe']} A")
