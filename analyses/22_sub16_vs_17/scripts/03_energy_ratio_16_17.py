"""Ligand internal-strain (PoseBusters energy_ratio) for submissions 16 and 17.

energy_ratio = pose UFF energy / mean energy of a generated conformer ensemble. Lower
means the pose is closer to a relaxed conformation; PoseBusters fails above ~7.

The two submissions differ on 8 of 20 structures. If strain is informative, picking the
lower-strain pose per ligand gives a third selection rule - and one built on a quantity
completely independent of iPTM and iPDE, which are both confidence estimates rather than
physical plausibility.
"""
import csv, sys

sys.path.insert(0, "submissions/01_boltz2_default/scripts/organizer_code")
from posebusters_utils import assign_bond_orders_from_template
from posebusters import PoseBusters
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
SUBS = {"16": "submissions/16_boltz2_pliptm/submission",
        "17": "submissions/17_boltz2_ipde/submission"}
smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
          csv.DictReader(open("input/cyp-challenge-TEST-BLINDED_structures.csv"))
          if r["Molecule_Name"].strip()}
differs = {r["ligand"] for r in csv.DictReader(
    open("analyses/22_sub16_vs_17/results/sub16_vs_17.csv")) if r["differs"] == "True"}
pb = PoseBusters(config="dock")


def strain(path, smi):
    pred = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    frags = Chem.SplitMolByPDBResidues(pred)
    lig = frags.pop("LIG")
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    tpl = Chem.MolFromSmiles(smi)
    try:
        lig = assign_bond_orders_from_template(tpl, lig)
    except Exception:
        pass
    r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot, full_report=True).iloc[0]
    g = lambda c: (None if c not in r.index or r[c] is None else float(r[c]))
    return g("energy_ratio"), g("mol_pred_energy"), g("ensemble_avg_energy")


rows = []
for n in sorted(smiles):
    rec = {"ligand": n, "differs": n in differs}
    for k, d in SUBS.items():
        er, e, ens = strain(f"{d}/{n}.pdb", smiles[n])
        rec[f"er_{k}"] = None if er is None else round(er, 3)
        rec[f"energy_{k}"] = None if e is None else round(e, 1)
    if rec["er_16"] is not None and rec["er_17"] is not None:
        rec["lower"] = "16" if rec["er_16"] <= rec["er_17"] else "17"
        rec["delta"] = round(rec["er_17"] - rec["er_16"], 3)
    rows.append(rec)
    print(f"{n:16} {'DIFF' if rec['differs'] else '    '} "
          f"er16 {str(rec['er_16']):>7}  er17 {str(rec['er_17']):>7}  "
          f"lower={rec.get('lower','-')}", flush=True)

with open("analyses/22_sub16_vs_17/results/energy_ratio_16_17.csv", "w",
          newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

d = [r for r in rows if r["differs"] and r.get("lower")]
print(f"\nof the {len(d)} differing structures: "
      f"16 lower in {sum(1 for r in d if r['lower']=='16')}, "
      f"17 lower in {sum(1 for r in d if r['lower']=='17')}")
import statistics as st
for k in ("16", "17"):
    v = [r[f"er_{k}"] for r in rows if r[f"er_{k}"] is not None]
    print(f"submission {k}: median energy_ratio {st.median(v):.3f}, max {max(v):.3f}")
best = [r["er_16"] if r["lower"] == "16" else r["er_17"] for r in rows if r.get("lower")]
print(f"min-strain hybrid: median {st.median(best):.3f}, max {max(best):.3f}")
