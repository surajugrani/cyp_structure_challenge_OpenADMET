import numpy as np, os, csv
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
cur=list(csv.DictReader(open("analyses/01_glide_score_vs_ic50/data/curated_set.csv")))
g5=list(csv.DictReader(open("analyses/05_full_pdb_redocking/results/redock_summary.csv")))
docked={r["pdb_id"] for r in g5 if r["status"]=="ok"}
out=[]
for r in cur:
    p=r["pdb_id"]
    f=f"analyses/06_icm_redocking/structures/{p}/{p}_crystal.sdf"
    rec=f"analyses/06_icm_redocking/structures/{p}/{p}_rec.pdb"
    if not (os.path.exists(f) and os.path.exists(rec)): continue
    lig=next(iter(Chem.SDMolSupplier(f,removeHs=True,sanitize=False)),None)
    m=Chem.MolFromPDBFile(rec,removeHs=True,sanitize=False)
    if lig is None or m is None: continue
    c=m.GetConformer(); fe=None
    for a in m.GetAtoms():
        ri=a.GetPDBResidueInfo()
        if ri and ri.GetName().strip()=="FE": fe=np.array(c.GetAtomPosition(a.GetIdx()))
    if fe is None: continue
    lc=lig.GetConformer()
    d=min(np.linalg.norm(np.array(lc.GetAtomPosition(i))-fe)
          for i,a in enumerate(lig.GetAtoms()) if a.GetSymbol()!="H")
    out.append({"pdb":p,"comp":r["comp_id"],"fe":round(float(d),2),
                "coord":bool(d<=2.6),"has_glide_grid":bool(p in docked)})
out.sort(key=lambda x:x["fe"])
non=[x for x in out if not x["coord"]]
co=[x for x in out if x["coord"] and x["has_glide_grid"]]
# spread the 5 coordinating across the distance range
idx=[0, len(co)//4, len(co)//2, 3*len(co)//4, len(co)-1]
pick_co=[co[i] for i in sorted(set(idx))][:5]
print("NON-coordinating (all 5):")
for x in non: print(f"  {x['pdb']} {x['comp']:>6}  Fe {x['fe']:>5}  grid={x['has_glide_grid']}")
print("\nCOORDINATING (5, spread over range):")
for x in pick_co: print(f"  {x['pdb']} {x['comp']:>6}  Fe {x['fe']:>5}  grid={x['has_glide_grid']}")
import json
json.dump({"coordinating":[x["pdb"] for x in pick_co],
           "non_coordinating":[x["pdb"] for x in non],
           "detail":non+pick_co},
          open("analyses/11_preprint_insights/results/metal_cons_set.json","w"), indent=1)
