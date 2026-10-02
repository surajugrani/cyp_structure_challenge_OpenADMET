import numpy as np, glob, os, csv
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
cur=list(csv.DictReader(open("analyses/01_glide_score_vs_ic50/data/curated_set.csv")))
n_coord=0; rows=[]
for r in cur:
    p=r["pdb_id"]
    f=f"analyses/06_icm_redocking/structures/{p}/{p}_crystal.sdf"
    rec=f"analyses/06_icm_redocking/structures/{p}/{p}_rec.pdb"
    if not (os.path.exists(f) and os.path.exists(rec)): 
        rows.append((p,None)); continue
    lig=next(iter(Chem.SDMolSupplier(f,removeHs=True,sanitize=False)),None)
    m=Chem.MolFromPDBFile(rec,removeHs=True,sanitize=False)
    if lig is None or m is None: rows.append((p,None)); continue
    fe=None; c=m.GetConformer()
    for a in m.GetAtoms():
        ri=a.GetPDBResidueInfo()
        if ri and ri.GetName().strip()=="FE": fe=np.array(c.GetAtomPosition(a.GetIdx()))
    if fe is None: rows.append((p,None)); continue
    lc=lig.GetConformer()
    d=min(np.linalg.norm(np.array(lc.GetAtomPosition(i))-fe)
          for i,a in enumerate(lig.GetAtoms()) if a.GetSymbol()!="H")
    rows.append((p,d)); n_coord += d<=2.6
ok=[d for _,d in rows if d is not None]
print(f"measured {len(ok)}/{len(cur)} | coordinating (<=2.6 A): {n_coord}")
print("distances:", " ".join(f"{d:.2f}" for _,d in sorted(rows,key=lambda x:(x[1] is None,x[1])) if d is not None))
