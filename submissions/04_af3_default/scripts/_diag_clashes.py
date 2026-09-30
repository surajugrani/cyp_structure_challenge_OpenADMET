import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
bad=["OCNT-2312554","OCNT-2312966","OCNT-2313300","OCNT-2315379","OCNT-2395342","OCNT-2535422","OCNT-2535519"]
good=["OCNT-2313038","OCNT-2395436"]
def rep(n,tag):
    m=Chem.MolFromPDBFile(f"submissions/04_af3_default/submission/{n}.pdb",removeHs=True,sanitize=False)
    c=m.GetConformer(); L=[];P=[];FE=None
    for a in m.GetAtoms():
        if a.GetSymbol()=="H": continue
        ri=a.GetPDBResidueInfo(); x=np.array(c.GetAtomPosition(a.GetIdx()))
        if ri.GetResidueName().strip()=="LIG": L.append(x)
        else:
            P.append((x,ri.GetResidueName().strip(),ri.GetName().strip()))
            if ri.GetName().strip()=="FE": FE=x
    L=np.array(L); Px=np.array([p[0] for p in P])
    d=np.linalg.norm(L[:,None]-Px[None],axis=2); i,j=np.unravel_index(d.argmin(),d.shape)
    dfe=np.linalg.norm(L-FE,axis=1).min() if FE is not None else float('nan')
    print(f"{tag:5}{n:16}{P[j][1]+' '+P[j][2]:>22}{d.min():>8.2f}{dfe:>9.2f}")
print(f"{'':5}{'ligand':16}{'closest contact':>22}{'dist':>8}{'lig-Fe':>9}")
for n in bad: rep(n,"FAIL")
for n in good: rep(n,"pass")
