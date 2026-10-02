import numpy as np, glob
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
print(f"{'ligand':16}{'min dist res1-22 to LIG':>26}{'to HEM Fe':>12}")
vals=[]
for f in sorted(glob.glob("submissions/06_af3_pliptm/submission/*.pdb"))[:6]:
    m=Chem.MolFromPDBFile(f,removeHs=True,sanitize=False); c=m.GetConformer()
    nt=[];lig=[];fe=None
    for a in m.GetAtoms():
        if a.GetSymbol()=="H": continue
        ri=a.GetPDBResidueInfo(); x=np.array(c.GetAtomPosition(a.GetIdx()))
        rn=ri.GetResidueName().strip()
        if rn=="LIG": lig.append(x)
        elif ri.GetName().strip()=="FE": fe=x
        elif rn!="HEM" and ri.GetResidueNumber()<=22: nt.append(x)
    if not nt:
        print("  no residues 1-22 found"); break
    nt=np.array(nt); lig=np.array(lig)
    d1=np.linalg.norm(nt[:,None]-lig[None],axis=2).min()
    d2=np.linalg.norm(nt-fe,axis=1).min()
    vals.append(d1)
    print(f"{f.split('/')[-1][:-4]:16}{d1:>26.1f}{d2:>12.1f}")
if vals: print(f"\nN-terminal 22 residues are {min(vals):.0f}-{max(vals):.0f} A from the ligand")
