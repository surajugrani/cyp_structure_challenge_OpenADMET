"""Which submitted ligands CAN coordinate Fe but are not doing so in the predicted pose?"""
import csv, os
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

# Fe-ligating motifs, same SMARTS set as analysis 07
M = {"imidazole": "c1cnc[nH0,nH]1", "triazole": "c1ncn[nH0,nH]1",
     "tetrazole": "c1nnn[nH0,nH]1", "pyridine": "n1ccccc1",
     "pyrimidine": "n1cnccc1", "thiazole": "c1scnc1", "oxazole": "c1ocnc1",
     "nitrile": "[NX1]#[CX2]", "primary_amine": "[NX3;H2;!$(NC=O)]",
     "thioether": "[#16X2]"}
COORD = 2.6     # generous upper bound for Fe-N/S coordination

rows = list(csv.DictReader(open("input/cyp-challenge-TEST-BLINDED_structures.csv")))
out = []
for r in rows:
    n, smi = r["Molecule_Name"].strip(), r["SMILES"].strip()
    if not n:
        continue
    mol = Chem.MolFromSmiles(smi)
    motifs = [k for k, s in M.items() if mol.HasSubstructMatch(Chem.MolFromSmarts(s))]

    f = f"submissions/07_af3_templates/submission/{n}.pdb"
    m = Chem.MolFromPDBFile(f, removeHs=True, sanitize=False)
    c = m.GetConformer()
    FE = None; best = (1e9, None)
    lig = []
    for a in m.GetAtoms():
        if a.GetSymbol() == "H":
            continue
        ri = a.GetPDBResidueInfo(); x = np.array(c.GetAtomPosition(a.GetIdx()))
        if ri.GetResidueName().strip() == "LIG":
            lig.append((x, a.GetSymbol(), ri.GetName().strip()))
        elif ri.GetName().strip() == "FE":
            FE = x
    for x, el, nm in lig:
        d = float(np.linalg.norm(x - FE))
        if d < best[0]:
            best = (d, f"{nm}({el})")
    # nearest N or S specifically - the atoms that actually ligate
    hetero = [(float(np.linalg.norm(x - FE)), f"{nm}({el})")
              for x, el, nm in lig if el in ("N", "S")]
    hb = min(hetero) if hetero else (float("nan"), "-")
    out.append({"ligand": n, "motifs": ",".join(motifs) or "-",
                "d_any": best[0], "atom_any": best[1],
                "d_NS": hb[0], "atom_NS": hb[1],
                "coordinating": best[0] <= COORD})

print(f"{'ligand':16}{'Fe-ligating motif':28}{'closest':>9}{'atom':>12}"
      f"{'closest N/S':>13}{'atom':>12}  coord?")
for r in sorted(out, key=lambda x: (x["motifs"] == "-", x["d_any"])):
    print(f"{r['ligand']:16}{r['motifs'][:27]:28}{r['d_any']:>9.2f}{r['atom_any']:>12}"
          f"{r['d_NS']:>13.2f}{r['atom_NS']:>12}  {'YES' if r['coordinating'] else 'no'}")

cap = [r for r in out if r["motifs"] != "-"]
miss = [r for r in cap if not r["coordinating"]]
print(f"\nligands with an Fe-ligating motif: {len(cap)}/20")
print(f"  of those, coordinating (<= {COORD} A): {len(cap)-len(miss)}")
print(f"  of those, NOT coordinating:            {len(miss)}")
for r in miss:
    print(f"    {r['ligand']}  motif={r['motifs']}  closest Fe contact {r['d_any']:.2f} A "
          f"({r['atom_any']}), nearest N/S {r['d_NS']:.2f} A")
