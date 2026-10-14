"""Base rate: in real CYP3A4 crystal structures, do all capable ligands coordinate Fe?

The question this answers: if a predicted pose has an Fe-capable ligand that does NOT
coordinate, how confident can we be that the prediction is wrong?

That inference is only valid if crystal structures essentially always coordinate when the
chemistry allows. If a sizeable fraction of real structures have a capable ligand sitting
3-5 A away, then "capable but not coordinating" is a normal binding mode and tells us
nothing about whether a prediction is right.

Measures, for each PDB structure in analysis 05: whether the crystal ligand carries a
nitrogen donor (same corrected SMARTS as analysis 16 script 04), and the minimum distance
from the haem iron to any crystal ligand heavy atom.
"""
import csv, glob, json, math, os

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
A05 = "analyses/05_full_pdb_redocking"
N_DONOR = Chem.MolFromSmarts(
    "[$([nX2;!$([nH])]),"
    "$([NX3;H0,H1,H2;!$(N[C,S]=[O,S,N]);!$(N[a])&!$(N=*)]),"
    "$([NX3;H0,H1;$(Nc);!$(N[C,S]=[O,S,N])]),"
    "$([NX1]#[CX2])]")

# ligand_smiles.json is keyed by 3-letter component id, not PDB id, and each value is a
# dict - the PDB -> comp_id mapping lives in redock_summary.csv.
_raw = json.load(open(f"{A05}/data/ligand_smiles.json"))
_comp = {k: (v.get("smiles") if isinstance(v, dict) else v) for k, v in _raw.items()}
smiles = {}
for r in csv.DictReader(open(f"{A05}/results/redock_summary.csv")):
    c = r.get("comp_id", "").strip()
    if c in _comp:
        smiles[r["pdb_id"].strip()] = _comp[c]


def fe_from_pdb(path):
    """All haem iron positions in the deposited structure."""
    out = []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[:2]).upper()
        if el == "FE":
            out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


rows = []
for d in sorted(glob.glob(f"{A05}/structures/*")):
    pdb = os.path.basename(d)
    sdf = os.path.join(d, f"{pdb}_crystal_lig.sdf")
    raw = os.path.join(d, f"{pdb}.pdb")
    if not (os.path.exists(sdf) and os.path.exists(raw)):
        continue
    m = next(iter(Chem.SDMolSupplier(sdf, removeHs=True, sanitize=False)), None)
    if m is None or m.GetNumConformers() == 0:
        continue
    fes = fe_from_pdb(raw)
    if not fes:
        continue
    conf = m.GetConformer()
    best = min(math.dist(fe, tuple(conf.GetAtomPosition(a.GetIdx())))
               for a in m.GetAtoms() if a.GetAtomicNum() > 1 for fe in fes)
    smi = smiles.get(pdb)
    ref = Chem.MolFromSmiles(smi) if smi else None
    nd = len(ref.GetSubstructMatches(N_DONOR)) if ref else None
    rows.append({"pdb": pdb, "smiles": smi or "", "n_donors": nd,
                 "fe_dist": round(best, 2)})

with open("analyses/19_crystal_fe_basterate/results/crystal_fe.csv", "w",
          newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

known = [r for r in rows if r["n_donors"] is not None]
cap = [r for r in known if r["n_donors"] > 0]
nocap = [r for r in known if r["n_donors"] == 0]
co = [r for r in cap if r["fe_dist"] <= 2.5]
far = [r for r in cap if r["fe_dist"] > 2.5]

print(f"{len(rows)} crystal structures with a ligand and a haem iron "
      f"({len(known)} with a known SMILES)\n")
print(f"ligands WITH an N donor:    {len(cap)}")
print(f"  coordinate (<=2.5 A):     {len(co)}  ({100*len(co)/len(cap):.0f} %)")
print(f"  do NOT coordinate:        {len(far)}  ({100*len(far)/len(cap):.0f} %)")
print(f"ligands with NO N donor:    {len(nocap)}")
print(f"  coordinate anyway:        {sum(r['fe_dist'] <= 2.5 for r in nocap)}")

if far:
    import statistics as st
    print(f"\nnon-coordinating but capable: median Fe dist "
          f"{st.median(r['fe_dist'] for r in far):.2f} A, range "
          f"{min(r['fe_dist'] for r in far):.2f}-{max(r['fe_dist'] for r in far):.2f}")
    print("  " + ", ".join(f"{r['pdb']}({r['fe_dist']})" for r in sorted(
        far, key=lambda x: x["fe_dist"])[:15]))
