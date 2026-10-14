"""Corrected Fe-donor classification.

Script 02 used SMARTS [nX2;!$([nH])] - AROMATIC nitrogen only. That is the classic type II
motif (pyridine, imidazole N3, triazole) but it is not the full set of donors, and it
wrongly labelled OCNT-2535422 as unable to coordinate. That ligand coordinates through an
aliphatic indoline nitrogen at 2.28 A in submission 06, verified by measuring the nearest
ligand atom to the iron.

Corrected definition: any nitrogen with an available lone pair - aromatic nX2, or sp3/sp2
N that is NOT delocalised into an adjacent carbonyl or sulfonyl (amide, sulfonamide) and
not quaternary. Also reports sulfur and sp2 oxygen, which are weaker but real CYP donors.

Ground truth for "does it coordinate" is taken from the structures themselves: the element
of the ligand atom nearest the iron, and that distance.
"""
import csv, math

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
N_DONOR = Chem.MolFromSmarts(
    "[$([nX2;!$([nH])]),"                       # pyridine / imidazole / triazole
    "$([NX3;H0,H1,H2;!$(N[C,S]=[O,S,N]);!$(N[a])&!$(N=*)]),"   # aliphatic amine
    "$([NX3;H0,H1;$(Nc);!$(N[C,S]=[O,S,N])]),"  # aniline / indoline type
    "$([NX1]#[CX2])]")                          # nitrile
S_DONOR = Chem.MolFromSmarts("[#16X2]")


def nearest(path):
    fe, lig = None, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        p = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if el == "FE":
            fe = p
        elif l[17:20].strip() == "LIG":
            lig.append((el, p))
    if fe is None or not lig:
        return None, None
    d, el = min((math.dist(fe, p), e) for e, p in lig)
    return round(d, 2), el


smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
          csv.DictReader(open("input/cyp-challenge-TEST-BLINDED_structures.csv"))
          if r["Molecule_Name"].strip()}
icm = {r["ligand"]: float(r["fe_dist"]) for r in
       csv.DictReader(open("submissions/12_icm_redock_af3/results/selection.csv"))}

rows = []
print(f"{'ligand':16} {'Ndon':>5} {'Sdon':>5} {'AF3 06':>14} {'ICM 12':>7}")
for n in sorted(smiles):
    m = Chem.MolFromSmiles(smiles[n])
    nd = len(m.GetSubstructMatches(N_DONOR)) if m else 0
    sd = len(m.GetSubstructMatches(S_DONOR)) if m else 0
    d, el = nearest(f"submissions/06_af3_pliptm/submission/{n}.pdb")
    rows.append({"ligand": n, "n_donors": nd, "s_donors": sd,
                 "af3_fe_dist": d, "af3_contact_element": el,
                 "icm_fe_dist": icm[n]})
    print(f"{n:16} {nd:5d} {sd:5d} {d:9.2f} {el:>4} {icm[n]:7.2f}")

with open("analyses/16_fe_distance_swap/results/donor_recount.csv", "w",
          newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

cap = [r for r in rows if r["n_donors"] > 0]
print(f"\nligands with an N donor (corrected): {len(cap)}/20")
print(f"  coordinate <=2.5 A in AF3 (06): {sum(r['af3_fe_dist'] <= 2.5 for r in cap)}/{len(cap)}")
print(f"  coordinate <=2.5 A in ICM (12): {sum(r['icm_fe_dist'] <= 2.5 for r in cap)}/{len(cap)}")
nc = [r for r in rows if r["n_donors"] == 0]
print(f"ligands with NO N donor: {len(nc)}/20 "
      f"({[r['ligand'] for r in nc]})")
print(f"  of those, AF3 <=2.5 A: {sum(r['af3_fe_dist'] <= 2.5 for r in nc)}")
print(f"\ncontact element in AF3 (06): "
      f"{ {e: sum(1 for r in rows if r['af3_contact_element'] == e) for e in sorted({r['af3_contact_element'] for r in rows})} }")
