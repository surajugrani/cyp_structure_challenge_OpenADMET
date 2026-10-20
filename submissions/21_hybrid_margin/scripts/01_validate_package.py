"""Submission 21: submission 19, but keep AF3 only where it wins Glide by >0.5 kcal/mol.

Submission 19 (0.4331 / 3.6754 / 0.7823 - best LDDT-PLI and BiSyRMSD of the project)
took whichever of submission 06 or 18 Glide scored better, giving 14 Boltz / 6 AF3.
Three of those six AF3 wins were on a margin small enough to be noise:

    OCNT-2395712   0.07 kcal/mol
    OCNT-2395503   0.08
    OCNT-2395457   0.31
    OCNT-2315379   1.01   <- kept
    OCNT-2312966   1.47   <- kept
    OCNT-2535519   2.15   <- kept

Reverting the first three to Boltz gives 17 Boltz / 3 AF3 and changes exactly 3 of 20
structures.

THE BET: AF3 structures cost pocket quality (LDDT-LP 0.6985 against Boltz's 0.8047),
so importing one for a 0.07 kcal/mol score advantage is paying real LDDT-LP for noise.
The mixing model - accurate to ~0.004 on both submissions 19 and 20 - predicts
(17 x 0.8047 + 3 x 0.6985)/20 = 0.789, i.e. about +0.016 LDDT-LP over 19.

THE RISK: submission 19 beat submission 18 by 0.0157 LDDT-PLI precisely because of its
AF3 content. Halving that content could hand some of the gain back, and LDDT-PLI is
the primary metric while LDDT-LP is not.

No new compute - structures are copied verbatim from submissions 06 and 18, both
already PrepWizard-minimised and validated.
"""

import csv, os, sys, zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "submission")
ZIP = os.path.join(HERE, "submission_hybrid_margin.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")


def struct_errors(lines, name, smi, tmp):
    e = []
    at = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if len({l[21] for l in at}) > 3:
        e.append("more than 3 chains")
    lig = [l for l in at if l[17:20].strip() == "LIG"]
    if not lig:
        e.append("no LIG residue")
    elif len({(l[21], l[22:27]) for l in lig}) > 1:
        e.append("more than one LIG residue")
    nhem = len([l for l in at if l[17:20].strip() == "HEM"])
    if nhem != 43:
        e.append(f"HEM has {nhem} heavy atoms")
    per = defaultdict(list)
    for l in at:
        per[(l[21], l[22:27])].append(l[12:16].strip())
    if any(len(v) != len(set(v)) for v in per.values()):
        e.append("duplicate atom names in a residue")
    s = [l[6:11] for l in at]
    if len(s) != len(set(s)):
        e.append("duplicate atom serials")
    p = os.path.join(tmp, f"{name}_lig.pdb")
    open(p, "w").writelines(lig + ["END\n"])
    ref = Chem.MolFromSmiles(smi)
    probe = Chem.MolFromPDBFile(p, removeHs=True, sanitize=False)
    if ref is None or probe is None:
        e.append("ligand parse failure")
    else:
        try:
            AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except ValueError:
            e.append("ligand connectivity != SMILES")
    return e


def pb_fail(path, smi, pb):
    pred = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    if pred is None:
        return ["unparseable"]
    frags = Chem.SplitMolByPDBResidues(pred)
    lig = frags.pop("LIG", None)
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    if lig is None or prot is None:
        return ["no LIG or protein"]
    tpl = Chem.MolFromSmiles(smi)
    try:
        lig = assign_bond_orders_from_template(tpl, lig)
    except Exception:
        pass
    r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot)
    bc = [c for c in r.columns if r[c].dtype == bool]
    return sorted(r[bc].iloc[0].index[~r[bc].iloc[0]].tolist())


def main():
    tmp = os.path.join(HERE, "results", "_tmp")
    os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    sel = {r["ligand"]: r["source"] for r in csv.DictReader(
        open(os.path.join(HERE, "results", "selection.csv")))}
    pb = PoseBusters(config="dock")
    ns = npb = 0
    for n in sorted(smiles):
        p = os.path.join(OUT, f"{n}.pdb")
        se = struct_errors(open(p).readlines(), n, smiles[n], tmp)
        pf = pb_fail(p, smiles[n], pb)
        ns += not se
        npb += not pf
        print(f"  {n}: {sel.get(n,'?')}" + (f"  struct {se} | PB {pf}" if (se or pf) else ""))
    print(f"\nFINAL structural {ns}/20   PoseBusters {npb}/20")
    if ns == 20 and npb == 20:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for n in sorted(smiles):
                z.write(os.path.join(OUT, f"{n}.pdb"), f"{n}.pdb")
        print(f"wrote {ZIP}")
    else:
        print("zip not written - failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
