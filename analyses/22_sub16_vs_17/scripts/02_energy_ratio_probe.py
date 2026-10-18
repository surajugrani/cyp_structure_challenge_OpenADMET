"""Does PoseBusters give a NUMERIC ligand internal-strain value for our structures?

The run logs show "Energy ratio module failed because UFF parameters missing" next to
"UFFTYPER: Unrecognized atom type: Fe4+2". But energy_ratio is a LIGAND-only check and
the ligand alone is passed as mol_pred, so the haem should be irrelevant - the warning
comes from elsewhere in the call. `full_report=True` exposes the underlying numbers
rather than just the pass/fail booleans.
"""
import csv, sys

sys.path.insert(0, "submissions/01_boltz2_default/scripts/organizer_code")
from posebusters_utils import assign_bond_orders_from_template
from posebusters import PoseBusters
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
          csv.DictReader(open("input/cyp-challenge-TEST-BLINDED_structures.csv"))
          if r["Molecule_Name"].strip()}
pb = PoseBusters(config="dock")
first = True
for n in ("OCNT-2312554", "OCNT-2313038", "OCNT-2395645", "OCNT-2395464"):
    pred = Chem.MolFromPDBFile(f"submissions/17_boltz2_ipde/submission/{n}.pdb",
                               removeHs=False, sanitize=False)
    frags = Chem.SplitMolByPDBResidues(pred)
    lig = frags.pop("LIG")
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    tpl = Chem.MolFromSmiles(smiles[n])
    try:
        lig = assign_bond_orders_from_template(tpl, lig)
    except Exception:
        pass
    r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot, full_report=True)
    row = r.iloc[0]
    if first:
        print("ALL columns with full_report=True:")
        for c in r.columns:
            print("   ", c)
        print()
        first = False
    hits = {c: row[c] for c in r.columns
            if any(k in c.lower() for k in ("energy", "strain", "ratio"))}
    print(f"{n}: {hits}")
