"""Submission 19: per ligand, whichever of submission 06 or 18 Glide scores better.

06 is AF3 + PL-iPTM selection (0.4291 / 4.0647 / 0.6985) - best LDDT-PLI of the project.
18 is tuned Boltz + Glide selection (0.4174 / 3.7523 / 0.8047) - best BiSyRMSD.

Result: 14 structures from 18, 6 from 06.

The structures are copied verbatim from each submission - both are already
PrepWizard-minimised and have passed validation - so this only re-runs the organizers'
structural checks and PoseBusters before zipping.

TWO CAVEATS worth stating plainly.

1. **The two gscores come from different receptors.** 06's is computed on AF3's protein,
   18's on Boltz's, each with its own grid. Glide scores are not calibrated for
   cross-receptor comparison, so "better gscore" here is a weaker claim than it looks.

2. **The two submissions disagree enormously** - median ligand RMSD 4.15 A, max 7.71 A,
   17 of 20 beyond 2 A - and analysis of that disagreement found NO relationship between
   how far apart two poses are and which one Glide prefers (Pearson +0.110, Spearman
   +0.018, n=19). So this rule is choosing between genuinely different answers using a
   score that does not track the difference.

Against that, Glide selection WITHIN the Boltz ensemble produced the best BiSyRMSD of the
project (submission 18), overturning an earlier recommendation not to try it. That is the
case for running this one too: predictions about what Glide can and cannot do have been
wrong in both directions here, and the leaderboard is the only arbiter.

OCNT-2395789 takes the Boltz structure by instruction: its chosen AF3 sample is one of
three that carried Glide's 10000.0 "no score" sentinel, so no comparable score exists.
"""
import csv, os, sys, zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "submission")
ZIP = os.path.join(HERE, "submission_hybrid_glide.zip")
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
