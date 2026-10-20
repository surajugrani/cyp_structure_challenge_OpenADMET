"""Submission 20: best Glide gscore over the UNION of all AF3 and Boltz samples.

Submission 19 compared only two structures per ligand, and they were not chosen the
same way: 06's was the best PL-iPTM AF3 sample, 18's the best Glide Boltz sample. This
applies one rule to everything - the most negative gscore across all 94 scored AF3
samples and all 100 Boltz samples.

It changes 10 of 20 structures relative to submission 19 and flips the source balance
to 12 AF3 / 8 Boltz (19 was 6/14), because the union sees AF3 samples that submission
06 never shipped - for six ligands a different AF3 sample outscores everything, e.g.
OCNT-2312554 seed-1_sample-1 at -7.01.

Rationale: submission 19 (0.4331 / 3.6754 / 0.7823) showed Glide resolves differences
between whole models well, even though it picks within one ensemble at chance
(analysis 23). The scale explanation is that candidates ~4 A apart are separable while
candidates 1-2 A apart are not. The union pool spans both regimes, so this tests
whether one global rule holds up.

RISK: submission 13 applied exactly this rule within AF3 alone and scored 0.3730, the
worst AF3-based attempt. Twelve of the twenty picks here are AF3 samples chosen on
gscore, which is that same losing rule - so this could easily land below 19.

No new compute: all 194 samples were already PrepWizard-minimised and Glide-scored
(submission 04 glide_rescoring, analysis 23). Hydrogens are stripped on write.
"""

import csv, os, sys, zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "submission")
ZIP = os.path.join(HERE, "submission_union_glide.zip")
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
    # this submission's selection.csv records model + sample, not a "source" column
    sel = {r["ligand"]: f'{r["model"]}/{r["sample"]}' for r in csv.DictReader(
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
