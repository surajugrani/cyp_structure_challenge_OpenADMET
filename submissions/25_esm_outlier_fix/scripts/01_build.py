"""Submission 25: submission 23 with its two chemically indefensible structures replaced.

Submission 23 (pure ESMFold2) is the best result of the project at 0.4527 LDDT-PLI.
Everything tried since has made it worse: submission 24 let Glide choose across three
co-folding models and landed 0.0179 BELOW simply shipping 23 unaltered. So the premise
here is to change as little as possible.

Two of 23's twenty structures are wrong on chemistry rather than by a small margin, and
both are replaced with submission 19's structure for that ligand:

    OCNT-2315379   6.89 -> 4.05 A   has an N donor, yet sits further from the haem than
                                    any other structure in the set by 2.4 A
    OCNT-2395460   1.79 -> 3.71 A   NO N donor, yet 1.79 A from the iron - too close to
                                    be a real contact, and analysis 19 found 0 of 7 such
                                    crystal ligands coordinate at all

These are not the Fe-distance heuristic that lost in submissions 11 and 22. That rule
moved ligands nearer the iron on the assumption that nearer is better, and it is not:
both models already place ligands in the right pocket (analysis 16) and the error is
orientation. These two swaps instead remove values that are impossible - one far outside
the distribution, one implying a coordination the ligand cannot form - and they happen to
move in OPPOSITE directions, which is what distinguishes them from a distance heuristic.

A third swap was considered and REJECTED: OCNT-2395328, where submission 19 sits at
3.01 A against ESMFold2's 4.03 A. That ligand has no nitrogen donor, so 4.03 A is the
more defensible number and swapping it in would move a non-coordinating ligand closer to
the iron - exactly the move that failed twice. It is built as a variant under
variant_3swap/ so the two can be compared, but it is not the primary submission.

Both inputs are already PrepWizard-minimised and validated 20/20, so structures are
copied byte-for-byte and only re-validated.
"""
import csv
import os
import shutil
import sys
import zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S23 = os.path.join(ROOT, "submissions", "23_esmfold2_default", "submission")
S19 = os.path.join(ROOT, "submissions", "19_hybrid_glide", "submission")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")

SWAP2 = {"OCNT-2315379": "too far: 6.89 A with an N donor",
         "OCNT-2395460": "too close: 1.79 A with no N donor"}
SWAP3 = dict(SWAP2, **{"OCNT-2395328": "user-requested; 4.03 -> 3.01 A, no N donor"})


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


def build(swaps, outdir, zippath, label, smiles, pb, tmp):
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for n in sorted(smiles):
        src = S19 if n in swaps else S23
        shutil.copy(os.path.join(src, f"{n}.pdb"), os.path.join(outdir, f"{n}.pdb"))
        rows.append({"ligand": n,
                     "source": "sub_19" if n in swaps else "sub_23_esmfold2",
                     "reason": swaps.get(n, "")})
    print(f"\n== {label}: {len(swaps)}/20 from submission 19, "
          f"{20 - len(swaps)}/20 unchanged from submission 23")
    ns = npb = 0
    for r in rows:
        p = os.path.join(outdir, f"{r['ligand']}.pdb")
        se = struct_errors(open(p).readlines(), r["ligand"],
                           smiles[r["ligand"]], tmp)
        pf = pb_fail(p, smiles[r["ligand"]], pb)
        r["struct_errors"] = "; ".join(se)
        r["pb_failed"] = "; ".join(pf)
        ns += not se
        npb += not pf
        if se or pf:
            print(f"  FINAL {r['ligand']}: struct {se} | PB {pf}")
    print(f"   structural {ns}/20   PoseBusters {npb}/20")
    if ns == 20 and npb == 20:
        with zipfile.ZipFile(zippath, "w", zipfile.ZIP_DEFLATED) as z:
            for r in rows:
                z.write(os.path.join(outdir, f"{r['ligand']}.pdb"),
                        f"{r['ligand']}.pdb")
        print(f"   wrote {os.path.relpath(zippath, ROOT)}")
    else:
        print("   zip not written - failures above")
    return rows


def main():
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    tmp = os.path.join(HERE, "results", "_tmp")
    os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = PoseBusters(config="dock")

    rows = build(SWAP2, os.path.join(HERE, "submission"),
                 os.path.join(HERE, "submission_esm_outlier_fix.zip"),
                 "PRIMARY, 2 swaps", smiles, pb, tmp)
    vdir = os.path.join(HERE, "variant_3swap")
    build(SWAP3, os.path.join(vdir, "submission"),
          os.path.join(vdir, "submission_esm_outlier_fix_3swap.zip"),
          "VARIANT, 3 swaps incl. OCNT-2395328", smiles, pb, tmp)

    with open(os.path.join(HERE, "results", "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
