"""Submission 22: submission 19, with 5 structures replaced by Glide-ranked ICM poses.

For each ligand, if the ICM redocked pose that GLIDE ranks best sits closer to the haem
iron than submission 19's structure, use it. Five qualify:

    OCNT-2395712   3.72 -> 2.29 A   crosses into coordination (has an N donor)
    OCNT-2395460   3.71 -> 2.19 A   crosses into coordination BUT HAS NO N DONOR
    OCNT-2535422   4.02 -> 3.56 A
    OCNT-2315379   4.05 -> 3.56 A
    OCNT-2312966   3.59 -> 3.41 A

Net effect on chemistry: coordination among the 15 Fe-capable ligands goes 11 -> 12,
median Fe distance 2.29 -> 2.21 A.

THE EVIDENCE IS AGAINST THIS HELPING.

  * Submission 11 ran almost exactly this logic - swap three structures for ones closer
    to the iron - and LOST 0.005 LDDT-PLI. "Closer to Fe" was not a proxy for a better
    pose, because both models already place ligands in the right pocket (analysis 16) and
    the error is orientation within it.
  * Three of the five swaps are small gains on ligands that remain non-coordinating
    (3.59->3.41, 4.05->3.56, 4.02->3.56), so they buy almost nothing while importing an
    ICM pose that analysis 25 measured at a median 4.75 A from the co-folded structure.
  * OCNT-2395460 has NO nitrogen donor, and analysis 19 found 0 of 7 such crystal ligands
    coordinate the iron. Its new 2.19 A contact is chemically suspect, not a gain.

The one defensible swap is OCNT-2395712, which gains a real coordination on a ligand that
can form one. Run at the user's request to test the heuristic directly rather than argue
it.

The ICM poses come from analysis 17 (docked into submission 06's receptors, ligand built
from SMILES) and were ranked by Glide mininplace in analysis 25.
"""

import argparse, csv, os, sys
import zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S06 = os.path.join(ROOT, "submissions", "06_af3_pliptm", "submission")
A17 = os.path.join(ROOT, "analyses", "17_icm_redock_predicted")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_icm_fe_swap.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")
EL2NAME = defaultdict(int)


def glide_pose(name):
    """Best-Score pose from the af3_06 hitlist, as an RDKit mol."""
    import csv as _csv
    best = None
    for r in _csv.DictReader(open(os.path.join(
            ROOT, "analyses", "25_glide_on_icm_poses", "results",
            "pick_comparison.csv"))):
        if r["ligand"] == name:
            best = {"pose_rank": r["glide_rank"], "Score": r["glide_rmsd"],
                    "fe": r["glide_fe"]}
    if best is None:
        raise RuntimeError(f"{name}: no Glide ranking for this ligand")
    hits = os.path.join(A17, "structures", f"af3_06__{name}", "hits.sdf")
    # ICM writes the hitlist worst-first: pose 10 of 10 carries the best Score.
    mols = [m for m in Chem.SDMolSupplier(hits, removeHs=True, sanitize=False)]
    idx = int(best["pose_rank"]) - 1
    if idx >= len(mols) or mols[idx] is None:
        raise RuntimeError(f"{name}: pose {idx+1} unreadable in {hits}")
    return mols[idx], best


def write_complex(name, mol, dst):
    """Protein + HEM from 06, then the ICM ligand as LIG with fresh names and serials."""
    keep = [l for l in open(os.path.join(S06, f"{name}.pdb"))
            if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() != "LIG"
            and (l[76:78].strip() or l[12:16].strip()[0]) != "H"]
    serial = max(int(l[6:11]) for l in keep) + 1
    chain = "C" if "C" not in {l[21] for l in keep} else "L"
    resseq = 1
    conf = mol.GetConformer()
    counts = defaultdict(int)
    lines = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            continue
        el = a.GetSymbol().upper()
        counts[el] += 1
        nm = f"{el}{counts[el]}"[:4]
        p = conf.GetAtomPosition(a.GetIdx())
        # PDB columns: 13-16 atom name, 17 altLoc, 18-20 resName, 22 chainID.
        # An earlier version omitted the altLoc space, pushing resName to 17-19 so
        # l[17:20] read "IG " and every structure failed with "no LIG residue".
        lines.append(f"HETATM{serial:5d} {nm:<4s} LIG {chain}{resseq:4d}    "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          "
                     f"{a.GetSymbol().rjust(2)}\n")
        serial += 1
    open(dst, "w").writelines(keep + lines + ["END\n"])
    return len(lines)


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


SWAP = ["OCNT-2312966", "OCNT-2315379", "OCNT-2395460",
        "OCNT-2395712", "OCNT-2535422"]
S19 = os.path.join(ROOT, "submissions", "19_hybrid_glide", "submission")


def main():
    """Copy submission 19, rebuilding only the 5 ligands whose Glide-ranked ICM pose
    sits closer to the haem iron."""
    import shutil
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    tmp = os.path.join(RES, "_tmp")
    os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    rows = []
    for n in sorted(smiles):
        if n in SWAP:
            mol, rec = glide_pose(n)
            na = write_complex(n, mol, os.path.join(OUT, f"{n}.pdb"))
            rows.append({"ligand": n, "source": "icm_glide_ranked",
                         "pose_rank": rec["pose_rank"], "fe_dist": rec["fe"],
                         "n_lig_atoms": na})
            print(f"  {n}: ICM pose rank {rec['pose_rank']}, Fe {rec['fe']} A"
                  f"   <-- SWAPPED")
        else:
            shutil.copy(os.path.join(S19, f"{n}.pdb"),
                        os.path.join(OUT, f"{n}.pdb"))
            rows.append({"ligand": n, "source": "sub_19", "pose_rank": "",
                         "fe_dist": "", "n_lig_atoms": ""})
    nsw = sum(r["source"] == "icm_glide_ranked" for r in rows)
    print(f"\n{nsw}/20 rebuilt from Glide-ranked ICM poses, "
          f"{20 - nsw}/20 copied from submission 19")

    pb = PoseBusters(config="dock")
    ns = npb = 0
    for r in rows:
        p = os.path.join(OUT, f"{r['ligand']}.pdb")
        se = struct_errors(open(p).readlines(), r["ligand"],
                           smiles[r["ligand"]], tmp)
        pf = pb_fail(p, smiles[r["ligand"]], pb)
        r["struct_errors"] = "; ".join(se)
        r["pb_failed"] = "; ".join(pf)
        ns += not se
        npb += not pf
        if se or pf:
            print(f"  FINAL {r['ligand']}: struct {se} | PB {pf}")
    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nFINAL structural {ns}/20   PoseBusters {npb}/20")
    if ns == 20 and npb == 20:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for r in rows:
                z.write(os.path.join(OUT, f"{r['ligand']}.pdb"),
                        f"{r['ligand']}.pdb")
        print(f"wrote {ZIP}")
    else:
        print("zip not written - failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
