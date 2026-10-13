"""Submission 12: submission 06's receptors, ICM's top-1 redocked ligand poses.

Deliberately testing how much worse ICM redocking is than co-folding. Analysis 17
measured it on the same 20 complexes and the answer looks bad:

  * ICM top-1 disagrees with 06's co-folded pose by a median 5.04 A
  * only 2/20 top-1 poses land within 2 A of the co-folded one
  * ICM keeps the Fe coordination for 4/20 ligands, against 12/20 co-folded

Sampling is not the problem - best-of-10 median is 1.69 A, and 12/20 pools contain a pose
within 2 A of the co-folded one. Ranking is. That is the same split analysis 09 measured
on crystal receptors (ICM top-1 27 % vs a much higher oracle).

Receptor is AF3's (submission 06), not Boltz's, because analysis 17 found AF3's pocket is
the better docking target on every measure despite scoring 0.098 LOWER on LDDT-LP - a
useful dissociation between backbone accuracy and dockability.

Built from the protein + HEM of 06 plus the ICM top-1 pose, with ligand bond orders
restored from SMILES. --allow-pb-fail writes the zip even if PoseBusters rejects some
structures, since the point here is to measure the damage, not to ship a clean set.
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
ZIP = os.path.join(HERE, "submission_icm_redock_af3.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")
EL2NAME = defaultdict(int)


def top1_pose(name):
    """Best-Score pose from the af3_06 hitlist, as an RDKit mol."""
    import csv as _csv
    best = None
    for r in _csv.DictReader(open(os.path.join(A17, "results", "icm_poses.csv"))):
        if r["receptor"] != "af3_06" or r["ligand"] != name or not r.get("Score"):
            continue
        if best is None or float(r["Score"]) < float(best["Score"]):
            best = r
    if best is None:
        raise RuntimeError(f"{name}: no ICM pose")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-pb-fail", action="store_true")
    a = ap.parse_args()
    for d in (OUT, RES):
        os.makedirs(d, exist_ok=True)
    tmp = os.path.join(RES, "_tmp"); os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    rows = []
    for n in sorted(smiles):
        mol, rec = top1_pose(n)
        na = write_complex(n, mol, os.path.join(OUT, f"{n}.pdb"))
        rows.append({"ligand": n, "pose_rank": rec["pose_rank"],
                     "icm_score": rec["Score"], "rtcnn": rec.get("RTCNN", ""),
                     "rmsd_to_cofolded": rec["rmsd_to_cofolded"],
                     "fe_dist": rec["fe_dist"], "n_lig_atoms": na})
        print(f"  {n}: ICM pose {rec['pose_rank']}  Score {rec['Score']}  "
              f"{rec['rmsd_to_cofolded']} A from co-folded  Fe {rec['fe_dist']}")

    pb = PoseBusters(config="dock")
    ns = npb = 0
    for r in rows:
        p = os.path.join(OUT, f"{r['ligand']}.pdb")
        se = struct_errors(open(p).readlines(), r["ligand"], smiles[r["ligand"]], tmp)
        pf = pb_fail(p, smiles[r["ligand"]], pb)
        r["struct_errors"] = "; ".join(se)
        r["pb_failed"] = "; ".join(pf)
        ns += not se
        npb += not pf
        if se or pf:
            print(f"  FINAL {r['ligand']}: struct {se} | PB {pf}")
    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nFINAL structural {ns}/20   PoseBusters {npb}/20")

    if (ns == 20 and npb == 20) or a.allow_pb_fail:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for r in rows:
                z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")
        print(f"wrote {ZIP}" + ("   (PoseBusters failures allowed)"
                                if npb < 20 else ""))
    else:
        print("zip not written - rerun with --allow-pb-fail to ship anyway")
    return 0


if __name__ == "__main__":
    sys.exit(main())
