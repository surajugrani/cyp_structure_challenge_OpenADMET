"""Can Glide mininplace flatten the ring PrepWizard could not?

OCNT-2313038's ESMFold2 pose is the one PoseBusters failure in submission 23
(non-aromatic_ring_non-flatness). PrepWizard could not fix it: at -rmsd 1.0 the ligand
moved only 0.19 A mean / 0.32 A max and the ring stayed flagged, so the restraint was
never the limiting factor. Submission 23 therefore shipped submission 19's structure for
this ligand.

Glide `mininplace` minimises the ligand in the receptor field, which is a different
operation from PrepWizard's restrained all-system minimisation, so it may relax the ring.
This scores the pose AND extracts the minimised ligand, grafts it back, and re-runs
PoseBusters to find out.

If it passes, submission 23 could be made purely ESMFold2 on all 20.
"""
import csv, glob, os, shutil, subprocess, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S23 = os.path.join(ROOT, "submissions", "23_esmfold2_default")
WORK = os.path.join(HERE, "work_2313038")
LIG = "OCNT-2313038"
SCHRODINGER = os.environ.get("SCHRODINGER", "")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))

SCORE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   mininplace
FORCEFIELD   OPLS_2005
POSTDOCK   True
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
POSES_PER_LIG   1
"""


def main():
    if not (os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE")):
        raise SystemExit("compute node only")
    os.makedirs(WORK, exist_ok=True)
    from schrodinger import structure
    from rdkit import Chem, RDLogger
    from posebusters_utils import assign_bond_orders_from_template
    from posebusters import PoseBusters
    RDLogger.DisableLog("rdApp.*")

    # the ESMFold2 pose for this ligand, PrepWizard-minimised, before the swap
    src = os.path.join(S23, "work", f"{LIG}_pw.pdb")
    if not os.path.exists(src):
        raise SystemExit(f"missing {src}")
    local = os.path.join(WORK, f"{LIG}.pdb")
    shutil.copy(src, local)

    st = next(iter(structure.StructureReader(local)))
    lig_idx = [a.index for a in st.atom if a.pdbres.strip() == "LIG"]
    lig_set, nfe = set(lig_idx), 0
    for a in st.atom:
        if a.element != "Fe":
            continue
        for b in list(a.bonded_atoms):
            if b.index in lig_set:
                st.deleteBond(a, b); nfe += 1
    st.extract(lig_idx).write(os.path.join(WORK, "lig.mae"))
    st.write(os.path.join(WORK, "complex.mae"))
    print(f"broke {nfe} Fe-ligand bonds")

    grid = glob.glob(os.path.join(WORK, "*grid*.zip"))
    if not grid:
        subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT", "-rec_file", "complex.mae",
                        "-lig_asl", "res. LIG", "-j", "grid"],
                       cwd=WORK, check=False, capture_output=True, text=True,
                       timeout=3600)
        grid = glob.glob(os.path.join(WORK, "*grid*.zip"))
    if not grid:
        raise SystemExit("no grid produced")

    pv = os.path.join(WORK, "score_pv.mae")
    if not os.path.exists(pv):
        with open(os.path.join(WORK, "score.in"), "w") as f:
            f.write(SCORE_IN.format(grid=os.path.basename(grid[0]), lig="lig.mae"))
        subprocess.run([GLIDE, "-NOJOBID", "-WAIT", "score.in"], cwd=WORK,
                       check=False, capture_output=True, text=True, timeout=3600)
    if not os.path.exists(pv):
        raise SystemExit("no poseviewer output - mininplace failed")

    # poseviewer: entry 0 is the receptor, entry 1 the minimised ligand
    poses = list(structure.StructureReader(pv))
    mini = poses[1]
    print(f"minimised ligand: {mini.atom_total} atoms")

    # graft the minimised coordinates back onto the shipped-format PDB
    base = [l for l in open(local)]
    coords = [a.xyz for a in mini.atom if a.element != "H"]
    out, k, moved = [], 0, []
    for l in base:
        if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "LIG" \
                and (l[76:78].strip() or l[12:16].strip()[0]) != "H":
            old = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            x, y, z = coords[k]; k += 1
            moved.append(sum((a - b) ** 2 for a, b in zip(old, (x, y, z))) ** 0.5)
            l = f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}"
        out.append(l)
    dst = os.path.join(WORK, f"{LIG}_mininplace.pdb")
    open(dst, "w").writelines(out)
    print(f"grafted {k} ligand atoms; displacement mean "
          f"{sum(moved)/len(moved):.2f} A, max {max(moved):.2f} A")

    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = PoseBusters(config="dock")
    pred = Chem.MolFromPDBFile(dst, removeHs=False, sanitize=False)
    frags = Chem.SplitMolByPDBResidues(pred)
    lg = frags.pop("LIG")
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    tpl = Chem.MolFromSmiles(smiles[LIG])
    try:
        lg = assign_bond_orders_from_template(tpl, lg)
    except Exception:
        pass
    r = pb.bust(mol_pred=lg, mol_true=tpl, mol_cond=prot)
    bc = [c for c in r.columns if r[c].dtype == bool]
    failed = sorted(r[bc].iloc[0].index[~r[bc].iloc[0]].tolist())
    print(f"\nPoseBusters after mininplace: {'PASS' if not failed else failed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
