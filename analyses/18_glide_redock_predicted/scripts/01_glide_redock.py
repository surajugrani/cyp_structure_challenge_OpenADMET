"""Analysis 18 - redock all 20 blinded ligands into AF3's receptor with Glide.

The Glide counterpart to analysis 17, so the two docking programs can be compared on the
same receptors and the same prospective footing. Analysis 17 found ICM top-1 sits a
median 5.04 A from the co-folded pose, keeps Fe coordination for only 4 of the 12 ligands
chemically able to form it (AF3 gets 11/12), and ejects two ligands from the pocket
entirely (OCNT-2313429 at 10.23 A, OCNT-2315379 at 8.62 A).

Receptor: the PrepWizard-minimised AF3 sample that submission 06 ships, taken from
submission 04's glide_rescoring work directory so the hydrogens PrepWizard added are
still present (06's shipped PDBs are heavy-atom only, and Glide wants the H).

Prospective by construction: the ligand is built from SMILES with LigPrep, never from the
co-folded coordinates. The co-folded ligand defines the grid box and nothing else.

Settings chosen to match analysis 17's ICM run as closely as the two programs allow:
POSES_PER_LIG 10 against ICM's confs=10, SP precision, flexible ligand docking.

Reuses the grid recipe proven in submission 04:
  * command-line -rec_file / -lig_asl / -j, NOT a keyword .in file. In schrodinger/2025-3
    an .in file is parsed as CSV and fails with "Input file header is missing 'rec_file'
    column".
  * output is <job>-gridgen.zip, not <job>.zip.
  * Fe-ligand bonds are broken before the ligand is removed. PrepWizard (-rehtreat)
    bonds a coordinating ligand atom to the haem iron; leaving that bond in place would
    hand the grid a dangling valence.

Two complexes at a time (--workers, default 2) in a single job - the licence allows two
concurrent Glide runs. Each complex still runs its own grid, LigPrep and docking
end-to-end; only the outer loop is parallel. Resumable. Compute node only.
"""
import argparse, csv, glob, os, shutil, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
PWORK = os.path.join(ROOT, "submissions", "04_af3_default", "glide_rescoring", "work")
SEL06 = os.path.join(ROOT, "submissions", "06_af3_pliptm", "results", "selection.csv")
WORK = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")
OUT_CSV = os.path.join(RESULTS, "glide_poses.csv")

SCHRODINGER = os.environ.get("SCHRODINGER", "")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
LIGPREP = os.path.join(SCHRODINGER, "ligprep")
GLIDE = os.path.join(SCHRODINGER, "glide")

DOCK_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
PRECISION   SP
POSES_PER_LIG   10
POSTDOCK_NPOSE   10
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
FORCEFIELD   OPLS_2005
"""
TERMS = ["gscore", "emodel", "energy", "evdw", "ecoul", "einternal",
         "lipo", "hbond", "rewards", "erotb", "esite", "metal"]


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def graph_arrays(st):
    idx = [a.index for a in st.atom if a.element != "H"]
    pos = {j: i for i, j in enumerate(idx)}
    xyz = np.array([st.atom[j].xyz for j in idx])
    z = np.array([st.atom[j].atomic_number for j in idx])
    adj = np.zeros((len(idx), len(idx)), dtype=int)
    for b in st.bond:
        i, j = b.atom1.index, b.atom2.index
        if i in pos and j in pos:
            adj[pos[i], pos[j]] = adj[pos[j], pos[i]] = 1
    return xyz, z, adj


def process(n, smiles, chosen, structure, spy):
    """Grid, LigPrep, dock and score one complex. Returns the rows to append."""
    d = os.path.join(WORK, n); os.makedirs(d, exist_ok=True)
    tag = chosen[n]
    src = os.path.join(PWORK, f"{n}__{tag}", f"{n}__{tag}_pw.pdb")
    if not os.path.exists(src):
        raise RuntimeError(f"missing prepped receptor {src}")
    local = os.path.join(d, "complex_pw.pdb")
    if not os.path.exists(local):
        shutil.copy(src, local)

    st = next(iter(structure.StructureReader(local)))
    lig_idx = [a.index for a in st.atom if a.pdbres.strip() == "LIG"]
    if not lig_idx:
        raise RuntimeError("no LIG residue")
    lig_set = set(lig_idx)
    nfe = 0
    for a in st.atom:
        if a.element != "Fe":
            continue
        for b in list(a.bonded_atoms):
            if b.index in lig_set:
                st.deleteBond(a, b); nfe += 1

    cof = os.path.join(d, "cofolded_lig.mae")
    if not os.path.exists(cof):
        st.extract(lig_idx).write(cof)
    cmae = os.path.join(d, "complex.mae")
    if not os.path.exists(cmae):
        st.write(cmae)

    grid = glob.glob(os.path.join(d, "grid*.zip"))
    if not grid:
        subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT", "-rec_file", "complex.mae",
                        "-lig_asl", "res. LIG", "-j", "grid"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=3600)
        grid = glob.glob(os.path.join(d, "grid*.zip"))
    if not grid:
        raise RuntimeError("no grid produced")

    ligmae = os.path.join(d, "ligprep-out.maegz")
    if not os.path.exists(ligmae):
        with open(os.path.join(d, "lig.smi"), "w") as f:
            f.write(f"{smiles[n]} {n}\n")
        subprocess.run([LIGPREP, "-NOJOBID", "-WAIT", "-ismi", "lig.smi",
                        "-omae", "ligprep-out.maegz", "-epik", "-ph", "7.4", "-s", "1"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=3600)
    if not os.path.exists(ligmae):
        raise RuntimeError("ligprep produced no output")

    pv = os.path.join(d, "dock_pv.mae")
    if not os.path.exists(pv):
        with open(os.path.join(d, "dock.in"), "w") as f:
            f.write(DOCK_IN.format(grid=os.path.basename(grid[0]),
                                   lig="ligprep-out.maegz"))
        subprocess.run([GLIDE, "-NOJOBID", "-WAIT", "dock.in"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=7200)
    if not os.path.exists(pv):
        raise RuntimeError("glide produced no poseviewer output")

    cx, cz, ca = graph_arrays(next(iter(structure.StructureReader(cof))))
    fe = None
    for a in next(iter(structure.StructureReader(cmae))).atom:
        if a.element == "Fe":
            fe = np.array(a.xyz); break
    rows = []
    for k, pose in enumerate(structure.StructureReader(pv)):
        if k == 0:                      # poseviewer entry 0 is the receptor
            continue
        px, pz, pa = graph_arrays(pose)
        try:
            r = float(spy.symmrmsd(cx, px, cz, pz, ca, pa, minimize=False))
        except Exception:
            r = np.nan
        row = {"ligand": n, "pose_rank": k, "fe_lig_bonds_broken": nfe,
               "rmsd_to_cofolded": None if np.isnan(r) else round(r, 3),
               "fe_dist": None if fe is None
               else round(float(np.min(np.linalg.norm(px - fe, axis=1))), 2),
               "status": "ok"}
        for t in TERMS:
            v = pose.property.get(f"r_i_glide_{t}")
            if v is not None:
                row[t] = v
        rows.append(row)
    if not rows:
        raise RuntimeError("no poses in poseviewer file")
    return rows


def main():
    require_compute_node()
    os.makedirs(WORK, exist_ok=True); os.makedirs(RESULTS, exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2,
                    help="concurrent Glide runs; 2 is what the licence allows")
    args = ap.parse_args()
    from schrodinger import structure
    from spyrmsd import rmsd as spy

    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    chosen = {r["ligand"]: r["chosen"] for r in csv.DictReader(open(SEL06))}

    done = set()
    if os.path.exists(OUT_CSV):
        done = {r["ligand"] for r in csv.DictReader(open(OUT_CSV))}

    names = [n for n in sorted(smiles) if n not in done]
    print(f"{len(smiles)} ligands, {len(done)} already done, {len(names)} to run "
          f"on {args.workers} workers", flush=True)

    cols = ["ligand", "pose_rank", "status", "fe_lig_bonds_broken",
            "rmsd_to_cofolded", "fe_dist"] + TERMS
    lock = threading.Lock()
    counter = {"k": 0}

    def write(rows):
        with lock:
            new_file = not os.path.exists(OUT_CSV)
            with open(OUT_CSV, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
                if new_file:
                    w.writeheader()
                w.writerows(rows)

    def work(n):
        try:
            rows = process(n, smiles, chosen, structure, spy)
            b = rows[0]
            msg = (f"top-1 gscore {b.get('gscore')}  "
                   f"rmsd_to_cofolded {b['rmsd_to_cofolded']} A  Fe {b['fe_dist']} A")
        except Exception as e:
            rows = [{"ligand": n, "pose_rank": 0, "status": f"error: {e}"[:200]}]
            msg = f"ERROR: {e}"
        write(rows)
        with lock:
            counter["k"] += 1
            print(f"[{counter['k']}/{len(names)}] {n}  {msg}", flush=True)

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        list(ex.map(work, names))

    print("\ndone")
    return 0


if __name__ == "__main__":
    sys.exit(main())
