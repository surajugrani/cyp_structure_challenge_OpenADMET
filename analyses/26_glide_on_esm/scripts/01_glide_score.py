"""Glide-score the 20 ESMFold2 structures, for use in a cross-model hybrid.

Submission 19 showed Glide can choose between whole co-folding models - it beat both its
parents (0.4331 against 0.4291 and 0.4174) picking per ligand from AF3 and Boltz. It
cannot rank within one model's ensemble (analysis 23: agreement with PL-iPTM at chance;
submission 13: -0.056; submission 20: -0.048), and the explanation is scale - candidates
~4 A apart are separable, 1-2 A apart are not.

ESMFold2 (submission 23) is a genuinely independent third model with the best Fe geometry
of the three (12/15 coordinating, median 2.16 A), so adding it as a third candidate is
the natural extension of the one rule that works.

Input is the shipped submission 23 PDBs, re-prepared here. Two earlier attempts failed
and the reason is worth recording:

  1. Feeding the shipped PDBs straight to generate_glide_grids -> "no grid produced".
  2. Feeding submission 23's own PrepWizard output (work/*_pw.pdb) -> the same error.

Neither was about hydrogens. The gridgen log shows the real cause:

     FATAL mmlewis error: Problems remain with the following atoms

Glide could not assign a Lewis structure to the receptor, and the flagged atoms are TYR
ring/hydroxyl and ARG guanidinium atoms. PDB carries no bond orders, so every read of a
.pdb re-perceives aromatic rings and charged groups from geometry alone, and on the
ESMFold2 backbone that perception fails.

The fix is to never round-trip through PDB: PrepWizard is re-run here writing **.mae**,
which stores bond orders and formal charges explicitly, and the grid is built from that.
Options are identical to submission 23 and to analysis 23, so the prepared geometry is
the same - only the file format differs.

All 20 are scored, including OCNT-2313038. Submission 23 shipped submission 19's
structure for that ligand because its ESMFold2 pose fails PoseBusters on ring flatness,
but the ESMFold2 pose is scored here anyway: mininplace may relax the ring, which would
make a pure-ESMFold2 submission possible.

`mininplace` not `inplace`: `inplace` returns emodel and einternal as 0.00 (submission
02), and emodel is Glide's own pose-ranking term.
"""

import argparse, csv, glob, os, shutil, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SAMPLES = os.path.join(ROOT, "submissions", "23_esmfold2_default", "submission")
SKIP = set()
WORK = os.path.join(HERE, "work")
OUT_CSV = os.path.join(HERE, "results", "glide_scores.csv")

SCHRODINGER = os.environ.get("SCHRODINGER", "")
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")

SCORE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   mininplace
FORCEFIELD   OPLS_2005
POSTDOCK   True
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
POSES_PER_LIG   1
"""
TERMS = ["gscore", "emodel", "energy", "evdw", "ecoul", "einternal",
         "lipo", "hbond", "rewards", "erotb", "esite", "metal"]
COLS = ["ligand", "model", "status", "fe_lig_bonds_broken"] + TERMS


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def one(src, structure, pd):
    base = os.path.basename(src)[:-4]
    lig_name, model = base, "esmfold2"
    d = os.path.join(WORK, base)
    os.makedirs(d, exist_ok=True)
    row = {"ligand": lig_name, "model": model}

    local = os.path.join(d, f"{base}.pdb")
    if not os.path.exists(local):
        shutil.copy(src, local)
    # .mae output, not .pdb: PDB drops bond orders and Glide's Lewis perception then
    # fails on the TYR and ARG side chains (see module docstring).
    prepped = os.path.join(d, f"{base}_pw.mae")
    if not os.path.exists(prepped):
        subprocess.run([PREPWIZARD, "-NOJOBID", "-rehtreat", "-captermini",
                        "-propka_pH", "7.4", "-epik_pH", "7.4", "-f", "OPLS_2005",
                        "-rmsd", "0.3", os.path.basename(local),
                        os.path.basename(prepped)],
                       cwd=d, check=False, capture_output=True, text=True, timeout=3600)
    if not os.path.exists(prepped):
        raise RuntimeError("prepwizard produced no .mae")

    st = next(iter(structure.StructureReader(prepped)))
    lig_idx = [a.index for a in st.atom if a.pdbres.strip() == "LIG"]
    if not lig_idx:
        raise RuntimeError("no LIG residue after prepwizard")
    # Break any covalent Fe-ligand bond before extracting the ligand: PrepWizard
    # (-rehtreat) bonds a coordinating atom to the iron, and Glide cannot score a
    # covalently attached ligand in place. Same fix as submission 04.
    lig_set, nfe = set(lig_idx), 0
    for a in st.atom:
        if a.element != "Fe":
            continue
        for b in list(a.bonded_atoms):
            if b.index in lig_set:
                st.deleteBond(a, b); nfe += 1
    row["fe_lig_bonds_broken"] = nfe

    ligf = os.path.join(d, "lig.mae")
    if not os.path.exists(ligf):
        st.extract(lig_idx).write(ligf)
    cmae = os.path.join(d, "complex.mae")
    if not os.path.exists(cmae):
        st.write(cmae)

    grid = glob.glob(os.path.join(d, "grid*.zip"))
    if not grid:
        # command-line -rec_file/-lig_asl, NOT a keyword .in file: in 2025-3 an .in is
        # parsed as CSV and fails with "missing 'rec_file' column". Output is
        # <job>-gridgen.zip.
        subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT", "-rec_file", "complex.mae",
                        "-lig_asl", "res. LIG", "-j", "grid"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=3600)
        grid = glob.glob(os.path.join(d, "grid*.zip"))
    if not grid:
        raise RuntimeError("no grid produced")

    csvp = os.path.join(d, "score.csv")   # score_pv.mae holds the minimised pose
    if not os.path.exists(csvp):
        with open(os.path.join(d, "score.in"), "w") as f:
            f.write(SCORE_IN.format(grid=os.path.basename(grid[0]), lig="lig.mae"))
        subprocess.run([GLIDE, "-NOJOBID", "-WAIT", "score.in"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=3600)
    if not os.path.exists(csvp):
        raise RuntimeError("no glide results csv")
    g = pd.read_csv(csvp)
    if g.empty:
        raise RuntimeError("glide returned no pose")
    r0 = g.iloc[0]
    for t in TERMS:
        c = f"r_i_glide_{t}"
        if c in g.columns:
            try:
                row[t] = float(r0[c])
            except (TypeError, ValueError):
                pass
    row["status"] = "ok"
    return row


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    from schrodinger import structure
    import pandas as pd

    done = set()
    if os.path.exists(OUT_CSV):
        done = {(r["ligand"], r["model"]) for r in csv.DictReader(open(OUT_CSV))}
    files = [f for f in sorted(glob.glob(os.path.join(SAMPLES, "*.pdb")))
             if os.path.basename(f)[:-4] not in SKIP
             and (os.path.basename(f)[:-4], "esmfold2") not in done]
    print(f"{len(files)} samples to score ({len(done)} already done), "
          f"{a.workers} concurrent", flush=True)

    lock = threading.Lock()
    n = {"k": 0}

    def work(src):
        base = os.path.basename(src)[:-4]
        try:
            row = one(src, structure, pd)
            msg = f"gscore={row.get('gscore')} emodel={row.get('emodel')}"
        except Exception as e:
            row = {"ligand": base, "model": "esmfold2",
                   "status": f"error: {e}"[:200]}
            msg = f"ERROR: {e}"
        with lock:
            new = not os.path.exists(OUT_CSV)
            with open(OUT_CSV, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerow(row)
            n["k"] += 1
            print(f"[{n['k']}/{len(files)}] {base}  {msg}", flush=True)

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(work, files))
    print(f"\nwrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
