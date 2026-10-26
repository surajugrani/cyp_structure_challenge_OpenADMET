"""Score all three ESMFold2 submissions with BOTH Glide in-place modes.

Every Glide-based selection this project has shipped - submissions 13, 18, 19, 20, 21,
24 - used DOCKING_METHOD mininplace. The choice rests on one observation from submission
02: plain `inplace` returns emodel and einternal as 0.00, and emodel is Glide's own
pose-ranking term. It has never been tested against a leaderboard metric. The only
head-to-head anywhere here is analysis 01, which scored 20 crystal poses both ways
against IC50 and found `inplace` the better of the two (Spearman -0.265 against -0.093),
though neither was significant and that analysis concluded no Glide score predicts
potency.

Testing the modes by shipping two submissions would cost two attempts. This costs none.
Score the three ESMFold2 sets both ways, then ask whether the two modes would even pick
different structures.

    23  default endpoint + our ColabFold MSA     0.4527  (best of the project)
    26  max_accuracy endpoint, its own MSA       0.4365
    27  default endpoint, no MSA                 0.4384

If the 20 selections are identical, the resulting submissions are byte-identical and
their scores must be too - submissions 18 and 21 already demonstrated that, scoring the
same to four decimals - so the question closes with no attempt spent. Only if the
selections differ is there anything to submit.

TEMPER EXPECTATIONS ABOUT THE SELECTION ITSELF. These three sets are far closer together
than the models in submission 19: analysis 27 measured a median pocket-superposed ligand
RMSD of 1.53 A between sets 23 and 27, which is the within-ensemble regime where Glide
has failed every time (submissions 13, 18, 20, and 24 at the three-model level). A
submission built from this is not expected to beat 0.4527. The point is the mode
comparison, which is free.

Prep is run fresh and identically for all 60 structures rather than reusing analysis
26's work directory for set 23, so no set is advantaged by a different path to its grid.

PrepWizard writes .mae, not .pdb: a PDB carries no bond orders, and re-reading one makes
Glide re-perceive aromatic rings and charged groups from geometry, which fails on these
structures with "FATAL mmlewis error" on the TYR and ARG side chains (analysis 26).
"""
import argparse, csv, glob, os, shutil, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "work")
OUT_CSV = os.path.join(HERE, "results", "scores_both_modes.csv")
SCHRODINGER = os.environ.get("SCHRODINGER", "/opt/applications/schrodinger/2025-3")
GLIDE = os.path.join(SCHRODINGER, "glide")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")

SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
MODES = ["mininplace", "inplace"]

SCORE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   {mode}
FORCEFIELD   OPLS_2005
POSTDOCK   True
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
POSES_PER_LIG   1
"""
TERMS = ["gscore", "emodel", "energy", "evdw", "ecoul", "einternal",
         "lipo", "hbond", "rewards", "erotb", "esite", "metal"]
COLS = (["ligand", "set", "mode", "status", "fe_lig_bonds_broken"] + TERMS)


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def one(setname, src, structure, pd):
    """Prep, grid, then score the same ligand twice - once per mode."""
    lig_name = os.path.basename(src)[:-4]
    d = os.path.join(WORK, setname, lig_name)
    os.makedirs(d, exist_ok=True)

    local = os.path.join(d, f"{lig_name}.pdb")
    if not os.path.exists(local):
        shutil.copy(src, local)
    prepped = os.path.join(d, f"{lig_name}_pw.mae")
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
    # PrepWizard (-rehtreat) bonds a coordinating atom to the iron; Glide cannot score a
    # covalently attached ligand in place, so the bond is broken first (submission 04).
    lig_set, nfe = set(lig_idx), 0
    for a in st.atom:
        if a.element != "Fe":
            continue
        for b in list(a.bonded_atoms):
            if b.index in lig_set:
                st.deleteBond(a, b); nfe += 1

    ligf = os.path.join(d, "lig.mae")
    if not os.path.exists(ligf):
        st.extract(lig_idx).write(ligf)
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

    rows = []
    for mode in MODES:
        row = {"ligand": lig_name, "set": setname, "mode": mode,
               "fe_lig_bonds_broken": nfe}
        sub = os.path.join(d, mode)
        os.makedirs(sub, exist_ok=True)
        csvp = os.path.join(sub, "score.csv")
        if not os.path.exists(csvp):
            with open(os.path.join(sub, "score.in"), "w") as f:
                f.write(SCORE_IN.format(grid=os.path.join("..",
                                        os.path.basename(grid[0])),
                                        lig=os.path.join("..", "lig.mae"),
                                        mode=mode))
            subprocess.run([GLIDE, "-NOJOBID", "-WAIT", "score.in"],
                           cwd=sub, check=False, capture_output=True, text=True,
                           timeout=3600)
        if not os.path.exists(csvp):
            row["status"] = "error: no glide results csv"
            rows.append(row)
            continue
        g = pd.read_csv(csvp)
        if g.empty:
            row["status"] = "error: glide returned no pose"
            rows.append(row)
            continue
        r0 = g.iloc[0]
        for t in TERMS:
            c = f"r_i_glide_{t}"
            if c in g.columns:
                try:
                    row[t] = float(r0[c])
                except (TypeError, ValueError):
                    pass
        row["status"] = "ok"
        rows.append(row)
    return rows


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
        done = {(r["ligand"], r["set"]) for r in csv.DictReader(open(OUT_CSV))}

    jobs = []
    for setname, subdir in SETS.items():
        sdir = os.path.join(ROOT, "submissions", subdir, "submission")
        for f in sorted(glob.glob(os.path.join(sdir, "*.pdb"))):
            if (os.path.basename(f)[:-4], setname) not in done:
                jobs.append((setname, f))
    print(f"{len(jobs)} structures to score x {len(MODES)} modes "
          f"({len(done)} already done), {a.workers} concurrent", flush=True)

    lock = threading.Lock()
    n = {"k": 0}

    def work(job):
        setname, src = job
        base = os.path.basename(src)[:-4]
        try:
            rows = one(setname, src, structure, pd)
            msg = "  ".join(f"{r['mode']}={r.get('gscore')}" for r in rows)
        except Exception as e:
            rows = [{"ligand": base, "set": setname, "mode": m,
                     "status": f"error: {e}"[:200]} for m in MODES]
            msg = f"ERROR: {e}"
        with lock:
            new = not os.path.exists(OUT_CSV)
            with open(OUT_CSV, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows(rows)
            n["k"] += 1
            print(f"[{n['k']}/{len(jobs)}] {setname} {base}  {msg}", flush=True)

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(work, jobs))
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
