"""Analysis 30 - Glide `inplace` scores for the AF3 and Boltz-2 structures.

Submission 29 selected among three ESMFold2 runs by `inplace` gscore and scored
0.4638 / 3.7277 / 0.8001 - the project's best LDDT-PLI and the first submission to beat
its best parent on all three metrics. Every Glide selection that FAILED here
(submissions 13, 18, 20, 24) had used `mininplace`, so the mode looks like a confound we
introduced ourselves rather than evidence that physics scores cannot select.

Submission 24 is the sharpest test of that. It was the three-model hybrid - per-ligand
best `mininplace` gscore across AF3 (submission 06), Boltz-2 (submission 18) and
ESMFold2 (submission 23) - and it scored 0.4348, which is 0.0179 BELOW simply shipping
ESMFold2 unaltered. Rerunning that selection under `inplace` asks whether the three-model
cross-model regime also flips, or whether submission 29 was a one-off.

ESMFold2's `inplace` scores already exist, from analysis 28. This supplies the other two:
the 20 shipped structures of submission 06 and the 20 of submission 18, scored `inplace`
only. Those are exactly the structures submission 19 and submission 24 chose between, so
the comparison is like-for-like.

One ligand, OCNT-2395789, has no AF3 structure in submission 06 - submission 19 noted it
and took Boltz there - so that set contributes 19, not 20.

Prep is run fresh here rather than reusing earlier work directories, so all sets reach
their grid the same way. PrepWizard writes .mae: a PDB carries no bond orders and Glide's
Lewis perception then fails on the TYR and ARG side chains (analysis 26).
"""
import argparse, csv, glob, os, shutil, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "work")
OUT_CSV = os.path.join(HERE, "results", "inplace_scores.csv")
SCHRODINGER = os.environ.get("SCHRODINGER", "/opt/applications/schrodinger/2025-3")
GLIDE = os.path.join(SCHRODINGER, "glide")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")

SETS = {"06_af3_pliptm": "06_af3_pliptm",
        "18_boltz2_glide": "18_boltz2_glide"}
MODES = ["inplace"]

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
