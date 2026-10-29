"""Score analysis 33's ICM poses with Glide `inplace` - the pose as it stands.

This is the correction to analysis 25. That analysis ranked ICM poses with
`DOCKING_METHOD mininplace`, which minimises the pose in the Glide field before
scoring, so the number described a pose Glide had created rather than the ICM pose
being chosen between. Submission 22 then shipped the raw poses out of `hits.sdf` while
ranking them on those minimised scores: the object scored and the object shipped were
not the same object. Whatever the -0.0171 loss measured, it was not the ranking of the
structures that went in the zip.

`inplace` removes that gap. The score attaches to the exact coordinates that would be
packaged. It also preserves the spread being selected over - minimisation pulls poses
sitting in one basin toward a common minimum, which flattens the differences the
selection depends on.

The known cost, stated so it is not discovered later as a surprise: ICM poses sit in
ICM's own minima and can read as mildly strained under OPLS_2005, so `inplace` may
charge a pose for a clash that is an artefact of changing force field rather than a
defect of the pose. That is the real argument for minimising before cross-program
rescoring. It is accepted here because the alternative scores something we would not
ship.

No quality claim is made by any number this produces. gscore is a scoring function, not
ground truth; the leaderboard is the only measurement of quality this project has. What
these scores are for is ranking candidates within one receptor, which is the single use
of Glide that has a leaderboard result behind it (submission 29, +0.0111).

Grids are reused from analysis 28 rather than rebuilt: every one of submission 29's 20
picks already has a PrepWizard-prepped receptor and grid there, from the identical
protocol. Nothing is re-prepped, so no receptor differs between this and the scores
submission 29 was built from.

CAVEAT on the receptor: ICM docked into the submission PDB as shipped, while the Glide
grid comes from the PrepWizard-prepped version of that same structure. Heavy atoms are
the same - PrepWizard was run at -rmsd 0.3 - but hydrogens and some side-chain flips
differ. The poses are therefore scored in a very slightly different receptor from the
one they were docked in. This is unavoidable without re-prepping, and it is the same
mismatch analysis 25 carried.
"""
import argparse, csv, glob, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "work")
OUT_CSV = os.path.join(HERE, "results", "inplace_on_icm_poses.csv")
SEL = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                   "results", "selection.csv")
A28 = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "work")
A33 = os.path.join(ROOT, "analyses", "33_icm_redock_sub29", "structures")
SCHRODINGER = os.environ.get("SCHRODINGER", "/opt/applications/schrodinger/2025-3")
GLIDE = os.path.join(SCHRODINGER, "glide")

SCORE_IN = """GRIDFILE   {grid}
LIGANDFILE   poses.mae
DOCKING_METHOD   inplace
FORCEFIELD   OPLS_2005
POSTDOCK   True
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
POSES_PER_LIG   1
"""
TERMS = ["gscore", "energy", "evdw", "ecoul", "lipo", "hbond", "rewards",
         "erotb", "esite", "metal"]
COLS = ["ligand", "picked_set", "pose_rank", "status"] + TERMS


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def one(lig, pset, structure, pd):
    d = os.path.join(WORK, lig)
    os.makedirs(d, exist_ok=True)
    grid = glob.glob(os.path.join(A28, pset, lig, "grid*.zip"))
    if not grid:
        raise RuntimeError(f"no analysis 28 grid for {pset}/{lig}")
    hits = os.path.join(A33, f"sub29__{lig}", "hits.sdf")
    if not os.path.exists(hits):
        raise RuntimeError("no hits.sdf from analysis 33")

    # Retitle each pose so the Glide output can be matched back to its ICM rank:
    # ICM writes the same title for all ten.
    posef = os.path.join(d, "poses.mae")
    ranks = []
    if not os.path.exists(posef):
        w = structure.StructureWriter(posef)
        for i, st in enumerate(structure.StructureReader(hits)):
            st.title = f"{lig}__pose{i + 1:02d}"
            w.append(st)
        w.close()
    for i, st in enumerate(structure.StructureReader(posef)):
        ranks.append(st.title)

    csvp = os.path.join(d, "score.csv")
    if not os.path.exists(csvp):
        with open(os.path.join(d, "score.in"), "w") as f:
            f.write(SCORE_IN.format(grid=grid[0]))
        subprocess.run([GLIDE, "-NOJOBID", "-WAIT", "score.in"], cwd=d, check=False,
                       capture_output=True, text=True, timeout=7200)
    if not os.path.exists(csvp):
        raise RuntimeError("glide produced no results csv")

    g = pd.read_csv(csvp)
    title_col = next((c for c in ("title", "SMILES", "s_m_title") if c in g.columns),
                     g.columns[0])
    rows = []
    for t in ranks:
        rank = int(t.split("pose")[-1])
        hit = g[g[title_col].astype(str) == t]
        row = {"ligand": lig, "picked_set": pset, "pose_rank": rank}
        if hit.empty:
            row["status"] = "no score returned"
        else:
            r0 = hit.iloc[0]
            for term in TERMS:
                c = f"r_i_glide_{term}"
                if c in g.columns:
                    try:
                        row[term] = float(r0[c])
                    except (TypeError, ValueError):
                        pass
            row["status"] = "ok"
        rows.append(row)
    return rows


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    from schrodinger import structure
    import pandas as pd

    picks = [(r["ligand"], r["picked_set"]) for r in csv.DictReader(open(SEL))]
    if a.only:
        picks = [p for p in picks if p[0] == a.only]
    done = set()
    if os.path.exists(OUT_CSV):
        done = {r["ligand"] for r in csv.DictReader(open(OUT_CSV))}

    print(f"{len(picks)} ligands, {len(done)} already scored", flush=True)
    for k, (lig, pset) in enumerate(picks, 1):
        if lig in done:
            print(f"[{k}/{len(picks)}] {lig}  already done", flush=True)
            continue
        try:
            rows = one(lig, pset, structure, pd)
            got = [r for r in rows if r.get("gscore") is not None]
            msg = (f"{len(got)}/{len(rows)} scored, best "
                   f"{min((r['gscore'] for r in got), default=float('nan')):.3f}")
        except Exception as e:
            rows = [{"ligand": lig, "picked_set": pset, "pose_rank": 0,
                     "status": f"error: {e}"[:200]}]
            msg = f"ERROR: {e}"
        new = not os.path.exists(OUT_CSV)
        with open(OUT_CSV, "a", newline="") as fh:
            w = csv.DictWriter(fh, COLS, extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerows(rows)
        print(f"[{k}/{len(picks)}] {lig} ({pset})  {msg}", flush=True)
    print("\ndone")
    return 0


if __name__ == "__main__":
    sys.exit(main())
