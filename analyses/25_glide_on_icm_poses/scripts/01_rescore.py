"""Rescore ICM's 10 redocked poses per ligand with Glide, and see if it picks differently.

Analysis 17 found ICM's own Score ships a badly drifted pose for two ligands while much
better ones sit in the same pool:

  OCNT-2313429  8 of 10 poses at 3.70-4.92 A from the iron; ICM shipped rank 10 at
                10.23 A because it carried the best Score (-24.36). Rank 7 is 0.24 A
                from the co-folded pose and was discarded.
  OCNT-2315379  8 of 10 poses at 3.25-4.34 A; ICM shipped rank 10 at 8.62 A
                (Score -22.99). Rank 6 is 0.84 A from the co-folded pose.

So ICM's sampling is fine and its ranking actively prefers poses that leave the catalytic
site. The question here: does Glide rank them differently?

This is the one regime where Glide has shown real skill - submission 19 beat both its
parents choosing between structures ~4 A apart, while within a 1-2 A ensemble it picks at
chance (analysis 23). ICM's 10 poses span up to 10 A, so they sit in the favourable
regime.

Cheap because the grids already exist: submission 04's glide_rescoring built one per AF3
sample, and analysis 17 docked into those same PrepWizard-minimised receptors, so the
coordinate frames match and only the scoring needs running.

`mininplace`, not `inplace`: `inplace` returns emodel and einternal as 0.00 (submission
02), and emodel is Glide's own pose-ranking term.
"""
import argparse, csv, glob, os, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A17 = os.path.join(ROOT, "analyses", "17_icm_redock_predicted")
AF3W = os.path.join(ROOT, "submissions", "04_af3_default", "glide_rescoring", "work")
WORK = os.path.join(HERE, "work")
OUT_CSV = os.path.join(HERE, "results", "glide_on_icm_poses.csv")
OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
SCHRODINGER = os.environ.get("SCHRODINGER", "")
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
COLS = ["ligand", "icm_rank", "status", "icm_score", "icm_rtcnn",
        "fe_dist", "rmsd_to_cofolded"] + TERMS


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    import pandas as pd

    sel = {r["ligand"]: r["chosen"] for r in csv.DictReader(open(os.path.join(
        ROOT, "submissions", "06_af3_pliptm", "results", "selection.csv")))}
    icm = {}
    for r in csv.DictReader(open(os.path.join(A17, "results", "icm_poses.csv"))):
        if r["receptor"] == "af3_06" and r.get("Score"):
            icm[(r["ligand"], int(r["pose_rank"]))] = r

    done = set()
    if os.path.exists(OUT_CSV):
        done = {(r["ligand"], int(r["icm_rank"])) for r in csv.DictReader(open(OUT_CSV))}
    ligs = [l for l in sorted(sel) if (l, 1) not in done]
    print(f"{len(ligs)} ligands to rescore, {a.workers} concurrent", flush=True)

    lock = threading.Lock()

    def work(lg):
        rows = []
        try:
            tag = sel[lg]
            src = os.path.join(AF3W, f"{lg}__{tag}")
            grid = glob.glob(os.path.join(src, "*grid*.zip"))
            if not grid:
                raise RuntimeError("no grid for this receptor")
            hits = os.path.join(A17, "structures", f"af3_06__{lg}", "hits.sdf")
            if not os.path.exists(hits):
                raise RuntimeError("no ICM hitlist")
            d = os.path.join(WORK, lg)
            os.makedirs(d, exist_ok=True)
            # one Glide run per pose keeps the rank mapping unambiguous
            subprocess.run([OBABEL, hits, "-osdf", "-O", "poses.sdf", "-m"],
                           cwd=d, check=False, capture_output=True)
            files = sorted(glob.glob(os.path.join(d, "poses*.sdf")))
            for i, f in enumerate(files, 1):
                jb = f"s{i}"
                csvp = os.path.join(d, f"{jb}.csv")
                if not os.path.exists(csvp):
                    with open(os.path.join(d, f"{jb}.in"), "w") as fh:
                        fh.write(SCORE_IN.format(grid=grid[0],
                                                 lig=os.path.basename(f)))
                    subprocess.run([GLIDE, "-NOJOBID", "-WAIT", f"{jb}.in"],
                                   cwd=d, check=False, capture_output=True,
                                   text=True, timeout=1800)
                row = {"ligand": lg, "icm_rank": i}
                src_icm = icm.get((lg, i), {})
                row["icm_score"] = src_icm.get("Score")
                row["icm_rtcnn"] = src_icm.get("RTCNN")
                row["fe_dist"] = src_icm.get("fe_dist")
                row["rmsd_to_cofolded"] = src_icm.get("rmsd_to_cofolded")
                if os.path.exists(csvp):
                    g = pd.read_csv(csvp)
                    if not g.empty:
                        for t in TERMS:
                            c = f"r_i_glide_{t}"
                            if c in g.columns:
                                try:
                                    row[t] = float(g.iloc[0][c])
                                except (TypeError, ValueError):
                                    pass
                        row["status"] = "ok"
                row.setdefault("status", "no glide output")
                rows.append(row)
            msg = f"{sum(r['status']=='ok' for r in rows)}/{len(rows)} scored"
        except Exception as e:
            rows = [{"ligand": lg, "icm_rank": 0, "status": f"error: {e}"[:200]}]
            msg = f"ERROR: {e}"
        with lock:
            new = not os.path.exists(OUT_CSV)
            with open(OUT_CSV, "a", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerows(rows)
            print(f"  {lg}: {msg}", flush=True)

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(work, ligs))
    print(f"\nwrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
