"""
Step 1: Redock each crystal ligand into its own receptor with Glide SP confgen.

Reuses the prepped receptors and grids already built by analysis 01 — nothing is
re-prepared here, so this only runs Glide.

Unlike analysis 01 (which scored the crystal pose in place), this generates poses from
scratch and keeps the top N, so we can separate two different questions:

    top-1  RMSD : did the scoring function *pick* the right pose?
    best-of-N RMSD : did sampling *find* it at all?

Usage (via scripts/run_redock.sh):
    $SCHRODINGER/run python3 scripts/01_redock.py [--only 2V0M,4K9T] [--poses 10]

Outputs (git-ignored): ../01_glide_score_vs_ic50/structures/<PDB>/<PDB>_redock*
Outputs (committed):   results/redock_status.csv
"""

import argparse
import os
import subprocess

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
STRUCT_DIR = os.path.join(A01, "structures")
RESULTS_DIR = os.path.join(HERE, "results")

SCHRODINGER = os.environ["SCHRODINGER"]
GLIDE = os.path.join(SCHRODINGER, "glide")

GLIDE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   confgen
PRECISION   SP
FORCEFIELD   OPLS_2005
POSTDOCK   True
POSES_PER_LIG   {poses}
POSTDOCK_NPOSE   {poses}
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
CALC_INPUT_RMS   True
"""


def require_compute_node(what):
    """Refuse to launch CPU-heavy Schrodinger tools outside a SLURM allocation.

    PrepWizard, grid generation and Glide must not run on the shared login node.
    Set ALLOW_LOGIN_NODE=1 to override (don't).
    """
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(
        f"REFUSING to run {what} outside SLURM - this would load the login node.\n"
        f"Submit it instead, e.g.:  sbatch scripts/run_*.sh"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated PDB IDs")
    ap.add_argument("--poses", type=int, default=10, help="poses to keep per ligand")
    args = ap.parse_args()

    cur = pd.read_csv(os.path.join(A01, "data", "curated_set.csv"))
    if args.only:
        cur = cur[cur["pdb_id"].isin(args.only.split(","))]
    os.makedirs(RESULTS_DIR, exist_ok=True)

    rows = []
    for _, row in cur.iterrows():
        pdb = row["pdb_id"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        grid = os.path.join(wdir, f"{pdb}_grid-gridgen.zip")
        lig = os.path.join(wdir, f"{pdb}_lig.sdf")
        job = f"{pdb}_redock"
        print(f"\n== {pdb} ({row['comp_id']})")
        out = {"pdb_id": pdb, "comp_id": row["comp_id"]}
        try:
            for path, what in ((grid, "grid"), (lig, "ligand")):
                if not os.path.exists(path):
                    raise RuntimeError(f"missing {what}: run analysis 01 first")

            in_file = os.path.join(wdir, f"{job}.in")
            with open(in_file, "w") as f:
                f.write(GLIDE_IN.format(grid=os.path.basename(grid),
                                        lig=os.path.basename(lig), poses=args.poses))
            require_compute_node("Glide")
            cmd = [GLIDE, "-NOJOBID", "-WAIT", os.path.basename(in_file)]
            print("    $", " ".join(cmd))
            subprocess.run(cmd, cwd=wdir, check=False)

            csv_path = os.path.join(wdir, f"{job}.csv")
            if not os.path.exists(csv_path):
                raise RuntimeError("no Glide output CSV")
            df = pd.read_csv(csv_path)
            out["n_poses"] = len(df)
            out["status"] = "ok" if len(df) else "no poses generated"
            if len(df):
                out["best_score"] = float(df["r_i_docking_score"].min())
            print(f"    {out['n_poses']} poses")
        except Exception as e:
            print(f"    ERROR: {e}")
            out["status"] = f"error: {e}"
        rows.append(out)

    df = pd.DataFrame(rows)
    out_csv = os.path.join(RESULTS_DIR, "redock_status.csv")
    if args.only and os.path.exists(out_csv):
        old = pd.read_csv(out_csv)
        df = pd.concat([old[~old["pdb_id"].isin(df["pdb_id"])], df], ignore_index=True)
    df.to_csv(out_csv, index=False)
    print(f"\n{(df['status'] == 'ok').sum()}/{len(df)} docked. Wrote {out_csv}")


if __name__ == "__main__":
    main()
