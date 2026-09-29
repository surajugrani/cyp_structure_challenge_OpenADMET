"""
Analysis 09, step 5: score each structure's CRYSTAL pose with Glide `mininplace`, to give
every ligand an 11th near-native candidate alongside its docked poses - the Glide
counterpart of step 4's `_dockScan -n` rows.

Why `mininplace` and not `inplace`, and why not both:
  - The docked poses being compared against are minimised, so `mininplace` is the matching
    route. Measured on the 7 structures where analysis 01 already had both, `mininplace`
    reproduces the score of a near-identical docked pose to within 0.10-0.18 GlideScore
    units (against a within-ligand spread of ~0.77), while `inplace` is off by 0.29-0.44
    and by +4.92 in one case.
  - `inplace` and `mininplace` are the SAME geometry with and without minimisation. Adding
    both would put two near-duplicate rows per ligand into the pool, double-weighting the
    crystal geometry without adding independent information.

All 30 are re-run here rather than only the 23 that analysis 01 lacks: analysis 01 built
its own grids, and scoring some crystal poses on analysis 01 grids and the rest on
analysis 05 grids would reintroduce exactly the between-route offset this step exists to
avoid. Every row must come from the same grid as the docked poses it will sit beside.

CAVEAT to carry forward: the docked poses were generated from SMILES via LigPrep, so their
protomer/tautomer may differ from the deposited ligand scored here. This is a smaller
effect than the ICM route offset but it is not zero.

Serial, resumable, compute node only.

Usage:
    sbatch analyses/09_pose_selection/scripts/05_glide_crystal.sh

Outputs:
    results/glide_crystal_scores.csv
"""

import glob
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A05_STRUCT = os.path.join(ANALYSES, "05_full_pdb_redocking", "structures")
RESULTS = os.path.join(HERE, "results")

SCHRODINGER = os.environ["SCHRODINGER"]
GLIDE = os.path.join(SCHRODINGER, "glide")
METHOD = "mininplace"

GLIDE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   {method}
FORCEFIELD   OPLS_2005
POSTDOCK   True
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
POSES_PER_LIG   1
CALC_INPUT_RMS   True
"""

CSV_FIELDS = {
    "r_i_docking_score": "docking_score", "r_i_glide_gscore": "gscore",
    "r_i_glide_lipo": "lipo", "r_i_glide_hbond": "hbond",
    "r_i_glide_metal": "metal", "r_i_glide_rewards": "rewards",
    "r_i_glide_evdw": "evdw", "r_i_glide_ecoul": "ecoul",
    "r_i_glide_erotb": "erotb", "r_i_glide_esite": "esite",
    "r_i_glide_emodel": "emodel", "r_i_glide_energy": "energy",
    "r_i_glide_einternal": "einternal", "r_i_glide_rmsd_to_input": "rmsd_to_input",
}


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Glide outside SLURM. "
                     "Use: sbatch analyses/09_pose_selection/scripts/05_glide_crystal.sh")


def main():
    require_compute_node()
    os.makedirs(RESULTS, exist_ok=True)
    out_csv = os.path.join(RESULTS, "glide_crystal_scores.csv")

    pdbs = sorted(pd.read_csv(os.path.join(RESULTS, "rtcnn_rescore_poses.csv"))
                  .pdb_id.unique())
    done = set()
    if os.path.exists(out_csv):
        done = set(pd.read_csv(out_csv).pdb_id)

    rows = []
    for k, pdb in enumerate(pdbs, 1):
        if pdb in done:
            print(f"[{k}/{len(pdbs)}] {pdb}: already scored", flush=True)
            continue
        wdir = os.path.join(A05_STRUCT, pdb)
        rec = {"pdb_id": pdb}
        print(f"[{k}/{len(pdbs)}] {pdb}: glide {METHOD}", flush=True)
        try:
            grid = glob.glob(os.path.join(wdir, f"{pdb}_grid*.zip"))
            lig = os.path.join(wdir, f"{pdb}_crystal_lig.sdf")
            if not grid:
                raise RuntimeError("no grid zip")
            if not os.path.exists(lig):
                raise RuntimeError("no crystal ligand sdf")

            job = f"{pdb}_{METHOD}"
            with open(os.path.join(wdir, f"{job}.in"), "w") as f:
                f.write(GLIDE_IN.format(grid=os.path.basename(grid[0]),
                                        lig=os.path.basename(lig), method=METHOD))
            subprocess.run([GLIDE, "-NOJOBID", "-WAIT", f"{job}.in"],
                           cwd=wdir, check=False, capture_output=True,
                           text=True, timeout=3600)

            csv_path = os.path.join(wdir, f"{job}.csv")
            if not os.path.exists(csv_path):
                raise RuntimeError("no results csv")
            df = pd.read_csv(csv_path)
            if df.empty:
                raise RuntimeError("no pose (ligand skipped)")
            row = df.iloc[0]
            for col, name in CSV_FIELDS.items():
                if col in df.columns:
                    try:
                        rec[name] = float(row[col])
                    except (TypeError, ValueError):
                        pass
            rec["status"] = str(row.get("docking_status", "?"))
            print(f"    gscore={rec.get('gscore')} emodel={rec.get('emodel')}", flush=True)
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"    ERROR: {e}", flush=True)
        rows.append(rec)
        # append after every structure so a walltime kill does not lose progress
        d = pd.DataFrame(rows)
        if os.path.exists(out_csv):
            d = pd.concat([pd.read_csv(out_csv), d], ignore_index=True)
            d = d.drop_duplicates("pdb_id", keep="last")
        d.to_csv(out_csv, index=False, float_format="%.4f")

    fin = pd.read_csv(out_csv)
    ok = (~fin.status.astype(str).str.startswith("error")).sum()
    print(f"\n{ok}/{len(pdbs)} scored. Wrote {out_csv}")
    return 0 if ok == len(pdbs) else 1


if __name__ == "__main__":
    sys.exit(main())
