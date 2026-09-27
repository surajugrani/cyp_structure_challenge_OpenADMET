"""
Step 4: Build a Glide grid per structure and score the crystal pose.

For each prepped complex from step 3:
  1. Split the prepped .mae into receptor (protein + HEM + any extra ligand copies)
     and the scoring ligand copy (written as .sdf).
  2. generate_glide_grids centred on the scoring ligand copy (removed from the grid).
  3. Glide with DOCKING_METHOD inplace    -> score of the crystal pose as-is
  4. Glide with DOCKING_METHOD mininplace -> score after local minimization in the site
  5. Collect docking score / GlideScore / component terms from the output .rept files.

Both scoring modes are run because `inplace` uses the exact crystal coordinates (which
can clash after H addition) while `mininplace` relaxes the pose locally in the pocket.

Usage (normally via scripts/run_prep_and_score.sh):
    $SCHRODINGER/run python3 scripts/04_grid_and_score.py [--only 2V0M,4K9T]

Outputs (git-ignored): structures/<PDB>/grid.zip, *_inplace*, *_mininplace*
Outputs (committed):   results/glide_scores.csv
"""

import argparse
import glob
import os
import subprocess

import pandas as pd
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(HERE, "data")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS_DIR = os.path.join(HERE, "results")

SCHRODINGER = os.environ["SCHRODINGER"]
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")

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

# Columns to pull out of the Glide results CSV (<jobname>.csv)
CSV_FIELDS = {
    "r_i_docking_score": "docking_score",
    "r_i_glide_gscore": "gscore",
    "r_i_glide_lipo": "lipo",
    "r_i_glide_hbond": "hbond",
    "r_i_glide_metal": "metal",
    "r_i_glide_rewards": "rewards",
    "r_i_glide_evdw": "evdw",
    "r_i_glide_ecoul": "ecoul",
    "r_i_glide_erotb": "erotb",
    "r_i_glide_esite": "esite",
    "r_i_glide_emodel": "emodel",
    "r_i_glide_energy": "energy",
    "r_i_glide_einternal": "einternal",
}


def split_receptor_ligand(prepped_mae, info, wdir, pdb):
    """Write receptor (.mae, ligand copy removed) and scoring ligand (.sdf)."""
    st = structure.StructureReader.read(prepped_mae)
    lig_idx = [a.index for a in st.atom
               if a.chain == info["chain"] and a.resnum == info["score_resnum"]
               and a.pdbres.strip() == info["comp_id"]]
    if not lig_idx:
        raise RuntimeError("scoring ligand copy not found in prepped structure")

    lig_st = st.extract(lig_idx, copy_props=True)
    lig_st.title = f"{pdb}_{info['comp_id']}"
    lig_sdf = os.path.join(wdir, f"{pdb}_lig.sdf")
    lig_st.write(lig_sdf)

    rec_st = st.copy()
    rec_st.deleteAtoms(lig_idx)
    rec_mae = os.path.join(wdir, f"{pdb}_rec.mae")
    rec_st.write(rec_mae)
    return rec_mae, lig_sdf, len(lig_idx)


def make_grid(rec_mae, lig_sdf, wdir, pdb):
    """Grid centred on the crystal ligand, which is merged in then removed by -lig_asl."""
    # generate_glide_grids appends "-gridgen" to the job name, so the output is
    # <jobname>-gridgen.zip rather than <jobname>.zip.
    grid_zip = os.path.join(wdir, f"{pdb}_grid-gridgen.zip")
    if os.path.exists(grid_zip):
        return grid_zip

    # Merge receptor + ligand so the ligand defines the box centre, then -lig_asl
    # removes those atoms from the grid.
    rec = structure.StructureReader.read(rec_mae)
    lig = structure.StructureReader.read(lig_sdf)
    for a in lig.atom:
        a.pdbres = "LIG "
        a.chain = "L"
    merged = rec.merge(lig)
    merged_mae = os.path.join(wdir, f"{pdb}_reclig.mae")
    merged.write(merged_mae)

    cmd = [GEN_GRIDS, "-NOJOBID", "-WAIT",
           "-rec_file", os.path.basename(merged_mae),
           "-lig_asl", "res. LIG",
           "-j", f"{pdb}_grid"]
    print("    $", " ".join(cmd))
    subprocess.run(cmd, cwd=wdir, check=True)
    if not os.path.exists(grid_zip):
        found = glob.glob(os.path.join(wdir, f"{pdb}_grid*.zip"))
        if not found:
            raise RuntimeError("grid not produced")
        grid_zip = found[0]
    return grid_zip


def run_glide(grid_zip, lig_sdf, method, wdir, pdb):
    job = f"{pdb}_{method}"
    in_file = os.path.join(wdir, f"{job}.in")
    with open(in_file, "w") as f:
        f.write(GLIDE_IN.format(grid=os.path.basename(grid_zip),
                                lig=os.path.basename(lig_sdf), method=method))
    cmd = [GLIDE, "-NOJOBID", "-WAIT", os.path.basename(in_file)]
    print("    $", " ".join(cmd))
    subprocess.run(cmd, cwd=wdir, check=False)
    return parse_glide_csv(os.path.join(wdir, f"{job}.csv"), method)


def parse_glide_csv(csv_path, method):
    """Read the Glide results CSV (one row per pose; we keep the first)."""
    if not os.path.exists(csv_path):
        return {f"{method}_status": "no csv"}
    df = pd.read_csv(csv_path)
    if df.empty:
        return {f"{method}_status": "no pose (ligand skipped)"}

    row = df.iloc[0]
    out = {f"{method}_status": str(row.get("docking_status", "?"))}
    for col, name in CSV_FIELDS.items():
        if col in df.columns:
            try:
                out[f"{method}_{name}"] = float(row[col])
            except (TypeError, ValueError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated PDB IDs")
    args = ap.parse_args()

    cur = pd.read_csv(os.path.join(DATA_DIR, "curated_set.csv"))
    if args.only:
        cur = cur[cur["pdb_id"].isin(args.only.split(","))]
    prep = pd.read_csv(os.path.join(RESULTS_DIR, "prep_summary.csv"))
    prep = prep.set_index("pdb_id")

    rows = []
    for _, row in cur.iterrows():
        pdb = row["pdb_id"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        prepped = os.path.join(wdir, f"{pdb}_prepped.mae")
        print(f"\n== {pdb} ({row['comp_id']})")
        out = {"pdb_id": pdb, "comp_id": row["comp_id"]}
        try:
            if not os.path.exists(prepped):
                raise RuntimeError("no prepped structure (run step 3 first)")
            info = {"chain": prep.loc[pdb, "chain"],
                    "score_resnum": int(prep.loc[pdb, "score_resnum"]),
                    "comp_id": row["comp_id"]}
            rec_mae, lig_sdf, n_lig = split_receptor_ligand(prepped, info, wdir, pdb)
            out["n_lig_atoms"] = n_lig
            grid_zip = make_grid(rec_mae, lig_sdf, wdir, pdb)
            for method in ("inplace", "mininplace"):
                out.update(run_glide(grid_zip, lig_sdf, method, wdir, pdb))
            print("   ", {k: v for k, v in out.items() if "docking_score" in k or "status" in k})
        except Exception as e:
            print(f"    ERROR: {e}")
            out["error"] = str(e)
        rows.append(out)

    df = pd.DataFrame(rows)
    # When --only is used, merge into any existing table rather than clobbering it
    # with just the subset that was re-run.
    out_csv = os.path.join(RESULTS_DIR, "glide_scores.csv")
    if args.only and os.path.exists(out_csv):
        old = pd.read_csv(out_csv)
        old = old[~old["pdb_id"].isin(df["pdb_id"])]
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(out_csv, index=False)
    ok = df.get("inplace_docking_score", pd.Series(dtype=float)).notna().sum()
    print(f"\n{ok}/{len(df)} scored. Wrote results/glide_scores.csv")


if __name__ == "__main__":
    main()
