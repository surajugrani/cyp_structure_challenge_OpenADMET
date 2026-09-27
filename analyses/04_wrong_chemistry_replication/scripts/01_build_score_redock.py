"""
Analysis 04: repeat analyses 01 and 02 under the heme representation that ACTIVATES
Glide's metal term, even though that representation is chemically wrong.

Analysis 03 found the metal term is 0.00 for every chemically defensible heme state,
including the correct ferric Fe(III) + porphyrin dianion. Only `fe3_monoanion`
(Fe +3, one pyrrole N neutralised) activates it, giving a fixed -2.30.

The question here is purely empirical: **does the wrong chemistry dock and score better?**
If it does, the result is right for the wrong reasons — that is understood going in and is
stated in the README, not discovered afterwards.

Controlled comparison. Everything that is chemically correct and *independent of the metal*
is kept exactly as in analyses 01/02 — same PrepWizard run, same deprotonated propionates,
same chain and ligand selection. Only the porphyrin/Fe state changes. This matters because
PrepWizard is non-deterministic about the heme (it gave 21 structures the dianion and 8SPD
the monoanion), so simply accepting its raw output would activate the metal term for an
arbitrary subset and confound the comparison.

Known confound, tracked not ignored: `fe3_monoanion` also shifts electrostatics hard
(8SPD ecoul -1.00 -> -7.12, ~1.2 GlideScore units at weight 0.19). All 13 Glide terms are
recorded for every pose so the metal and electrostatic contributions can be separated.

Usage (via scripts/run_all.sh):
    $SCHRODINGER/run python3 scripts/01_build_score_redock.py [--only 8SPD] [--poses 10]

Outputs (git-ignored): structures/<PDB>/
Outputs (committed):
    results/inplace_scores.csv    all terms, crystal pose, inplace + mininplace
    results/redock_poses.csv      all terms + RMSD for every redocked pose
    results/redock_summary.csv    per-ligand top-1 / best-of-N RMSD
"""

import argparse
import glob
import json
import os
import shutil
import subprocess

import numpy as np
import pandas as pd
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
A01_STRUCT = os.path.join(A01, "structures")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS_DIR = os.path.join(HERE, "results")

SCHRODINGER = os.environ["SCHRODINGER"]
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")

VARIANT = "fe3_monoanion"

GLIDE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   {method}
FORCEFIELD   OPLS_2005
POSTDOCK   True
POSES_PER_LIG   {poses}
POSTDOCK_NPOSE   {poses}
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
CALC_INPUT_RMS   True
"""

CSV_FIELDS = {
    "r_i_docking_score": "docking_score", "r_i_glide_gscore": "gscore",
    "r_i_glide_lipo": "lipo", "r_i_glide_hbond": "hbond", "r_i_glide_metal": "metal",
    "r_i_glide_rewards": "rewards", "r_i_glide_evdw": "evdw", "r_i_glide_ecoul": "ecoul",
    "r_i_glide_erotb": "erotb", "r_i_glide_esite": "esite", "r_i_glide_emodel": "emodel",
    "r_i_glide_energy": "energy", "r_i_glide_einternal": "einternal",
}


def require_compute_node(what):
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(f"REFUSING to run {what} outside SLURM. Use: sbatch scripts/run_all.sh")


def apply_variant(src_mae, dst_mae):
    """Fe -> +3 and neutralise one anionic pyrrole N (the metal-term-activating state)."""
    st = structure.StructureReader.read(src_mae)
    fe = [a for a in st.atom if a.element == "Fe" and a.pdbres.strip() == "HEM"]
    if len(fe) != 1:
        raise RuntimeError(f"expected 1 heme Fe, found {len(fe)}")
    fe = fe[0]
    pyrrole_n = [b for b in fe.bonded_atoms if b.element == "N" and b.pdbres.strip() == "HEM"]
    if len(pyrrole_n) != 4:
        raise RuntimeError(f"expected 4 pyrrole N, found {len(pyrrole_n)}")
    fe.formal_charge = 3
    anionic = sorted([n for n in pyrrole_n if n.formal_charge == -1], key=lambda a: a.pdbname)
    if not anionic:
        raise RuntimeError("no anionic pyrrole N to neutralise")
    anionic[0].formal_charge = 0
    st.write(dst_mae)
    heme = [a for a in st.atom if a.pdbres.strip() == "HEM"]
    return {"fe_formal_charge": 3, "heme_formal_charge": sum(a.formal_charge for a in heme)}


def split(complex_mae, pdb, wdir, tag):
    """Receptor (.mae) and ligand (.sdf); ligand identified from analysis 01's prep_info."""
    with open(os.path.join(A01_STRUCT, pdb, "prep_info.json")) as f:
        info = json.load(f)
    st = structure.StructureReader.read(complex_mae)
    lig_idx = [a.index for a in st.atom
               if a.chain == info["chain"] and a.resnum == info["score_resnum"]
               and a.pdbres.strip() == info["comp_id"]]
    if not lig_idx:
        raise RuntimeError("ligand not located in variant complex")
    lig_sdf = os.path.join(wdir, f"{tag}_lig.sdf")
    st.extract(lig_idx, copy_props=True).write(lig_sdf)
    rec = st.copy()
    rec.deleteAtoms(lig_idx)
    rec_mae = os.path.join(wdir, f"{tag}_rec.mae")
    rec.write(rec_mae)
    return rec_mae, lig_sdf


def build_grid(rec_mae, lig_sdf, wdir, tag):
    grid = os.path.join(wdir, f"{tag}_grid-gridgen.zip")
    if os.path.exists(grid):
        return grid
    rec = structure.StructureReader.read(rec_mae)
    lig = structure.StructureReader.read(lig_sdf)
    for a in lig.atom:
        a.pdbres, a.chain = "LIG ", "L"
    merged = os.path.join(wdir, f"{tag}_reclig.mae")
    rec.merge(lig).write(merged)
    require_compute_node("Glide grid generation")
    subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT", "-rec_file", os.path.basename(merged),
                    "-lig_asl", "res. LIG", "-j", f"{tag}_grid"], cwd=wdir, check=True)
    if not os.path.exists(grid):
        found = glob.glob(os.path.join(wdir, f"{tag}_grid*.zip"))
        if not found:
            raise RuntimeError("grid not produced")
        grid = found[0]
    return grid


def run_glide(grid, lig_sdf, wdir, job, method, poses):
    # Glide resolves LIGANDFILE relative to cwd, so the ligand must live in wdir.
    local_lig = os.path.join(wdir, os.path.basename(lig_sdf))
    if os.path.abspath(local_lig) != os.path.abspath(lig_sdf):
        shutil.copy2(lig_sdf, local_lig)
    in_file = os.path.join(wdir, f"{job}.in")
    with open(in_file, "w") as f:
        f.write(GLIDE_IN.format(grid=os.path.basename(grid),
                                lig=os.path.basename(local_lig), method=method, poses=poses))
    require_compute_node("Glide")
    subprocess.run([GLIDE, "-NOJOBID", "-WAIT", os.path.basename(in_file)],
                   cwd=wdir, check=False)
    csv_path = os.path.join(wdir, f"{job}.csv")
    if not os.path.exists(csv_path):
        return pd.DataFrame()
    return pd.read_csv(csv_path)


def terms(row, prefix=""):
    out = {}
    for col, name in CSV_FIELDS.items():
        if col in row.index:
            try:
                out[f"{prefix}{name}"] = float(row[col])
            except (TypeError, ValueError):
                pass
    return out


def heavy_xyz(st):
    return np.array([a.xyz for a in st.atom if a.element != "H"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    ap.add_argument("--poses", type=int, default=10)
    args = ap.parse_args()

    cur = pd.read_csv(os.path.join(A01, "data", "curated_set.csv"))
    if args.only:
        cur = cur[cur["pdb_id"].isin(args.only.split(","))]
    os.makedirs(RESULTS_DIR, exist_ok=True)

    inplace_rows, pose_rows, summary_rows = [], [], []
    for _, row in cur.iterrows():
        pdb, comp = row["pdb_id"], row["comp_id"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        os.makedirs(wdir, exist_ok=True)
        tag = f"{pdb}_{VARIANT}"
        print(f"\n== {pdb} ({comp})")
        rec = {"pdb_id": pdb, "comp_id": comp}
        try:
            src = os.path.join(A01_STRUCT, pdb, f"{pdb}_prepped.mae")
            if not os.path.exists(src):
                raise RuntimeError("run analysis 01 first")
            complex_mae = os.path.join(wdir, f"{tag}_complex.mae")
            rec.update(apply_variant(src, complex_mae))
            rec_mae, lig_sdf = split(complex_mae, pdb, wdir, tag)
            grid = build_grid(rec_mae, lig_sdf, wdir, tag)

            # --- analysis 01 equivalent: score the crystal pose ---
            for method in ("inplace", "mininplace"):
                df = run_glide(grid, lig_sdf, wdir, f"{tag}_{method}", method, 1)
                if len(df):
                    rec.update(terms(df.iloc[0], prefix=f"{method}_"))
                    rec[f"{method}_status"] = str(df.iloc[0].get("docking_status", "?"))
            print(f"    inplace metal={rec.get('inplace_metal')} "
                  f"score={rec.get('inplace_docking_score')}")

            # --- analysis 02 equivalent: redock ---
            rd = run_glide(grid, lig_sdf, wdir, f"{tag}_redock", "confgen", args.poses)
            pv = os.path.join(wdir, f"{tag}_redock_pv.mae")
            if len(rd) and os.path.exists(pv):
                ref_xyz = heavy_xyz(structure.StructureReader.read(lig_sdf))
                rmsds, naive = [], []
                for i, st in enumerate(structure.StructureReader(pv)):
                    if i == 0:
                        continue
                    rmsds.append(st.property.get("r_i_glide_rmsd_to_input"))
                    pz = heavy_xyz(st)
                    naive.append(float(np.sqrt(((ref_xyz - pz) ** 2).sum(axis=1).mean()))
                                 if pz.shape == ref_xyz.shape else np.nan)
                if any(r is None for r in rmsds):
                    rmsds = naive
                for k in range(len(rmsds)):
                    pr = {"pdb_id": pdb, "comp_id": comp, "rank": k + 1,
                          "rmsd": rmsds[k], "rmsd_naive": naive[k]}
                    if k < len(rd):
                        pr.update(terms(rd.iloc[k]))
                    pose_rows.append(pr)
                best = int(np.nanargmin(rmsds))
                summary_rows.append({"pdb_id": pdb, "comp_id": comp, "n_poses": len(rmsds),
                                     "rmsd_top1": rmsds[0], "rmsd_best": rmsds[best],
                                     "rank_of_best": best + 1,
                                     "redock_metal_top1": float(rd.iloc[0]["r_i_glide_metal"])
                                     if "r_i_glide_metal" in rd.columns else None})
                print(f"    redock top1={rmsds[0]:.2f} best={rmsds[best]:.2f} "
                      f"(rank {best + 1}/{len(rmsds)})")
            rec["status"] = "ok"
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"    ERROR: {e}")
        inplace_rows.append(rec)

    pd.DataFrame(inplace_rows).to_csv(os.path.join(RESULTS_DIR, "inplace_scores.csv"),
                                      index=False, float_format="%.4f")
    pd.DataFrame(pose_rows).to_csv(os.path.join(RESULTS_DIR, "redock_poses.csv"),
                                   index=False, float_format="%.4f")
    s = pd.DataFrame(summary_rows)
    s.to_csv(os.path.join(RESULTS_DIR, "redock_summary.csv"), index=False, float_format="%.3f")

    if len(s):
        print(f"\n=== Pose recovery ({VARIANT}), n = {len(s)} ===")
        for label, col in (("top-1", "rmsd_top1"), ("best-of-N", "rmsd_best")):
            for cut in (2.0, 1.0):
                n = (s[col] <= cut).sum()
                print(f"  {label:<10} RMSD <= {cut:.1f} A: {n:2d}/{len(s)} ({100*n/len(s):.0f} %)")
            print(f"  {label:<10} median: {s[col].median():.2f} A")


if __name__ == "__main__":
    main()
