"""
Analysis 03: does the heme representation control Glide's metal term?

Analyses 01/02 found the metal term is identically 0.00 under the chemically correct heme
(Fe as PrepWizard leaves it, +2; porphyrin dianion), and that Glide then displaces ligands
off the iron by a median +0.73 A. The term is decisive at the relevant scale: vdW is
weighted ~0.046 in GlideScore while metal enters at ~1.0.

This script rebuilds the receptor under several heme variants, regenerates the grid for
each, rescores the crystal pose in place, and records **every** Glide term.

Variants (applied to the already-prepped complex from analysis 01; no PrepWizard re-run):

  fe2_dianion   baseline, = analyses 01/02. Fe +2, porphyrin dianion (2 N at -1, no N-H).
  fe3_dianion   Fe +3. CYP3A4's resting state is ferric Fe(III); PrepWizard leaves +2.
                The chemically motivated lever — does a ferric iron activate the term?
  fe2_monoanion Fe +2, one pyrrole N set to neutral (no H added). Isolates whether the
                pyrrole *charge* matters, separately from the N-H itself.
  fe3_monoanion both changes together.

Note that changing Fe's formal charge also shifts grid electrostatics, so the metal term is
recorded explicitly rather than inferred from the total score. Formal charge on a metal in
Schrodinger's zero-order-bond scheme is partly bookkeeping, so the question here is
strictly empirical: *does it change the metal term?*

Usage (via scripts/run_variants.sh):
    $SCHRODINGER/run python3 scripts/01_variants_and_score.py --pdbs 8SPD,4K9X,2V0M

Outputs (git-ignored): structures/<PDB>/<variant>/
Outputs (committed):   results/variant_scores.csv   (all Glide terms, one row per variant)
"""

import argparse
import glob
import json
import os
import shutil
import subprocess

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

VARIANTS = ["fe2_dianion", "fe3_dianion", "fe2_monoanion", "fe3_monoanion"]

GLIDE_IN = """GRIDFILE   {grid}
LIGANDFILE   {lig}
DOCKING_METHOD   inplace
FORCEFIELD   OPLS_2005
POSTDOCK   True
POSES_PER_LIG   1
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
CALC_INPUT_RMS   True
"""

# Every Glide term, not just the docking score (needed to see the metal term directly)
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
    raise SystemExit(f"REFUSING to run {what} outside SLURM. Use: sbatch scripts/run_variants.sh")


def apply_variant(src_mae, dst_mae, variant):
    """Write a copy of the prepped complex with the heme in the requested state."""
    st = structure.StructureReader.read(src_mae)
    fe = [a for a in st.atom if a.element == "Fe" and a.pdbres.strip() == "HEM"]
    if len(fe) != 1:
        raise RuntimeError(f"expected 1 heme Fe, found {len(fe)}")
    fe = fe[0]
    pyrrole_n = [b for b in fe.bonded_atoms if b.element == "N" and b.pdbres.strip() == "HEM"]
    if len(pyrrole_n) != 4:
        raise RuntimeError(f"expected 4 pyrrole N, found {len(pyrrole_n)}")

    fe.formal_charge = 3 if variant.startswith("fe3") else 2

    if variant.endswith("monoanion"):
        # Neutralise one of the two anionic pyrrole N (no H added - this isolates the
        # charge effect from the N-H itself).
        anionic = sorted([n for n in pyrrole_n if n.formal_charge == -1],
                         key=lambda a: a.pdbname)
        if not anionic:
            raise RuntimeError("no anionic pyrrole N to neutralise")
        anionic[0].formal_charge = 0

    st.write(dst_mae)
    heme = [a for a in st.atom if a.pdbres.strip() == "HEM"]
    return {"fe_formal_charge": fe.formal_charge,
            "heme_formal_charge": sum(a.formal_charge for a in heme),
            "pyrrole_N_charges": ",".join(str(n.formal_charge) for n in pyrrole_n)}


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


def score(grid, lig_sdf, wdir, tag):
    # Glide runs with cwd=wdir and resolves LIGANDFILE relative to it, so the ligand must
    # be copied in - the source SDF lives in analysis 01's directory.
    local_lig = os.path.join(wdir, os.path.basename(lig_sdf))
    if not os.path.exists(local_lig):
        shutil.copy2(lig_sdf, local_lig)
    in_file = os.path.join(wdir, f"{tag}_score.in")
    with open(in_file, "w") as f:
        f.write(GLIDE_IN.format(grid=os.path.basename(grid), lig=os.path.basename(local_lig)))
    require_compute_node("Glide")
    subprocess.run([GLIDE, "-NOJOBID", "-WAIT", os.path.basename(in_file)],
                   cwd=wdir, check=False)
    csv_path = os.path.join(wdir, f"{tag}_score.csv")
    if not os.path.exists(csv_path):
        return {"status": "no csv"}
    df = pd.read_csv(csv_path)
    if df.empty:
        return {"status": "no pose"}
    row = df.iloc[0]
    out = {"status": str(row.get("docking_status", "?"))}
    for col, name in CSV_FIELDS.items():
        if col in df.columns:
            try:
                out[name] = float(row[col])
            except (TypeError, ValueError):
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbs", required=True, help="comma-separated PDB IDs")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    args = ap.parse_args()
    os.makedirs(RESULTS_DIR, exist_ok=True)

    rows = []
    for pdb in args.pdbs.split(","):
        src = os.path.join(A01_STRUCT, pdb, f"{pdb}_prepped.mae")
        lig_src = os.path.join(A01_STRUCT, pdb, f"{pdb}_lig.sdf")
        for variant in args.variants.split(","):
            wdir = os.path.join(STRUCT_DIR, pdb, variant)
            os.makedirs(wdir, exist_ok=True)
            tag = f"{pdb}_{variant}"
            print(f"\n== {pdb} / {variant}")
            rec = {"pdb_id": pdb, "variant": variant}
            try:
                for path in (src, lig_src):
                    if not os.path.exists(path):
                        raise RuntimeError(f"missing {path} - run analysis 01 first")
                complex_mae = os.path.join(wdir, f"{tag}_complex.mae")
                rec.update(apply_variant(src, complex_mae, variant))

                # Split receptor / ligand from the variant complex.
                # NOTE: identify the ligand from analysis 01's prep_info.json, NOT from the
                # SDF - SDF format carries no chain/residue information, so matching on
                # (chain, resnum, pdbres) read back from the SDF always fails.
                with open(os.path.join(A01_STRUCT, pdb, "prep_info.json")) as f:
                    info = json.load(f)
                st = structure.StructureReader.read(complex_mae)
                lig_idx = [a.index for a in st.atom
                           if a.chain == info["chain"]
                           and a.resnum == info["score_resnum"]
                           and a.pdbres.strip() == info["comp_id"]]
                if not lig_idx:
                    raise RuntimeError(
                        f"ligand {info['comp_id']} {info['chain']}:{info['score_resnum']} "
                        "not located in variant complex")
                rec_st = st.copy()
                rec_st.deleteAtoms(lig_idx)
                rec_mae = os.path.join(wdir, f"{tag}_rec.mae")
                rec_st.write(rec_mae)

                grid = build_grid(rec_mae, lig_src, wdir, tag)
                rec.update(score(grid, lig_src, wdir, tag))
                print(f"    heme={rec['heme_formal_charge']:+d} Fe={rec['fe_formal_charge']:+d} "
                      f"metal={rec.get('metal')} score={rec.get('docking_score')}")
            except Exception as e:
                rec["status"] = f"error: {e}"
                print(f"    ERROR: {e}")
            rows.append(rec)

    df = pd.DataFrame(rows)
    out_csv = os.path.join(RESULTS_DIR, "variant_scores.csv")
    if os.path.exists(out_csv):
        old = pd.read_csv(out_csv)
        key = old["pdb_id"].astype(str) + "|" + old["variant"].astype(str)
        new = set(df["pdb_id"].astype(str) + "|" + df["variant"].astype(str))
        df = pd.concat([old[~key.isin(new)], df], ignore_index=True)
    df.to_csv(out_csv, index=False, float_format="%.4f")

    print("\n=== Metal term by variant ===")
    ok = df[df["status"].astype(str).str.startswith(("Done", "ok"))]
    if len(ok):
        print(ok.pivot_table(index="pdb_id", columns="variant", values="metal").to_string())
    print(f"\nWrote {out_csv}")


if __name__ == "__main__":
    main()
