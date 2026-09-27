"""
Analysis 05 step 2: prepare every CYP3A4 structure and redock its ligand FROM SMILES.

Difference from analysis 02, and the point of this analysis: the ligand is rebuilt from a
SMILES string with LigPrep rather than taken from the crystal structure. Analysis 02's
`confgen` regenerated torsions but inherited crystal bond lengths, angles and
stereochemistry, which flatters the result. Here the only crystal-derived inputs are the
receptor and the grid centre — the same situation as the 20 blinded challenge compounds,
where we have a SMILES string and nothing else.

Expect worse numbers than analysis 02's 18 % top-1. The gap is the point: it measures how
much the crystal-geometry shortcut was inflating the estimate.

Consequence for RMSD. LigPrep output has a different atom ordering from the crystal
ligand, so Glide's `r_i_glide_rmsd_to_input` is useless here — it measures distance from
the LigPrep starting conformer, not from the crystal pose. RMSD is therefore computed with
**spyrmsd** (graph-matched, symmetry-corrected, no superposition). Validated against
Glide's own values on analysis 05's crystal-geometry poses: exact agreement on 6/8 test
poses, and <= 0.16 A on the two that differ (both on TPF, a ligand with genuine symmetry,
where spyrmsd is the stricter of the two).

Runs SERIALLY - one Schrodinger licence token. Resumable and incremental: a structure with
completed results is skipped, and results are appended after every structure, so an
interrupted run continues where it stopped.

Usage (via scripts/run_all.sh):
    $SCHRODINGER/run python3 scripts/02_prep_dock.py [--poses 10] [--limit N] [--only IDS]

Outputs (git-ignored): structures/<PDB>/
Outputs (committed):
    results/prep_status.csv    per-structure prep record and heme checks
    results/redock_poses.csv   all Glide terms + spyrmsd RMSD for every pose
    results/redock_summary.csv per-structure top-1 / best-of-N RMSD
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from schrodinger import structure
from spyrmsd import rmsd as spy

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
DATA = os.path.join(HERE, "data")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")

# Reuse analysis 01's preparation verbatim - identical PrepWizard settings, normalize_heme
# and ligand selection, so receptors are comparable with analyses 01-04.
sys.path.insert(0, os.path.join(A01, "scripts"))
from importlib import import_module
_p = import_module("03_prep_structures")
download_pdb, clean_structure = _p.download_pdb, _p.clean_structure
run_prepwizard, normalize_heme, check_prepped = _p.run_prepwizard, _p.normalize_heme, _p.check_prepped

SCHRODINGER = os.environ["SCHRODINGER"]
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")
LIGPREP = os.path.join(SCHRODINGER, "ligprep")

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


def append_csv(path, rows):
    if not rows:
        return
    pd.DataFrame(rows).to_csv(path, mode="a", header=not os.path.exists(path),
                              index=False, float_format="%.4f")


def graph_arrays(st):
    """Heavy-atom coordinates, atomic numbers and adjacency matrix for spyrmsd."""
    idx = [a.index for a in st.atom if a.element != "H"]
    pos = {j: i for i, j in enumerate(idx)}
    xyz = np.array([st.atom[j].xyz for j in idx])
    z = np.array([st.atom[j].atomic_number for j in idx])
    adj = np.zeros((len(idx), len(idx)), dtype=int)
    for b in st.bond:
        i, j = b.atom1.index, b.atom2.index
        if i in pos and j in pos:
            adj[pos[i], pos[j]] = adj[pos[j], pos[i]] = 1
    return xyz, z, adj


def ligprep_from_smiles(smiles, wdir, pdb, comp):
    """3D ligand states from SMILES. Epik ionization/tautomers at pH 7.4.

    -g respects the stereochemistry given in the SMILES rather than enumerating it; all
    92 ligands in this set have stereo-defined SMILES from the RCSB chemical component
    dictionary, so there is nothing to enumerate.
    """
    smi = os.path.join(wdir, f"{pdb}_lig.smi")
    with open(smi, "w") as f:
        f.write(f"{smiles} {pdb}_{comp}\n")
    out = os.path.join(wdir, f"{pdb}_ligprep.maegz")
    if os.path.exists(out):
        return out
    require_compute_node("LigPrep")
    subprocess.run([LIGPREP, "-NOJOBID", "-WAIT", "-ismi", os.path.basename(smi),
                    "-omae", os.path.basename(out), "-epikx", "-ph", "7.4", "-pht", "2.0",
                    "-g", "-bff", "14"], cwd=wdir, check=False)
    if not os.path.exists(out):
        raise RuntimeError("LigPrep produced no output")
    return out


def prepare(pdb, comp, chain, wdir):
    pdb_path = os.path.join(wdir, f"{pdb}.pdb")
    download_pdb(pdb, pdb_path)
    clean, info = clean_structure(pdb_path, chain, comp)
    clean_mae = os.path.join(wdir, f"{pdb}_clean.mae")
    clean.write(clean_mae)
    prepped = os.path.join(wdir, f"{pdb}_prepped.mae")
    if not os.path.exists(prepped):
        require_compute_node("PrepWizard")
        info["prep_mode"] = run_prepwizard(clean_mae, prepped, wdir)
    info.update({f"heme_fixed_{k}": v for k, v in normalize_heme(prepped).items()})
    info.update(check_prepped(prepped, info))
    return prepped, info


def split(prepped, info, wdir, pdb):
    """Receptor, plus the crystal ligand kept only as the RMSD reference."""
    st = structure.StructureReader.read(prepped)
    lig_idx = [a.index for a in st.atom
               if a.chain == info["chain"] and a.resnum == info["score_resnum"]
               and a.pdbres.strip() == info["comp_id"]]
    if not lig_idx:
        raise RuntimeError("crystal ligand copy not found in prepped structure")
    ref_sdf = os.path.join(wdir, f"{pdb}_crystal_lig.sdf")
    st.extract(lig_idx, copy_props=True).write(ref_sdf)
    rec = st.copy()
    rec.deleteAtoms(lig_idx)
    rec_mae = os.path.join(wdir, f"{pdb}_rec.mae")
    rec.write(rec_mae)
    return rec_mae, ref_sdf, len(lig_idx)


def build_grid(rec_mae, ref_sdf, wdir, pdb):
    """Grid centred on the crystal ligand position (site is given, as in analysis 02)."""
    grid = os.path.join(wdir, f"{pdb}_grid-gridgen.zip")
    if os.path.exists(grid):
        return grid
    rec = structure.StructureReader.read(rec_mae)
    lig = structure.StructureReader.read(ref_sdf)
    for a in lig.atom:
        a.pdbres, a.chain = "LIG ", "L"
    merged = os.path.join(wdir, f"{pdb}_reclig.mae")
    rec.merge(lig).write(merged)
    require_compute_node("Glide grid generation")
    subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT", "-rec_file", os.path.basename(merged),
                    "-lig_asl", "res. LIG", "-j", f"{pdb}_grid"], cwd=wdir, check=True)
    if not os.path.exists(grid):
        found = glob.glob(os.path.join(wdir, f"{pdb}_grid*.zip"))
        if not found:
            raise RuntimeError("grid not produced")
        grid = found[0]
    return grid


def dock(grid, ligprep_mae, wdir, pdb, poses):
    local = os.path.join(wdir, os.path.basename(ligprep_mae))
    if os.path.abspath(local) != os.path.abspath(ligprep_mae):
        shutil.copy2(ligprep_mae, local)
    in_file = os.path.join(wdir, f"{pdb}_dock.in")
    with open(in_file, "w") as f:
        f.write(GLIDE_IN.format(grid=os.path.basename(grid),
                                lig=os.path.basename(local), poses=poses))
    require_compute_node("Glide")
    subprocess.run([GLIDE, "-NOJOBID", "-WAIT", os.path.basename(in_file)],
                   cwd=wdir, check=False)
    csv_path = os.path.join(wdir, f"{pdb}_dock.csv")
    return pd.read_csv(csv_path) if os.path.exists(csv_path) else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poses", type=int, default=10)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--only")
    args = ap.parse_args()

    s = pd.read_csv(os.path.join(DATA, "structure_set.csv"))
    smiles = json.load(open(os.path.join(DATA, "ligand_smiles.json")))
    if args.only:
        s = s[s["pdb_id"].isin(args.only.split(","))]
    os.makedirs(RESULTS, exist_ok=True)

    prep_csv = os.path.join(RESULTS, "prep_status.csv")
    pose_csv = os.path.join(RESULTS, "redock_poses.csv")
    summ_csv = os.path.join(RESULTS, "redock_summary.csv")

    done = set()
    if os.path.exists(summ_csv):
        done = set(pd.read_csv(summ_csv)["pdb_id"])
        print(f"resuming: {len(done)} structures already complete")
    todo = [r for _, r in s.iterrows() if r["pdb_id"] not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(todo)} structures to process (serial - one Schrodinger licence)\n")

    for i, row in enumerate(todo, 1):
        pdb, comp, chain = row["pdb_id"], row["comp_id"], row["chain"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        os.makedirs(wdir, exist_ok=True)
        print(f"[{i}/{len(todo)}] {pdb} ({comp}, chain {chain})", flush=True)
        rec = {"pdb_id": pdb, "comp_id": comp,
               "has_affinity": bool(row["has_affinity"]), "resolution": row["resolution"]}
        try:
            smi = (smiles.get(comp) or {}).get("smiles_stereo")
            if not smi:
                raise RuntimeError(f"no SMILES for {comp}")
            rec["smiles"] = smi

            prepped, info = prepare(pdb, comp, chain, wdir)
            rec.update({k: v for k, v in info.items() if not isinstance(v, list)})
            rec_mae, ref_sdf, n_lig = split(prepped, info, wdir, pdb)
            rec["n_crystal_lig_atoms"] = n_lig
            grid = build_grid(rec_mae, ref_sdf, wdir, pdb)

            lp = ligprep_from_smiles(smi, wdir, pdb, comp)
            n_states = sum(1 for _ in structure.StructureReader(lp))
            rec["n_ligprep_states"] = n_states
            print(f"    LigPrep: {n_states} state(s)", flush=True)

            rd = dock(grid, lp, wdir, pdb, args.poses)
            pv = os.path.join(wdir, f"{pdb}_dock_pv.mae")
            if not len(rd) or not os.path.exists(pv):
                raise RuntimeError("no poses produced")

            rx, rz, ra = graph_arrays(structure.StructureReader.read(ref_sdf))
            rmsds = []
            for k, st in enumerate(structure.StructureReader(pv)):
                if k == 0:
                    continue
                px, pz, pa = graph_arrays(st)
                try:
                    rmsds.append(float(spy.symmrmsd(rx, px, rz, pz, ra, pa, minimize=False)))
                except Exception as e:
                    print(f"      spyrmsd failed on pose {k}: {e}", flush=True)
                    rmsds.append(np.nan)
            if not rmsds or np.all(np.isnan(rmsds)):
                raise RuntimeError("no RMSD could be computed")

            pose_rows = []
            for k in range(len(rmsds)):
                pr = {"pdb_id": pdb, "comp_id": comp, "rank": k + 1, "rmsd": rmsds[k]}
                if k < len(rd):
                    for col, name in CSV_FIELDS.items():
                        if col in rd.columns:
                            try:
                                pr[name] = float(rd.iloc[k][col])
                            except (TypeError, ValueError):
                                pass
                pose_rows.append(pr)
            append_csv(pose_csv, pose_rows)

            best = int(np.nanargmin(rmsds))
            append_csv(summ_csv, [{
                "pdb_id": pdb, "comp_id": comp, "has_affinity": bool(row["has_affinity"]),
                "resolution": row["resolution"], "n_crystal_lig_atoms": n_lig,
                "n_ligprep_states": n_states, "n_poses": len(rmsds),
                "rmsd_top1": rmsds[0], "rmsd_best": rmsds[best], "rank_of_best": best + 1,
                "score_top1": float(rd.iloc[0]["r_i_docking_score"])
                if "r_i_docking_score" in rd.columns else None,
                "status": "ok"}])
            rec["status"] = "ok"
            print(f"    top1={rmsds[0]:.2f} best={rmsds[best]:.2f} "
                  f"(rank {best + 1}/{len(rmsds)})", flush=True)
        except Exception as e:
            rec["status"] = f"error: {e}"
            append_csv(summ_csv, [{"pdb_id": pdb, "comp_id": comp,
                                   "has_affinity": bool(row["has_affinity"]),
                                   "status": f"error: {e}"}])
            print(f"    ERROR: {e}", flush=True)
        append_csv(prep_csv, [rec])

    if os.path.exists(summ_csv):
        d = pd.read_csv(summ_csv)
        ok = d[d["status"] == "ok"]
        print(f"\n=== {len(ok)}/{len(d)} succeeded ===")
        for label, col in (("top-1", "rmsd_top1"), ("best-of-N", "rmsd_best")):
            for cut in (2.0, 1.0):
                n = (ok[col] <= cut).sum()
                print(f"  {label:<10} RMSD <= {cut:.1f} A: {n:3d}/{len(ok)} "
                      f"({100*n/len(ok):.0f} %)")
            print(f"  {label:<10} median: {ok[col].median():.2f} A")


if __name__ == "__main__":
    main()
