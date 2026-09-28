"""
Analysis 06 step 2: dock with ICM-Pro and measure pose recovery against the crystal pose.

Comparison design (see README):
  * receptor  - ICM prepares it itself from the PRE-PrepWizard cleaned structure
  * ligand    - the SAME LigPrep states analysis 05 gave Glide, so ligand preparation is
                held constant and differences are attributable to the docking pipeline
  * site      - defined from the crystal ligand, as Glide's grid was
  * RMSD      - spyrmsd, identical code path to analysis 05

Per structure:
  1. obabel converts the SDFs to V2000 - Schrodinger writes V3000 MOL blocks, which ICM
     cannot read ("no atoms in the object")
  2. `_dockBatch <rec>.pdb ligmol=<site_lig>.sdf proj=D_<PDB>` builds the receptor and maps.
     _dockBatch runs convertObject with hydrogen optimisation, His optimisation and
     hetero charge assignment - the last is what charges the heme, and is exactly what the
     PXR-era invocation (`convertObject a_rec. no no no ...`) disabled.
  3. `_dockScan ... confs=10 effort=3 seed=<N>` docks (10 to match Glide's POSES_PER_LIG,
     so best-of-N is comparable; _dockBatch alone hardcodes confs=3). Called separately from _dockBatch
     purely so `seed=` can be passed: ICM docking is Monte Carlo and therefore stochastic,
     unlike Glide which was verified deterministic. A fixed seed makes runs reproducible.
  4. `_dockScan output=hits.sdf` writes one molecule per pose with the scores as SD tags,
     so each score is bound to its own pose. RMSD is then computed with spyrmsd.

Note on ranking. Scores must be read from `hits.sdf`, NOT parsed from `_dockScan`'s stdout
`SCORES>` lines: those appear in a third ordering and indexing them positionally silently
mis-assigns every score. Verified on 4D6Z — the stdout mapping gave RMSD 2.81 A a score of
-13.94 where the true value is -12.85. `hits.sdf` is ordered worst-to-best by `Score` and
is the exact reverse of the .ob conformation stack.

"Top-1" means the **best-scoring** pose, matching how Glide's poseviewer output is ordered
and what one would select prospectively.

Usage (via scripts/run_icm.sh):
    $SCHRODINGER/run python3 scripts/02_icm_dock.py [--only 4D6Z] [--confs 3] [--seed 42]

Outputs (git-ignored): structures/<PDB>/
Outputs (committed):
    results/icm_poses.csv    per-pose RMSD and ICM scores
    results/icm_summary.csv  per-structure top-1 / best-of-N RMSD
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from schrodinger import structure
from spyrmsd import rmsd as spy

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A05 = os.path.join(ANALYSES, "05_full_pdb_redocking")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")

OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
ICM = "icm64"
ICM_HOME = "/opt/applications/icm/3.9.4"
DOCKBATCH = os.path.join(ICM_HOME, "_dockBatch")
DOCKSCAN = os.path.join(ICM_HOME, "_dockScan")


def require_compute_node(what):
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(f"REFUSING to run {what} outside SLURM. Use: sbatch scripts/run_icm.sh")


def append_csv(path, rows):
    """Append rows, aligning to the existing header.

    Rows written for failures carry fewer keys than success rows. Appending them raw put
    the error text under `resolution` and left `status` empty, so a `status != "ok"`
    filter silently returned nothing and 22 failures looked like successes. Reindexing to
    the existing columns (and refusing to silently drop new ones) prevents that.
    """
    if not rows:
        return
    df = pd.DataFrame(rows)
    if os.path.exists(path):
        cols = list(pd.read_csv(path, nrows=0).columns)
        extra = [c for c in df.columns if c not in cols]
        if extra:
            # New fields would be lost by reindexing - rewrite the file with the union.
            old = pd.read_csv(path)
            df = pd.concat([old, df], ignore_index=True)
            df.to_csv(path, index=False, float_format="%.4f")
            return
        df = df.reindex(columns=cols)
    df.to_csv(path, mode="a", header=not os.path.exists(path), index=False,
              float_format="%.4f")


def to_v2000(src, dst, wdir):
    """ICM cannot read Schrodinger's V3000 MOL blocks; obabel rewrites them as V2000."""
    if os.path.exists(dst):
        return dst
    subprocess.run([OBABEL, os.path.basename(src), "-osdf", "-O", os.path.basename(dst)],
                   cwd=wdir, check=True, capture_output=True)
    if not os.path.exists(dst):
        raise RuntimeError("obabel produced no V2000 output")
    return dst


def graph_arrays(st):
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


def icm_run(args, wdir, script=None):
    cmd = [ICM] + (["-s", script] if script else args)
    return subprocess.run(cmd, cwd=wdir, capture_output=True, text=True, timeout=7200)


def setup_project(pdb, wdir, proj):
    if os.path.exists(os.path.join(wdir, f"{proj}.dtb")):
        return
    require_compute_node("ICM receptor setup")
    r = icm_run([DOCKBATCH, f"{pdb}_rec.pdb", f"ligmol={pdb}_site_lig_v2000.sdf",
                 f"proj={proj}"], wdir)
    if not os.path.exists(os.path.join(wdir, f"{proj}.dtb")):
        raise RuntimeError(f"ICM project not created: {r.stdout[-400:]}")


def dock(pdb, wdir, proj, confs, seed, effort):
    """Dock and return the path to the hitlist SDF (one molecule per pose, scores as SD tags)."""
    lig = f"{pdb}_ligand_v2000"
    hits = os.path.join(wdir, f"{pdb}_hits.sdf")
    if os.path.exists(hits):
        return hits
    require_compute_node("ICM docking")
    r = icm_run([DOCKSCAN, proj, f"input={lig}.sdf", f"confs={confs}",
                 f"effort={effort}", f"seed={seed}",
                 f"output={os.path.basename(hits)}", "-s", "-a", "-S"], wdir)
    if not os.path.exists(hits):
        raise RuntimeError(f"no ICM hitlist produced: {r.stdout[-400:]}")
    return hits


def read_hits(hits, ref_sdf):
    """(rmsds, scores) for every pose in the hitlist, in file order."""
    rx, rz, ra = graph_arrays(structure.StructureReader.read(ref_sdf))
    rmsds, scores = [], []
    for st in structure.StructureReader(hits):
        px, pz, pa = graph_arrays(st)
        try:
            rmsds.append(float(spy.symmrmsd(rx, px, rz, pz, ra, pa, minimize=False)))
        except Exception as e:
            print(f"      spyrmsd failed: {e}", flush=True)
            rmsds.append(np.nan)
        scores.append({k.replace("r_sd_", "").replace("s_sd_", ""): v
                       for k, v in st.property.items()
                       if k.startswith(("r_sd_", "s_sd_"))
                       and k.split("_sd_")[1] in ("Score", "mfScore", "RTCNNscore",
                                                  "Hbond", "Hphob", "Dsolv", "Eintl",
                                                  "Natom", "Nflex")})
    return rmsds, scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    ap.add_argument("--confs", type=int, default=10)  # matches Glide POSES_PER_LIG for a comparable best-of-N
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--effort", default="3.")
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()

    inputs = pd.read_csv(os.path.join(RESULTS, "icm_inputs.csv"))
    inputs = inputs[inputs["status"] == "ok"]
    if args.only:
        inputs = inputs[inputs["pdb_id"].isin(args.only.split(","))]
    os.makedirs(RESULTS, exist_ok=True)

    pose_csv = os.path.join(RESULTS, "icm_poses.csv")
    summ_csv = os.path.join(RESULTS, "icm_summary.csv")
    done = set(pd.read_csv(summ_csv)["pdb_id"]) if os.path.exists(summ_csv) else set()
    todo = [r for _, r in inputs.iterrows() if r["pdb_id"] not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(todo)} structures to dock with ICM "
          f"(confs={args.confs}, effort={args.effort}, seed={args.seed})\n")

    for i, row in enumerate(todo, 1):
        pdb, comp = row["pdb_id"], row["comp_id"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        proj = f"D_{pdb}"
        print(f"[{i}/{len(todo)}] {pdb} ({comp})", flush=True)
        try:
            for base in ("site_lig", "ligand"):
                to_v2000(os.path.join(wdir, f"{pdb}_{base}.sdf"),
                         os.path.join(wdir, f"{pdb}_{base}_v2000.sdf"), wdir)
            setup_project(pdb, wdir, proj)
            hits = dock(pdb, wdir, proj, args.confs, args.seed, args.effort)
            rmsds, scores = read_hits(hits, os.path.join(wdir, f"{pdb}_crystal.sdf"))
            if not rmsds:
                raise RuntimeError("no poses in hitlist")

            sc = [float(d.get("Score")) if d.get("Score") is not None else np.inf
                  for d in scores]
            top = int(np.argmin(sc))            # best (most negative) ICM Score
            best = int(np.nanargmin(rmsds))     # lowest RMSD available

            append_csv(pose_csv, [
                {"pdb_id": pdb, "comp_id": comp, "pose": k + 1, "rmsd": rmsds[k],
                 "is_top_scoring": (k == top), **scores[k]}
                for k in range(len(rmsds))])
            append_csv(summ_csv, [{
                "pdb_id": pdb, "comp_id": comp, "n_poses": len(rmsds),
                "rmsd_top1": rmsds[top], "rmsd_best": rmsds[best],
                "pose_of_top_scoring": top + 1, "pose_of_best": best + 1,
                "score_top1": sc[top] if np.isfinite(sc[top]) else None,
                "confs": args.confs, "seed": args.seed, "effort": args.effort,
                "status": "ok"}])
            print(f"    top1={rmsds[top]:.2f} (best-scoring, pose {top + 1}) "
                  f"best={rmsds[best]:.2f} (pose {best + 1}/{len(rmsds)})", flush=True)
        except Exception as e:
            append_csv(summ_csv, [{"pdb_id": pdb, "comp_id": comp,
                                   "status": f"error: {e}"}])
            print(f"    ERROR: {e}", flush=True)

    if os.path.exists(summ_csv):
        d = pd.read_csv(summ_csv)
        ok = d[d["status"] == "ok"]
        print(f"\n=== ICM: {len(ok)}/{len(d)} succeeded ===")
        for label, col in (("top-1", "rmsd_top1"), ("best-of-N", "rmsd_best")):
            for cut in (2.0, 1.0):
                n = (ok[col] <= cut).sum()
                print(f"  {label:<10} RMSD <= {cut:.1f} A: {n:3d}/{len(ok)} "
                      f"({100*n/len(ok):.0f} %)")
            print(f"  {label:<10} median: {ok[col].median():.2f} A")


if __name__ == "__main__":
    main()
