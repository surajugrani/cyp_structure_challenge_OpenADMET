"""
Analysis 09, step 10: rescore ICM's poses with Glide `inplace`, and ask how many
structures Glide's Emodel recovers.

Motivation. Step 9b found Emodel is the strongest-looking Glide term (within-ligand
rho = +0.378, 25/30 top-1) but that all of it is route detection: Emodel is 13.4 units
better on `mininplace` crystal rows than on docked rows, against a within-ligand spread of
5.2, and on docked rows alone Emodel scores 9/30 - exactly the gscore baseline, saving 2
and breaking 2.

This design removes that artefact by construction. Every candidate - all 10 ICM docked
poses AND the ICM `-n` near-native pose - is scored through the SAME Glide `inplace` call,
so there is no mixed provenance for Emodel to detect. Any Emodel advantage that survives
is pose discrimination.

`mininplace`, not `inplace`. An earlier version of this script used `inplace` to keep the
RMSD labels from analyses 06 and 09 exact. That does not work, for two reasons found by
running it:

  1. **Emodel is identically 0.00 under `inplace`** (as is einternal) - Glide only computes
     it during the minimisation/post-docking stage, which `inplace` skips. `POSTDOCK True`
     does not change this. The one term this analysis exists to test was absent.
  2. ICM poses were optimised in an ICM-prepared receptor and hard-clash in the
     PrepWizard-prepared one: gscore hit +10000 and evdw +25800 on some rows, Glide's
     clash sentinel. Those rows are unusable rather than merely penalised.

`mininplace` fixes both - it computes Emodel and relaxes the clashes. The RMSD labels are
then **recomputed** from the minimised pose written to `<job>_pv.mae`, using the same
symmetry-corrected spyrmsd measurement used for the `-n` poses, rather than carried over.
Labels therefore describe the geometry actually scored, which is the correct target.
The artefact control is unaffected: all 11 rows of a ligand go through one `mininplace`
call, so there is no provenance difference for Emodel to detect.

Coordinate frames were checked to match (crystal-ligand centroids agree to 0.3-0.9 A
between the two pipelines, i.e. protomer/H differences, not a transform).

Serial, resumable, compute node only.

Usage:
    sbatch analyses/09_pose_selection/scripts/08_glide_rescore_icm.sh

Outputs:
    results/glide_on_icm_poses.csv   one row per ICM pose with Glide terms + recomputed RMSD
"""

import glob
import os
import re
import subprocess
import sys

import pandas as pd
from schrodinger import structure

# .pylibs lives at the repo root: scripts -> 09_pose_selection -> analyses -> root
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))), ".pylibs"))
from spyrmsd import rmsd as srm, io as sio  # noqa: E402


def load_ref(pdb):
    """Deposited ligand, as the RMSD reference. Same file analysis 05 used."""
    f = os.path.join(A05, pdb, f"{pdb}_crystal_lig.sdf")
    return next(iter(structure.StructureReader(f)))


def rmsd_to_ref(st, ref):
    """Symmetry-corrected heavy-atom RMSD, via spyrmsd on Schrodinger structures."""
    try:
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            a, b = os.path.join(td, "a.sdf"), os.path.join(td, "b.sdf")
            with structure.StructureWriter(a) as w:
                w.append(ref)
            with structure.StructureWriter(b) as w:
                w.append(st)
            ra, rb = sio.loadmol(a), sio.loadmol(b)
            ra.strip()
            rb.strip()
            return float(srm.symmrmsd(ra.coordinates, rb.coordinates,
                                      ra.atomicnums, rb.atomicnums,
                                      ra.adjacency_matrix, rb.adjacency_matrix))
    except Exception:
        return None

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A05 = os.path.join(ANALYSES, "05_full_pdb_redocking", "structures")
A06 = os.path.join(ANALYSES, "06_icm_redocking", "structures")
RESULTS = os.path.join(HERE, "results")

SCHRODINGER = os.environ["SCHRODINGER"]
GLIDE = os.path.join(SCHRODINGER, "glide")
MAX_DOCKED = 10

GLIDE_IN = """GRIDFILE   {grid}
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


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Glide outside SLURM. Use: sbatch "
                     "analyses/09_pose_selection/scripts/08_glide_rescore_icm.sh")


def sdf_scores(path):
    """Score tags per entry, in file order."""
    out = []
    for b in open(path).read().split("$$$$"):
        if "<Score>" not in b:
            continue
        out.append({m.group(1): m.group(2).strip()
                    for m in re.finditer(r"> *<([^>]+)>\s*\n([^\n]*)", b)})
    return out


def main():
    require_compute_node()
    os.makedirs(RESULTS, exist_ok=True)
    out_csv = os.path.join(RESULTS, "glide_on_icm_poses.csv")

    pool = pd.read_csv(os.path.join(RESULTS, "icm_pool11_poses.csv"))
    nat = pd.read_csv(os.path.join(RESULTS, "native_pose_rmsd.csv")).set_index("pdb_id")
    pdbs = sorted(pool.pdb_id.unique())
    done = set(pd.read_csv(out_csv).pdb_id) if os.path.exists(out_csv) else set()

    allrows = []
    for k, pdb in enumerate(pdbs, 1):
        if pdb in done:
            print(f"[{k}/{len(pdbs)}] {pdb}: already done", flush=True)
            continue
        print(f"[{k}/{len(pdbs)}] {pdb}: glide mininplace on ICM poses", flush=True)
        try:
            wdir = os.path.join(A05, pdb)
            grid = glob.glob(os.path.join(wdir, f"{pdb}_grid*.zip"))
            if not grid:
                raise RuntimeError("no grid")

            # the 10 best-Score ICM docked poses, in the order icm_pool11 selected them
            sel = pool[(pool.pdb_id == pdb) & (~pool.is_native_row)] \
                .nsmallest(MAX_DOCKED, "Score")
            want = {round(v, 4) for v in sel.Score}

            sts, meta = [], []
            docked = list(structure.StructureReader(
                os.path.join(A06, pdb, f"{pdb}_hits.sdf")))
            tags = sdf_scores(os.path.join(A06, pdb, f"{pdb}_hits.sdf"))
            for i, (st, t) in enumerate(zip(docked, tags)):
                if round(float(t["Score"]), 4) not in want:
                    continue
                rm = sel[sel.Score.round(4) == round(float(t["Score"]), 4)].rmsd
                st.title = f"{pdb}_d{i}"
                sts.append(st)
                meta.append({"pdb_id": pdb, "title": st.title,
                             "rmsd": float(rm.iloc[0]), "is_native_row": False})

            nst = next(iter(structure.StructureReader(
                os.path.join(A06, pdb, f"{pdb}_native_hits.sdf"))))
            nst.title = f"{pdb}_native"
            sts.append(nst)
            meta.append({"pdb_id": pdb, "title": nst.title,
                         "rmsd": float(nat.loc[pdb, "n_rmsd"]), "is_native_row": True})

            ligf = os.path.join(wdir, f"{pdb}_icmposes.sdf")
            with structure.StructureWriter(ligf) as w:
                for st in sts:
                    w.append(st)

            job = f"{pdb}_icmrescore"
            with open(os.path.join(wdir, f"{job}.in"), "w") as f:
                f.write(GLIDE_IN.format(grid=os.path.basename(grid[0]),
                                        lig=os.path.basename(ligf)))
            subprocess.run([GLIDE, "-NOJOBID", "-WAIT", f"{job}.in"],
                           cwd=wdir, check=False, capture_output=True,
                           text=True, timeout=7200)

            csvp = os.path.join(wdir, f"{job}.csv")
            if not os.path.exists(csvp):
                raise RuntimeError("no results csv")
            g = pd.read_csv(csvp)
            m = pd.DataFrame(meta)

            # Recompute each label from the MINIMISED pose, not the ICM input pose.
            pv = os.path.join(wdir, f"{job}_pv.mae")
            newr = {}
            if os.path.exists(pv):
                ref = load_ref(pdb)
                for st in list(structure.StructureReader(pv))[1:]:  # [0] is the receptor
                    v = rmsd_to_ref(st, ref)
                    if v is not None:
                        newr[st.title] = v
            m["rmsd_icm_input"] = m.rmsd
            m["rmsd"] = m.title.map(newr).fillna(m.rmsd)
            m["rmsd_recomputed"] = m.title.isin(newr)
            tcol = "title" if "title" in g.columns else g.columns[0]
            g = g.rename(columns={tcol: "title"})
            for t in TERMS:
                col = f"r_i_glide_{t}"
                if col in g.columns:
                    g[t] = pd.to_numeric(g[col], errors="coerce")
            merged = m.merge(g[["title"] + [t for t in TERMS if t in g.columns]],
                             on="title", how="left")
            miss = int(merged.gscore.isna().sum()) if "gscore" in merged else len(merged)
            print(f"    {len(merged)} poses in, {len(merged) - miss} scored", flush=True)
            allrows.append(merged)
        except Exception as e:
            print(f"    ERROR: {e}", flush=True)
            allrows.append(pd.DataFrame([{"pdb_id": pdb, "title": None,
                                          "status": f"error: {e}"}]))
        d = pd.concat(allrows, ignore_index=True)
        if os.path.exists(out_csv):
            d = pd.concat([pd.read_csv(out_csv), d], ignore_index=True)
            d = d.drop_duplicates(["pdb_id", "title"], keep="last")
        d.to_csv(out_csv, index=False, float_format="%.4f")

    print(f"\nWrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
