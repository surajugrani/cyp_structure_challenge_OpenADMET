"""
Analysis 08: score the CRYSTAL poses with ICM (Score + RTCNN), no docking.

This is the ICM counterpart to analysis 01, which scored crystal poses in place with Glide
and found no usable correlation with pIC50 (rho = -0.27). The question here is whether
ICM's scores - in particular the neural-network RTCNN rescoring - do better on the same
poses.

Method. `_dockScan -n` "generates and scores the near-native conformation as provided in
input", i.e. it scores the supplied pose through the normal VLS pipeline, which keeps the
values on the same footing as analysis 06's docked poses.

CAVEAT: it minimises first. The scored pose is 0.28-1.42 A from the deposited coordinates
(n=22, median 0.86 A) - near-native, but NOT the crystal pose. The like-for-like Glide
comparison is therefore `mininplace`, not `inplace`.

The alternative route, `scanScoreExternal`, also works (see scripts/03_sse_rigid7.sh; an
earlier version of this file wrongly called it broken - the real problem was passing the
receptor and ligand as two separate objects instead of merging them and passing molecule
selections). It minimises too, and ICM's docs note its scores are not on the same scale as
VLS docking scores, so the two routes are not interchangeable.

Inputs are reused from analysis 06 (ICM-prepared receptor project, crystal ligand SDF),
so nothing is re-prepared.

Usage (via scripts/run_score_only.sh - compute node only):
    $SCHRODINGER/run python3 scripts/01_score_crystal_poses.py

Outputs:
    results/icm_crystal_scores.csv   Score, RTCNN and other terms for each crystal pose
"""

import os
import subprocess

import pandas as pd
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
A06 = os.path.join(ANALYSES, "06_icm_redocking")
RESULTS = os.path.join(HERE, "results")

OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
ICM_HOME = "/opt/applications/icm/3.9.4"
DOCKSCAN = os.path.join(ICM_HOME, "_dockScan")

KEEP_TAGS = ("Score", "RTCNN", "mfScore", "Hbond", "Hphob", "Dsolv", "Eintl",
             "SolEl", "VwInt", "Natom", "Nflex")


def require_compute_node(what):
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(f"REFUSING to run {what} outside SLURM. "
                     "Use: sbatch scripts/run_score_only.sh")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    cur = pd.read_csv(os.path.join(A01, "data", "curated_set.csv"))

    rows = []
    for n, row in enumerate(cur.itertuples(), 1):
        pdb, comp = row.pdb_id, row.comp_id
        wdir = os.path.join(A06, "structures", pdb)
        proj = f"D_{pdb}"
        rec = {"pdb_id": pdb, "comp_id": comp, "pIC50_median": row.pIC50_median,
               "in_main_set": row.in_main_set}
        print(f"[{n}/{len(cur)}] {pdb} ({comp})", flush=True)
        try:
            if not os.path.exists(os.path.join(wdir, f"{proj}.dtb")):
                raise RuntimeError("no ICM project - run analysis 06 first")

            # ICM needs V2000; Schrodinger writes V3000
            v2000 = os.path.join(wdir, f"{pdb}_crystal_v2000.sdf")
            if not os.path.exists(v2000):
                subprocess.run([OBABEL, f"{pdb}_crystal.sdf", "-osdf",
                                "-O", os.path.basename(v2000)],
                               cwd=wdir, check=True, capture_output=True)

            hits = os.path.join(wdir, f"{pdb}_native_hits.sdf")
            if not os.path.exists(hits):
                require_compute_node("ICM scoring")
                subprocess.run(["icm64", DOCKSCAN, proj,
                                f"input={os.path.basename(v2000)}",
                                "confs=1", "effort=1.", "seed=42",
                                f"output={os.path.basename(hits)}",
                                "-n", "-s", "-a", "-S"],
                               cwd=wdir, capture_output=True, text=True, timeout=3600)
            if not os.path.exists(hits):
                raise RuntimeError("no native-scoring output produced")

            st = next(iter(structure.StructureReader(hits)))
            for k, v in st.property.items():
                if k.startswith(("r_sd_", "s_sd_")) and k.split("_sd_")[1] in KEEP_TAGS:
                    try:
                        rec[k.split("_sd_")[1]] = float(v)
                    except (TypeError, ValueError):
                        pass
            rec["status"] = "ok"
            print(f"    Score={rec.get('Score')} RTCNN={rec.get('RTCNN')}", flush=True)
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"    ERROR: {e}", flush=True)
        rows.append(rec)

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "icm_crystal_scores.csv")
    df.to_csv(out, index=False, float_format="%.4f")
    ok = df[df.status == "ok"]
    print(f"\n{len(ok)}/{len(df)} scored. Wrote {out}")
    if len(ok) >= 4:
        from scipy import stats
        for col in ("Score", "RTCNN", "mfScore"):
            sub = ok.dropna(subset=[col]) if col in ok.columns else ok.iloc[0:0]
            sub = sub[sub.in_main_set]
            if len(sub) >= 4:
                r = stats.spearmanr(sub[col], sub.pIC50_median)
                print(f"  {col:<8} vs pIC50: n={len(sub)} rho={r.statistic:+.3f} "
                      f"p={r.pvalue:.3f}")


if __name__ == "__main__":
    main()
