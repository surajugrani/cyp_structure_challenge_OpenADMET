"""
Analysis 09, step 4: generate a `_dockScan -n` near-native pose for every structure in the
30-structure pose-selection set, so each ligand has an 11th candidate that is genuinely
near-native (0.3-1.4 A from crystal) rather than only its 10-30 docked poses.

Why: in 10 of the 30 structures docking never sampled a pose under 2 A (best 2.2-6.4 A),
so those ligands contribute only negative examples to a pose-selector fit. `-n` minimises
the deposited ligand in place and scores it through ICM, giving a positive example for
every ligand.

KNOWN CONFOUND, recorded here so it is not forgotten downstream: `-n` scores are NOT on the
same scale as `_dockScan` docked-pose scores. Measured on the 7 structures where both
already existed, the `-n` route scores the SAME geometry 0.25-13.4 Score units better
(3UA1: 9.3 units better at only 0.66 A separation), and the offset is not constant so it
cannot be subtracted away. A model trained on a pooled set of docked + `-n` rows can
therefore identify the `-n` row from its score alone. Any use of these rows must account
for that.

Serial by design - one structure at a time, no array job, because ICM licences are
limited. Resumable: structures whose `<pdb>_native_hits.sdf` already exists are skipped.

Usage (compute node only):
    sbatch analyses/09_pose_selection/scripts/04_native_poses.sh

Outputs:
    <a06>/structures/<pdb>/<pdb>_native_hits.sdf   pose + scores, one per structure
    results/native_pose_status.csv                 per-structure status
"""

import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A06 = os.path.join(ANALYSES, "06_icm_redocking")
RESULTS = os.path.join(HERE, "results")

OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
DOCKSCAN = "/opt/applications/icm/3.9.4/_dockScan"
TIMEOUT = 1800


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run ICM outside SLURM. "
                     "Use: sbatch analyses/09_pose_selection/scripts/04_native_poses.sh")


def main():
    require_compute_node()
    os.makedirs(RESULTS, exist_ok=True)

    poses = pd.read_csv(os.path.join(RESULTS, "rtcnn_rescore_poses.csv"))
    pdbs = sorted(poses.pdb_id.unique())

    rows, done, failed = [], 0, 0
    for k, pdb in enumerate(pdbs, 1):
        wdir = os.path.join(A06, "structures", pdb)
        hits = os.path.join(wdir, f"{pdb}_native_hits.sdf")
        rec = {"pdb_id": pdb}
        if os.path.exists(hits):
            rec["status"] = "already_present"
            rows.append(rec)
            print(f"[{k}/{len(pdbs)}] {pdb}: already present", flush=True)
            continue
        print(f"[{k}/{len(pdbs)}] {pdb}: running _dockScan -n", flush=True)
        try:
            proj = f"D_{pdb}"
            if not os.path.exists(os.path.join(wdir, f"{proj}.dtb")):
                raise RuntimeError("no ICM project (.dtb)")

            # ICM reads V2000 only; Schrodinger writes V3000
            v2000 = os.path.join(wdir, f"{pdb}_crystal_v2000.sdf")
            if not os.path.exists(v2000):
                subprocess.run([OBABEL, f"{pdb}_crystal.sdf", "-osdf",
                                "-O", os.path.basename(v2000)],
                               cwd=wdir, check=True, capture_output=True)

            r = subprocess.run(
                ["icm64", DOCKSCAN, proj,
                 f"input={os.path.basename(v2000)}",
                 "confs=1", "effort=1.", "seed=42",
                 f"output={os.path.basename(hits)}",
                 "-n", "-s", "-a", "-S"],
                cwd=wdir, capture_output=True, text=True, timeout=TIMEOUT)
            if not os.path.exists(hits):
                raise RuntimeError(f"no output (rc={r.returncode}); "
                                   f"stderr tail: {r.stderr[-300:]}")
            rec["status"] = "ok"
            done += 1
        except Exception as e:
            rec["status"] = f"error: {e}"
            failed += 1
            print(f"    ERROR: {e}", flush=True)
        rows.append(rec)
        # write after every structure so progress survives a walltime kill
        pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "native_pose_status.csv"),
                                  index=False)

    print(f"\nnewly scored {done}, failed {failed}, "
          f"total present {sum(1 for r in rows if not str(r['status']).startswith('error'))}"
          f"/{len(pdbs)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
