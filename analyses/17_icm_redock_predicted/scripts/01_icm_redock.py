"""Analysis 17 - redock all 20 blinded ligands into BOTH predicted receptors with ICM.

Receptors:
  boltz_03  submission 03 (Boltz-2 + PrepWizard)  - LDDT-LP 0.7962, the better pocket
  af3_06    submission 06 (AF3 + PrepWizard)      - LDDT-LP 0.6985, our best LDDT-PLI

Why both: redocking replaces the ligand but keeps the receptor, so receptor quality is the
variable. Boltz's pocket scores 0.098 higher on LDDT-LP; whether that translates into
better docked poses is the question.

Prospective by construction: the ligand is built in 3D from SMILES (obabel --gen3d), NOT
from the co-folded coordinates, so no pose information leaks in. The co-folded ligand is
used only to define the binding site for the ICM project.

Settings match analysis 06 so the numbers are comparable: confs=10, effort=3., seed=42.

There is no ground truth for the blinded set, so this reports what can be measured:
  * RMSD between the ICM top-1 pose and the co-folded pose (how much it moved)
  * Fe-ligand distance of the ICM pose (analysis 16's type I / type II diagnostic)
  * ICM Score and RTCNN per pose

ONE COMPLEX AT A TIME, serial, in a single job - ICM licence limit.
Resumable: existing per-complex outputs are skipped on resubmission.
"""
import argparse, csv, os, subprocess, sys

import numpy as np
import pandas as pd
from schrodinger import structure
from spyrmsd import rmsd as spy

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")
OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
ICM = "icm64"
ICM_HOME = "/opt/applications/icm/3.9.4"
DOCKBATCH = os.path.join(ICM_HOME, "_dockBatch")
DOCKSCAN = os.path.join(ICM_HOME, "_dockScan")

RECEPTORS = {
    "boltz_03": os.path.join(ROOT, "submissions", "03_boltz2_prepwizard", "submission"),
    "af3_06": os.path.join(ROOT, "submissions", "06_af3_pliptm", "submission"),
}
# NOTE: analysis 06 looked for an SD tag named "RTCNNscore"; ICM actually writes "RTCNN",
# so that column came back empty there. Fixed here.
SD_TAGS = ("Score", "mfScore", "RTCNN", "Hbond", "Hphob", "Dsolv", "Eintl",
           "Natom", "Nflex")


def require_compute_node(what):
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(f"REFUSING to run {what} outside SLURM. Use: sbatch scripts/run.sh")


def icm_run(args, wdir):
    return subprocess.run([ICM] + args, cwd=wdir, capture_output=True, text=True,
                          timeout=7200)


def split_complex(src, wdir, tag):
    """Receptor (protein + HEM, no LIG) and the co-folded ligand, as separate files."""
    rec = os.path.join(wdir, f"{tag}_rec.pdb")
    ligpdb = os.path.join(wdir, f"{tag}_cofolded.pdb")
    if not (os.path.exists(rec) and os.path.exists(ligpdb)):
        rl, ll = [], []
        for l in open(src):
            if l.startswith(("ATOM", "HETATM")):
                (ll if l[17:20].strip() == "LIG" else rl).append(l)
            elif l.startswith(("TER", "END")):
                continue
        open(rec, "w").writelines(rl + ["END\n"])
        open(ligpdb, "w").writelines(ll + ["END\n"])
    return rec, ligpdb


def obabel(args, wdir):
    subprocess.run([OBABEL] + args, cwd=wdir, check=True, capture_output=True)


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


def fe_dist(rec_pdb, pose_xyz):
    fe = None
    for l in open(rec_pdb):
        if l.startswith(("ATOM", "HETATM")) and \
                (l[76:78].strip() or l[12:16].strip()[0]).upper() == "FE":
            fe = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            break
    if fe is None:
        return None
    return float(np.min(np.linalg.norm(pose_xyz - fe, axis=1)))


def run_one(name, smi, recname, confs, seed, effort):
    tag = f"{recname}__{name}"
    wdir = os.path.join(WORK, tag)
    os.makedirs(wdir, exist_ok=True)
    hits = os.path.join(wdir, "hits.sdf")
    src = os.path.join(RECEPTORS[recname], f"{name}.pdb")

    rec, ligpdb = split_complex(src, wdir, tag)
    site = os.path.join(wdir, "site_lig.sdf")
    if not os.path.exists(site):
        obabel([os.path.basename(ligpdb), "-osdf", "-O", "site_lig.sdf"], wdir)
    lig = os.path.join(wdir, "ligand.sdf")
    if not os.path.exists(lig):
        # built from SMILES, so no co-folded pose information reaches the docking
        obabel([f"-:{smi}", "--gen3d", "-osdf", "-O", "ligand.sdf"], wdir)

    proj = "proj"
    if not os.path.exists(os.path.join(wdir, f"{proj}.dtb")):
        require_compute_node("ICM receptor setup")
        r = icm_run([DOCKBATCH, os.path.basename(rec), "ligmol=site_lig.sdf",
                     f"proj={proj}"], wdir)
        if not os.path.exists(os.path.join(wdir, f"{proj}.dtb")):
            raise RuntimeError(f"no ICM project: {r.stdout[-400:]}")

    if not os.path.exists(hits):
        require_compute_node("ICM docking")
        r = icm_run([DOCKSCAN, proj, "input=ligand.sdf", f"confs={confs}",
                     f"effort={effort}", f"seed={seed}", "output=hits.sdf",
                     "-s", "-a", "-S"], wdir)
        if not os.path.exists(hits):
            raise RuntimeError(f"no hitlist: {r.stdout[-400:]}")

    cof = os.path.join(wdir, "cofolded.sdf")
    if not os.path.exists(cof):
        obabel([os.path.basename(ligpdb), "-osdf", "-O", "cofolded.sdf"], wdir)
    cx, cz, ca = graph_arrays(structure.StructureReader.read(cof))

    rows = []
    for i, st in enumerate(structure.StructureReader(hits)):
        px, pz, pa = graph_arrays(st)
        try:
            d = float(spy.symmrmsd(cx, px, cz, pz, ca, pa, minimize=False))
        except Exception:
            d = np.nan
        sc = {k.split("_sd_")[1]: v for k, v in st.property.items()
              if k.startswith(("r_sd_", "s_sd_")) and k.split("_sd_")[1] in SD_TAGS}
        rows.append({"receptor": recname, "ligand": name, "pose_rank": i + 1,
                     "rmsd_to_cofolded": None if np.isnan(d) else round(d, 3),
                     "fe_dist": None if fe_dist(rec, px) is None
                     else round(fe_dist(rec, px), 2), **sc})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--confs", type=int, default=10)
    ap.add_argument("--effort", default="3.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--only")
    ap.add_argument("--limit", type=int)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True); os.makedirs(RESULTS, exist_ok=True)

    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    names = sorted(smiles)
    if a.only:
        names = [n for n in names if n == a.only]
    if a.limit:
        names = names[:a.limit]

    out = os.path.join(RESULTS, "icm_poses.csv")
    done = set()
    if os.path.exists(out):
        d = pd.read_csv(out)
        done = set(zip(d["receptor"], d["ligand"]))

    total = len(names) * len(RECEPTORS)
    k = 0
    for recname in RECEPTORS:
        for n in names:
            k += 1
            if (recname, n) in done:
                print(f"[{k}/{total}] {recname} {n}  already done", flush=True)
                continue
            print(f"[{k}/{total}] {recname} {n}", flush=True)
            try:
                rows = run_one(n, smiles[n], recname, a.confs, a.seed, a.effort)
            except Exception as e:
                print(f"    FAILED: {e}", flush=True)
                rows = [{"receptor": recname, "ligand": n, "pose_rank": 0,
                         "error": str(e)[:200]}]
            df = pd.DataFrame(rows)
            df.to_csv(out, mode="a", header=not os.path.exists(out), index=False,
                      float_format="%.4f")
            if rows and "Score" in rows[0]:
                best = min(rows, key=lambda r: float(r.get("Score", 1e9)))
                print(f"    top-1 Score {best.get('Score')}  "
                      f"rmsd_to_cofolded {best.get('rmsd_to_cofolded')} A  "
                      f"Fe {best.get('fe_dist')} A", flush=True)
    print("\ndone")
    return 0


if __name__ == "__main__":
    sys.exit(main())
