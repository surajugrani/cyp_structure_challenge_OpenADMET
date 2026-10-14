"""Analysis 20 - score all 100 PrepWizard-minimised AF3 samples in place with ICM.

Test #3 of the agreed three: can ICM's scores pick the right sample out of AF3's five,
where Glide's could not? Submission 13 answered the Glide half - selecting on gscore
scored 0.3730 against submission 06's 0.4291, a loss of 0.056, confirming prospectively
what analysis 09 measured on crystals (within-ligand Spearman vs RMSD only +0.086).

ICM is the better candidate on the same measurement: RTCNN reached +0.334 in analysis 09,
nearly four times Glide's best term, and analysis 09 chose RTCNN + Hbond as the preferred
ICM model. Whether that transfers from crystal receptors to co-folded ones is the question.

Method, per sample, following analysis 08:
  * split the prepped complex into receptor (protein + haem) and ligand
  * `_dockBatch` builds an ICM project on THAT sample's own receptor conformation - each
    of the 100 is a different conformation, so each needs its own project
  * `_dockScan -n` "generates and scores the near-native conformation as provided in
    input", i.e. scores the supplied pose through the normal VLS pipeline

CAVEAT carried over from analysis 08: `-n` minimises before scoring, so the scored pose
sits ~0.3-1.4 A from the input (median 0.86 A on the crystal set). These are near-native
scores, not strictly in-place ones. That is the same footing as analysis 09's numbers, so
the comparison holds.

ONE COMPLEX AT A TIME - ICM licence limit. Resumable: a sample already in the CSV is
skipped. Compute node only.
"""
import csv, glob, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
PWORK = os.path.join(ROOT, "submissions", "04_af3_default", "glide_rescoring", "work")
WORK = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")
OUT_CSV = os.path.join(RESULTS, "icm_sample_scores.csv")

OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
ICM_HOME = "/opt/applications/icm/3.9.4"
DOCKBATCH = os.path.join(ICM_HOME, "_dockBatch")
DOCKSCAN = os.path.join(ICM_HOME, "_dockScan")
# analysis 06 looked for "RTCNNscore"; ICM writes "RTCNN".
KEEP = ("Score", "RTCNN", "mfScore", "Hbond", "Hphob", "Dsolv", "Eintl",
        "SolEl", "VwInt", "Natom", "Nflex")
COLS = ["ligand", "sample", "status", "fe_dist"] + list(KEEP)


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run ICM outside SLURM. Use sbatch.")


def split(src, d):
    rec, lig = os.path.join(d, "rec.pdb"), os.path.join(d, "lig.pdb")
    if not (os.path.exists(rec) and os.path.exists(lig)):
        rl, ll = [], []
        for l in open(src):
            if l.startswith(("ATOM", "HETATM")):
                (ll if l[17:20].strip() == "LIG" else rl).append(l)
        open(rec, "w").writelines(rl + ["END\n"])
        open(lig, "w").writelines(ll + ["END\n"])
    return rec, lig


def fe_distance(rec, lig):
    import math
    fe = None
    for l in open(rec):
        if l.startswith(("ATOM", "HETATM")) and \
                (l[76:78].strip() or l[12:16].strip()[0]).upper() == "FE":
            fe = (float(l[30:38]), float(l[38:46]), float(l[46:54])); break
    if fe is None:
        return None
    best = None
    for l in open(lig):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if (l[76:78].strip() or l[12:16].strip()[0]).upper() == "H":
            continue
        d = math.dist(fe, (float(l[30:38]), float(l[38:46]), float(l[46:54])))
        best = d if best is None else min(best, d)
    return None if best is None else round(best, 2)


def main():
    require_compute_node()
    os.makedirs(WORK, exist_ok=True); os.makedirs(RESULTS, exist_ok=True)
    from schrodinger import structure

    done = set()
    if os.path.exists(OUT_CSV):
        done = {(r["ligand"], r["sample"]) for r in csv.DictReader(open(OUT_CSV))}

    files = sorted(glob.glob(os.path.join(PWORK, "*", "*_pw.pdb")))
    print(f"{len(files)} prepped samples, {len(done)} already scored", flush=True)

    for i, src in enumerate(files, 1):
        base = os.path.basename(src)[:-7]
        lig_name, tag = base.split("__")
        if (lig_name, tag) in done:
            continue
        d = os.path.join(WORK, base); os.makedirs(d, exist_ok=True)
        row = {"ligand": lig_name, "sample": tag}
        print(f"[{i}/{len(files)}] {base}", flush=True)
        try:
            rec, ligpdb = split(src, d)
            row["fe_dist"] = fe_distance(rec, ligpdb)
            v2000 = os.path.join(d, "lig_v2000.sdf")
            if not os.path.exists(v2000):
                subprocess.run([OBABEL, "lig.pdb", "-osdf", "-O", "lig_v2000.sdf"],
                               cwd=d, check=True, capture_output=True)
            proj = "P"
            if not os.path.exists(os.path.join(d, f"{proj}.dtb")):
                subprocess.run(["icm64", DOCKBATCH, "rec.pdb",
                                "ligmol=lig_v2000.sdf", f"proj={proj}"],
                               cwd=d, capture_output=True, text=True, timeout=3600)
            if not os.path.exists(os.path.join(d, f"{proj}.dtb")):
                raise RuntimeError("no ICM project created")
            hits = os.path.join(d, "native_hits.sdf")
            if not os.path.exists(hits):
                subprocess.run(["icm64", DOCKSCAN, proj, "input=lig_v2000.sdf",
                                "confs=1", "effort=1.", "seed=42",
                                "output=native_hits.sdf", "-n", "-s", "-a", "-S"],
                               cwd=d, capture_output=True, text=True, timeout=3600)
            if not os.path.exists(hits):
                raise RuntimeError("no native-scoring output")
            st = next(iter(structure.StructureReader(hits)))
            for k, v in st.property.items():
                if k.startswith(("r_sd_", "s_sd_")) and k.split("_sd_")[1] in KEEP:
                    try:
                        row[k.split("_sd_")[1]] = float(v)
                    except (TypeError, ValueError):
                        pass
            row["status"] = "ok"
            print(f"    Score={row.get('Score')} RTCNN={row.get('RTCNN')} "
                  f"Hbond={row.get('Hbond')}", flush=True)
        except Exception as e:
            row["status"] = f"error: {e}"[:200]
            print(f"    ERROR: {e}", flush=True)
        new = not os.path.exists(OUT_CSV)
        with open(OUT_CSV, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerow(row)
    print(f"\nwrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
