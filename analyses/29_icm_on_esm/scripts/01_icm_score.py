"""Analysis 29 - ICM-score the three ESMFold2 submissions, to select best of three.

RTCNN is the one scoring function besides Glide that is comparable ACROSS programs: a
co-folding model's own confidence (iPTM, pLDDT) is calibrated to that model and means
nothing between models, whereas RTCNN is computed the same way whatever produced the
pose. It has only ever been run WITHIN AF3's ensemble (analysis 20, submission 14),
where it failed like every other within-ensemble selector - 0.3949 against PL-iPTM's
0.4291 - so its cross-run behaviour is untested.

The three ESMFold2 sets give a clean place to test it, and all three are already scored
by the organizers, so the baseline is known rather than assumed:

    23  default endpoint + our ColabFold MSA     0.4527   best of the project
    26  max_accuracy endpoint, its own MSA       0.4365
    27  default endpoint, no MSA                 0.4384

Selecting the best RTCNN of three per ligand is compared against those. Analysis 28 is
running the same comparison with Glide in both in-place modes, so the two selectors can
be read side by side on identical inputs.

EXPECT LITTLE. Analysis 27 measured a median pocket-superposed ligand RMSD of 1.53 A
between sets 23 and 27, which is the within-ensemble regime where every selector has
failed, not the ~4 A cross-model regime where submission 19 briefly appeared to work -
and submission 24 then showed that apparent success did not survive a third model.
A selection that beats 0.4527 would be a genuine surprise.

Input is each submission's PrepWizard output, work/*_pw.pdb, so all three sets reach ICM
by the same path they reached the leaderboard.

Method, per structure, unchanged from analysis 20:
  * split the prepped complex into receptor (protein + haem) and ligand
  * `_dockBatch` builds an ICM project on THAT structure's own receptor conformation
  * `_dockScan -n` scores the supplied pose through the normal VLS pipeline

CAVEAT from analysis 08: `-n` minimises before scoring, so the scored pose sits
~0.3-1.4 A from the input (median 0.86 A on the crystal set). Near-native scores, not
strictly in-place ones - the same footing as analyses 09 and 20, so comparisons hold.

ONE COMPLEX AT A TIME - ICM licence limit. Resumable: anything already in the CSV is
skipped. Compute node only.
"""
import csv, glob, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
WORK = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")
OUT_CSV = os.path.join(RESULTS, "icm_esm_scores.csv")

OBABEL = os.path.expanduser("~/micromamba/envs/obabel_env/bin/obabel")
ICM_HOME = "/opt/applications/icm/3.9.4"
DOCKBATCH = os.path.join(ICM_HOME, "_dockBatch")
DOCKSCAN = os.path.join(ICM_HOME, "_dockScan")
# analysis 06 looked for "RTCNNscore"; ICM writes "RTCNN".
KEEP = ("Score", "RTCNN", "mfScore", "Hbond", "Hphob", "Dsolv", "Eintl",
        "SolEl", "VwInt", "Natom", "Nflex")
COLS = ["ligand", "set", "status", "fe_dist"] + list(KEEP)


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
        done = {(r["ligand"], r["set"]) for r in csv.DictReader(open(OUT_CSV))}

    files = []
    for setname, subdir in SETS.items():
        for f in sorted(glob.glob(os.path.join(
                ROOT, "submissions", subdir, "work", "*_pw.pdb"))):
            files.append((setname, f))
    print(f"{len(files)} prepped structures across {len(SETS)} sets, "
          f"{len(done)} already scored", flush=True)

    for i, (setname, src) in enumerate(files, 1):
        lig_name = os.path.basename(src)[:-7]
        base = f"{setname}__{lig_name}"
        if (lig_name, setname) in done:
            continue
        d = os.path.join(WORK, base); os.makedirs(d, exist_ok=True)
        row = {"ligand": lig_name, "set": setname}
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
