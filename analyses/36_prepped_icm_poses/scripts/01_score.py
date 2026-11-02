"""Re-prepare each ICM pose in submission 29's receptor, then score it with Glide.

Why this exists. Analysis 34 scored analysis 33's 200 ICM poses with Glide `inplace`
and 177 came back at the 10000 rejection sentinel, van der Waals energies around 1e5.
Analysis 35 traced the cause: heavy-atom contacts are normal at 2.2-3.1 A, but the
closest ALL-atom contact is 0.73-1.6 A and is H...H in nine poses of ten. The receptor's
hydrogen network had been optimised by PrepWizard with the CO-FOLDED ligand in place,
so against an ICM pose the receptor's donors and the ligand's donors point into the same
space. Analysis 35 also confirmed ICM does modify the receptor it docks into - truncated
to 2969 of 4076 heavy atoms, 128 atoms moved more than 1 A - so the pose and the
receptor it would ship in were never relaxed against each other at all.

The fix is to let the receptor relax around each pose rather than asking Glide to score
a mismatch: rebuild the complex, run PrepWizard at `-rmsd 0.6`, regenerate the grid from
that receptor, then score `inplace`.

Why 0.6 and not the project default 0.3. The mismatch here is larger than the one
PrepWizard normally absorbs, because the two structures come from different programs.
0.6 gives the pocket room to accommodate the pose.

Why the ligand moving is acceptable now. `-rmsd` restrains heavy atoms of the whole
complex, so at 0.6 the ICM ligand can also move up to ~0.6 A. That is exactly the
objection that ruled out `mininplace` - the score would describe a pose Glide moved.
Submission 34 settled it: every ligand in submission 29 was displaced a median 0.497 A
and LDDT-PLI changed by 0.0001, with LDDT-LP identical to four decimals. Sub-angstrom
ligand movement is invisible to the metric, so it can no longer corrupt a selection.

Why the grid is rebuilt per pose. A Glide grid encodes the receptor including its
hydrogens. Re-optimising the hydrogen network changes nothing unless the grid is
regenerated from the relaxed receptor, which is the bulk of the cost here.

THE RISK THIS CARRIES, stated before the result. Relaxing the receptor by up to 0.6 A
reshapes the pocket, and LDDT-LP is a scored metric - submission 29 sits at 0.8001.
Any submission built from these poses ships the relaxed receptor along with the pose,
so the pocket score is exposed in a way none of the copy-only submissions were.

NO QUALITY CLAIM is made by any gscore here. Glide ranks candidates within one
receptor; whether that tracks LDDT-PLI is for the leaderboard to say.

Receptors come from analysis 28's prepped complexes, the same ones submission 29's
scores were computed in, with the co-folded ligand deleted and the ICM pose merged in
with explicit bonds - no bond perception from a bare PDB.
"""
import argparse, csv, glob, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(HERE, "work")
RES = os.path.join(HERE, "results")
SEL = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                   "results", "selection.csv")
A28 = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "work")
A34 = os.path.join(ROOT, "analyses", "34_inplace_on_icm_poses", "work")
SCHRODINGER = os.environ.get("SCHRODINGER", "/opt/applications/schrodinger/2025-3")
GLIDE = os.path.join(SCHRODINGER, "glide")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")
PREP_RMSD = "0.6"

SCORE_IN = """GRIDFILE   ../grid.zip
LIGANDFILE   ../lig.mae
DOCKING_METHOD   inplace
FORCEFIELD   OPLS_2005
POSTDOCK   True
COMPRESS_POSES   FALSE
POSE_OUTTYPE   poseviewer
POSES_PER_LIG   1
"""
TERMS = ["gscore", "energy", "evdw", "ecoul", "lipo", "hbond", "rewards",
         "erotb", "esite", "metal"]
COLS = ["ligand", "picked_set", "pose_rank", "status", "fe_lig_bonds_broken",
        "lig_shift_during_prep"] + TERMS


def require_compute_node():
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def one_pose(lig, pset, rank, pose_st, structure, pd):
    import math
    d = os.path.join(WORK, lig, f"pose{rank:02d}")
    os.makedirs(d, exist_ok=True)
    row = {"ligand": lig, "picked_set": pset, "pose_rank": rank}

    raw = os.path.join(d, "raw.mae")
    if not os.path.exists(raw):
        cx = next(iter(structure.StructureReader(
            os.path.join(A28, pset, lig, "complex.mae"))))
        rec = cx.copy()
        rec.deleteAtoms([a.index for a in rec.atom if a.pdbres.strip() == "LIG"])
        merged = rec.merge(pose_st)
        for a in merged.atom:
            if a.pdbres.strip() not in ("HEM",) and a.index > rec.atom_total:
                a.pdbres = "LIG "
        merged.write(raw)

    prepped = os.path.join(d, "prepped.mae")
    if not os.path.exists(prepped):
        subprocess.run([PREPWIZARD, "-NOJOBID", "-rehtreat", "-captermini",
                        "-propka_pH", "7.4", "-epik_pH", "7.4", "-f", "OPLS_2005",
                        "-rmsd", PREP_RMSD, "raw.mae", "prepped.mae"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=5400)
    if not os.path.exists(prepped):
        row["status"] = "error: prepwizard produced no .mae"
        return row

    st = next(iter(structure.StructureReader(prepped)))
    lig_idx = [a.index for a in st.atom if a.pdbres.strip() == "LIG"]
    if not lig_idx:
        row["status"] = "error: no LIG after prepwizard"
        return row

    # how far PrepWizard moved the ligand heavy atoms, so the cost is on record
    before = [a.xyz for a in pose_st.atom if a.atomic_number > 1]
    after = [st.atom[i].xyz for i in lig_idx if st.atom[i].atomic_number > 1]
    if len(before) == len(after):
        row["lig_shift_during_prep"] = round(
            (sum(math.dist(a, b) ** 2 for a, b in zip(before, after))
             / len(before)) ** 0.5, 3)

    # PrepWizard bonds a coordinating atom to the iron; Glide cannot score a
    # covalently attached ligand in place, so the bond is broken first.
    lig_set, nfe = set(lig_idx), 0
    for a in st.atom:
        if a.element != "Fe":
            continue
        for b in list(a.bonded_atoms):
            if b.index in lig_set:
                st.deleteBond(a, b); nfe += 1
    row["fe_lig_bonds_broken"] = nfe

    ligf = os.path.join(d, "..", "lig.mae")
    st.extract(lig_idx).write(os.path.join(d, "lig_single.mae"))
    cmae = os.path.join(d, "complex.mae")
    if not os.path.exists(cmae):
        st.write(cmae)

    # Schrodinger 2025-3 names the grid after the job stage, grid-gridgen.zip, not
    # grid.zip. Glob for it the way analysis 28 does rather than assuming the name.
    g_zip = glob.glob(os.path.join(d, "grid*.zip"))
    if not g_zip:
        subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT", "-rec_file", "complex.mae",
                        "-lig_asl", "res. LIG", "-j", "grid"],
                       cwd=d, check=False, capture_output=True, text=True, timeout=5400)
        g_zip = glob.glob(os.path.join(d, "grid*.zip"))
    if not g_zip:
        row["status"] = "error: no grid produced"
        return row
    grid = os.path.join("..", os.path.basename(g_zip[0]))

    sub = os.path.join(d, "score")
    os.makedirs(sub, exist_ok=True)
    csvp = os.path.join(sub, "score.csv")
    if not os.path.exists(csvp):
        with open(os.path.join(sub, "score.in"), "w") as f:
            f.write(SCORE_IN.replace("../grid.zip", grid)
                            .replace("../lig.mae", "../lig_single.mae"))
        subprocess.run([GLIDE, "-NOJOBID", "-WAIT", "score.in"], cwd=sub, check=False,
                       capture_output=True, text=True, timeout=5400)
    if not os.path.exists(csvp):
        row["status"] = "error: no glide results csv"
        return row

    g = pd.read_csv(csvp)
    if g.empty:
        row["status"] = "error: glide returned no pose"
        return row
    r0 = g.iloc[0]
    for t in TERMS:
        c = f"r_i_glide_{t}"
        if c in g.columns:
            try:
                row[t] = float(r0[c])
            except (TypeError, ValueError):
                pass
    row["status"] = "ok"
    return row


def main():
    require_compute_node()
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True); os.makedirs(RES, exist_ok=True)
    from schrodinger import structure
    import pandas as pd

    out = os.path.join(RES, f"scores_shard{a.shard}.csv")
    done = set()
    if os.path.exists(out):
        done = {(r["ligand"], int(r["pose_rank"])) for r in csv.DictReader(open(out))
                if r["status"] == "ok"}

    picks = [(r["ligand"], r["picked_set"]) for r in csv.DictReader(open(SEL))]
    mine = [p for i, p in enumerate(sorted(picks)) if i % a.of == a.shard]
    print(f"shard {a.shard}/{a.of}: {len(mine)} ligands, {len(done)} poses already done",
          flush=True)

    for lig, pset in mine:
        poses = list(structure.StructureReader(
            os.path.join(A34, lig, "poses.mae")))
        for rank, ps in enumerate(poses, 1):
            if (lig, rank) in done:
                continue
            try:
                row = one_pose(lig, pset, rank, ps, structure, pd)
            except Exception as e:
                row = {"ligand": lig, "picked_set": pset, "pose_rank": rank,
                       "status": f"error: {e}"[:200]}
            new = not os.path.exists(out)
            with open(out, "a", newline="") as fh:
                w = csv.DictWriter(fh, COLS, extrasaction="ignore")
                if new:
                    w.writeheader()
                w.writerow(row)
            print(f"  {lig} pose {rank:2}  {row.get('status')}  "
                  f"gscore={row.get('gscore')}  "
                  f"moved={row.get('lig_shift_during_prep')}", flush=True)
    print("\ndone")
    return 0


if __name__ == "__main__":
    sys.exit(main())
