"""Analysis 12: does a Glide metal constraint fix Fe coordination in redocking?

Analysis 04 found Glide docking systematically pushes ligands off the haem iron (median
+0.73 A vs crystal, 20/22 structures, p = 0.0003) and that its metal reward only fires
below ~2.6 A, so only 13 of 214 docked poses ever earned any. Glide supports a metal
coordination constraint, which turns that from a weak reward into a hard filter.

Design - 10 structures, 2 arms, so the control is built in:

  5 that DO coordinate in the crystal (8SPD 1.94, 4I4G 2.10, 6MA7 2.20, 6DAJ 2.31,
    6UNI 2.40 A) - the constraint should help here
  5 that DO NOT (6UNJ 2.82, 3UA1 3.65, 7LXL 3.73, 6OO9 3.76, 4D6Z 6.09 A) - the
    constraint should HURT here

If it only ever helps, it is probably just dragging ligands to a site that happens to be
right. If it helps the first five and hurts the last five, it is doing something specific.

The unconstrained arm is NOT re-run: analysis 05 already docked all ten from SMILES with
the same grids and LigPrep inputs, so those RMSDs are reused directly.

Constraint setup: the grid is rebuilt with `-metal_cons "fe#<Fe atom number>"` and docking
adds USE_CONS / NREQUIRED_CONS so a pose must satisfy it.
"""
import csv, glob, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A05 = os.path.join(ROOT, "analyses", "05_full_pdb_redocking", "structures")
WORK = os.path.join(HERE, "work")
RESULTS = os.path.join(HERE, "results")
SCHRODINGER = os.environ["SCHRODINGER"]
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")

COORD = ["8SPD", "4I4G", "6MA7", "6DAJ", "6UNI"]
NONCOORD = ["6UNJ", "3UA1", "7LXL", "6OO9", "4D6Z"]

# Docked through XGlide, not glide directly.
#
# A metal constraint needs both a CONSTRAINT_GROUP referencing the grid's constraint
# label AND a matching FEATURE section defining the ligand feature that may satisfy it
# (default "Acceptor" for a metal constraint). Passing `USE_CONS fe:1` to glide without
# that section fails with "feature '1' is not defined", and passing bare `USE_CONS fe`
# silently returns zero poses. Schrodinger's own docs say the feature file "is not
# designed to be written from scratch".
#
# XGlide's DOCK_USE_CONS fills in the default feature for the constraint type, which is
# the supported way to do this without hand-writing SMARTS.
DOCK_IN = """GRID   {grid}
LIGAND   {lig}
DOCK_PRECISION   SP
DOCK_USE_CONS   fe
DOCK_NREQUIRED_CONS   1
DOCK_POSES_PER_LIG   10
"""


def main():
    if not (os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1"):
        raise SystemExit("REFUSING to run Schrodinger outside SLURM.")
    from schrodinger import structure
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    out_csv = os.path.join(RESULTS, "constrained_docking.csv")
    done = set()
    if os.path.exists(out_csv):
        done = {r["pdb_id"] for r in csv.DictReader(open(out_csv))}

    rows = []
    for pdb in COORD + NONCOORD:
        if pdb in done:
            print(f"{pdb}: already done")
            continue
        src = os.path.join(A05, pdb)
        d = os.path.join(WORK, pdb)
        os.makedirs(d, exist_ok=True)
        rec = {"pdb_id": pdb, "arm": "coordinating" if pdb in COORD else "non_coordinating"}
        print(f"\n=== {pdb} ({rec['arm']}) ===", flush=True)
        try:
            # receptor+ligand complex, as analysis 05 built it for grid generation
            comp = os.path.join(src, f"{pdb}_reclig.mae")
            if not os.path.exists(comp):
                comp = os.path.join(src, f"{pdb}_prepped.mae")
            if not os.path.exists(comp):
                raise RuntimeError("no prepped receptor found")
            st = next(iter(structure.StructureReader(comp)))

            fe = [a for a in st.atom if a.element == "Fe"]
            if len(fe) != 1:
                raise RuntimeError(f"found {len(fe)} Fe atoms, need exactly 1")
            fe_num = fe[0].index          # 1-based atom number, what -metal_cons wants
            print(f"  Fe atom number {fe_num}", flush=True)

            local = os.path.join(d, os.path.basename(comp))
            if not os.path.exists(local):
                st.write(local)

            gjob = f"{pdb}_mcons"
            found = glob.glob(os.path.join(d, f"{gjob}*.zip"))
            if not found:
                subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT",
                                "-rec_file", os.path.basename(local),
                                "-lig_asl", "res. LIG",
                                "-metal_cons", f"fe#{fe_num}",
                                "-j", gjob],
                               cwd=d, check=False, capture_output=True,
                               text=True, timeout=3600)
                found = glob.glob(os.path.join(d, f"{gjob}*.zip"))
            if not found:
                raise RuntimeError("no constrained grid produced")

            ligf = os.path.join(src, f"{pdb}_ligprep.maegz")
            if not os.path.exists(ligf):
                raise RuntimeError("no ligprep input (analysis 05)")
            djob = f"{pdb}_mdock"
            with open(os.path.join(d, f"{djob}.in"), "w") as f:
                f.write(DOCK_IN.format(grid=os.path.abspath(found[0]),
                                       lig=os.path.abspath(ligf)))
            subprocess.run([os.path.join(SCHRODINGER, "run"), "xglide.py",
                            "-NOJOBID", "-WAIT", f"{djob}.in"],
                           cwd=d, check=False, capture_output=True,
                           text=True, timeout=7200)
            pvs = glob.glob(os.path.join(d, "*_pv.mae")) + \
                  glob.glob(os.path.join(d, djob, "*_pv.mae"))
            pv = pvs[0] if pvs else os.path.join(d, f"{djob}_pv.mae")
            if not os.path.exists(pv):
                raise RuntimeError("no poses produced (constraint may be unsatisfiable)")
            rec["n_poses"] = len(list(structure.StructureReader(pv))) - 1
            rec["status"] = "ok"
            print(f"  {rec['n_poses']} poses", flush=True)
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"  ERROR: {e}", flush=True)
        rows.append(rec)
        allr = rows
        if os.path.exists(out_csv):
            allr = list(csv.DictReader(open(out_csv))) + rows
        keys = sorted({k for r in allr for k in r})
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader(); w.writerows(allr)
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    sys.exit(main())
