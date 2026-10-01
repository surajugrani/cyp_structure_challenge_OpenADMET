"""Glide rescoring of all 100 AF3 samples, to test it as a sample selector.

Idea: AF3's own confidences do not predict pose validity (analysis 10 - has_clash is
identically 0 across all 100 samples). Maybe a physics score does. For each of the 20
ligands x 5 samples: PrepWizard the complex, build a Glide grid on that sample's own
receptor conformation, and score the pose in place. Then pick the best-scoring sample per
ligand, and check whether Glide's score tracks any AF3 metric.

Each sample is a DIFFERENT receptor conformation, so each needs its own grid - hence 100
preps and 100 grids rather than 20.

`mininplace`, not `inplace`: submission 02 established that `inplace` returns emodel and
einternal identically 0.00 because Glide only computes them in the post-docking stage,
and analysis 09 found emodel the most informative Glide term for pose ranking. PrepWizard
has already relaxed the pose, so mininplace moves it only slightly further (0.15 A median,
measured in submission 03) while returning the full term set.

Expect a modest effect at best: analysis 09 measured Glide's within-ligand pose-ranking
power directly and found rho = +0.086 (gscore), +0.057 (emodel), +0.115 (lipo).

Serial and resumable - a sample whose row is already in the results CSV is skipped.
Compute node only.
"""
import csv
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SAMPLES = os.path.join(ROOT, "analyses", "10_af3_selection_vs_posebusters",
                       "work", "samples")
WORK = os.path.join(HERE, "work")
RESULTS = os.path.join(HERE, "results")
OUT_CSV = os.path.join(RESULTS, "glide_scores.csv")

SCHRODINGER = os.environ["SCHRODINGER"]
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")
GEN_GRIDS = os.path.join(SCHRODINGER, "utilities", "generate_glide_grids")
GLIDE = os.path.join(SCHRODINGER, "glide")

SCORE_IN = """GRIDFILE   {grid}
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
    raise SystemExit("REFUSING to run Schrodinger outside SLURM. Use sbatch.")


def done_keys():
    if not os.path.exists(OUT_CSV):
        return set()
    return {(r["ligand"], r["sample"]) for r in csv.DictReader(open(OUT_CSV))}


def append(row):
    new = not os.path.exists(OUT_CSV)
    with open(OUT_CSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "sample", "status",
                                          "fe_lig_bonds_broken"] + TERMS,
                           extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def main():
    require_compute_node()
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    from schrodinger import structure

    skip = done_keys()
    files = sorted(glob.glob(os.path.join(SAMPLES, "*.pdb")))
    print(f"{len(files)} samples, {len(skip)} already scored")

    for i, src in enumerate(files, 1):
        base = os.path.basename(src)[:-4]
        lig_name, tag = base.split("__")
        if (lig_name, tag) in skip:
            continue
        d = os.path.join(WORK, base)
        os.makedirs(d, exist_ok=True)
        row = {"ligand": lig_name, "sample": tag}
        print(f"[{i}/{len(files)}] {base}", flush=True)
        try:
            # 1. PrepWizard the whole complex, same options as submission 03.
            #    prepwizard resolves paths relative to cwd, so the input is copied in.
            import shutil
            local = os.path.join(d, os.path.basename(src))
            if not os.path.exists(local):
                shutil.copy(src, local)
            prepped = os.path.join(d, f"{base}_pw.pdb")
            if not os.path.exists(prepped):
                subprocess.run([PREPWIZARD, "-NOJOBID", "-rehtreat", "-captermini",
                                "-propka_pH", "7.4", "-epik_pH", "7.4",
                                "-f", "OPLS_2005", "-rmsd", "0.3",
                                os.path.basename(local), os.path.basename(prepped)],
                               cwd=d, check=False, capture_output=True,
                               text=True, timeout=3600)
            if not os.path.exists(prepped):
                raise RuntimeError("prepwizard produced no output")

            # 2. split receptor (protein + HEM) from ligand
            st = next(iter(structure.StructureReader(prepped)))
            lig_idx = [a.index for a in st.atom if a.pdbres.strip() == "LIG"]
            if not lig_idx:
                raise RuntimeError("no LIG residue after prepwizard")

            # Break any covalent Fe-ligand bond before splitting. 49 of the 100 AF3
            # samples place a ligand atom within 2.2 A of the haem iron, and PrepWizard
            # (-rehtreat treats metals) bonds them. Extracting the ligand across such a
            # bond would hand Glide a fragment with a dangling valence, and Glide cannot
            # score a covalently attached ligand in place anyway. Breaking it first turns
            # the interaction into a normal non-covalent one, which is what the Glide
            # metal term is meant to describe.
            lig_set = set(lig_idx)
            n_fe_bonds = 0
            for a in st.atom:
                if a.element != "Fe":
                    continue
                for b in list(a.bonded_atoms):
                    if b.index in lig_set:
                        st.deleteBond(a, b)
                        n_fe_bonds += 1
            row["fe_lig_bonds_broken"] = n_fe_bonds

            ligst = st.extract(lig_idx)
            ligf = os.path.join(d, f"{base}_lig.mae")
            ligst.write(ligf)

            # 3. Grid on this sample's own receptor conformation.
            #    Uses analysis 01's proven invocation: command-line -rec_file/-lig_asl/-j,
            #    NOT a keyword .in file - generate_glide_grids in schrodinger/2025-3
            #    treats an .in file as a CSV and fails with "Input file header is missing
            #    'rec_file' column". The whole prepped complex is passed as the receptor
            #    and -lig_asl removes the ligand atoms from the grid while they still
            #    define the box centre. Output is <job>-gridgen.zip, not <job>.zip.
            complex_mae = os.path.join(d, f"{base}_complex.mae")
            if not os.path.exists(complex_mae):
                st.write(complex_mae)
            gjob = f"{base}_grid"
            found = glob.glob(os.path.join(d, f"{gjob}*.zip"))
            if not found:
                subprocess.run([GEN_GRIDS, "-NOJOBID", "-WAIT",
                                "-rec_file", os.path.basename(complex_mae),
                                "-lig_asl", "res. LIG",
                                "-j", gjob],
                               cwd=d, check=False, capture_output=True,
                               text=True, timeout=3600)
                found = glob.glob(os.path.join(d, f"{gjob}*.zip"))
            if not found:
                raise RuntimeError("no grid produced")

            # 4. score the pose in place
            sjob = f"{base}_score"
            with open(os.path.join(d, f"{sjob}.in"), "w") as f:
                f.write(SCORE_IN.format(grid=os.path.basename(found[0]),
                                        lig=os.path.basename(ligf)))
            subprocess.run([GLIDE, "-NOJOBID", "-WAIT", f"{sjob}.in"],
                           cwd=d, check=False, capture_output=True,
                           text=True, timeout=3600)
            csvp = os.path.join(d, f"{sjob}.csv")
            if not os.path.exists(csvp):
                raise RuntimeError("no glide results csv")
            import pandas as pd
            g = pd.read_csv(csvp)
            if g.empty:
                raise RuntimeError("glide returned no pose")
            r0 = g.iloc[0]
            for t in TERMS:
                col = f"r_i_glide_{t}"
                if col in g.columns:
                    try:
                        row[t] = float(r0[col])
                    except (TypeError, ValueError):
                        pass
            row["status"] = "ok"
            print(f"    gscore={row.get('gscore')} emodel={row.get('emodel')}", flush=True)
        except Exception as e:
            row["status"] = f"error: {e}"
            print(f"    ERROR: {e}", flush=True)
        append(row)

    print(f"\nwrote {OUT_CSV}")


if __name__ == "__main__":
    sys.exit(main())
