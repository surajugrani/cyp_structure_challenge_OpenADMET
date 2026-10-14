"""Build one self-contained PyMOL session (.pse) per submission.

    pymol -cq tools/make_pymol_session.py

Processes every submission that has a submission/ directory. Set CYP_SUBMISSION to a
single directory to limit it:

    CYP_SUBMISSION=submissions/12_icm_redock_af3 pymol -cq tools/make_pymol_session.py

Writes <submission>/visualize/<submission>.pse holding all 20 complexes as named objects
(OCNT-2312554, ...), already styled. Double-click it, or `pymol <file>.pse` - nothing to
run, no dependency on the 20 individual PDBs once it is written.

A .pse is used rather than a combined PDB because PDB has no concept of named objects,
and MODEL/ENDMDL states require every state to share an identical atom composition -
which these do not, since each complex carries a different ligand. Loading them as states
makes PyMOL misassign coordinates and corrupt both ligand and backbone.

NOTE: PyMOL does not set __name__ == "__main__" when it executes a script, so this file
calls its entry point unconditionally at the bottom. Guarding it meant nothing ran and the
job produced no output at all. Argument passing through `pymol -- args` is likewise
unreliable, hence the environment variable.

Objects are superposed on chain A so the complexes overlay, with the ligands
colour-coded by element and the haem iron shown as a sphere.
"""
import glob
import os

from pymol import cmd

# PyMOL execs this file with __file__ pointing into its own site-packages, so the repo
# root cannot be derived from it. Use CYP_ROOT if set, else the working directory (the
# sbatch wrapper cds to SLURM_SUBMIT_DIR, i.e. the repo root).
ROOT = os.environ.get("CYP_ROOT") or os.getcwd()
if not os.path.isdir(os.path.join(ROOT, "submissions")):
    raise SystemExit(f"no submissions/ under {ROOT} - run from the repo root "
                     f"or set CYP_ROOT")


def build(subdir):
    src = os.path.join(subdir, "submission")
    pdbs = sorted(glob.glob(os.path.join(src, "*.pdb")))
    if not pdbs:
        print(f"{os.path.basename(subdir):28} no submission/*.pdb, skipped")
        return
    out = os.path.join(subdir, "visualize")
    os.makedirs(out, exist_ok=True)
    dst = os.path.join(out, f"{os.path.basename(subdir)}.pse")

    cmd.reinitialize()
    names = []
    for p in pdbs:
        n = os.path.basename(p)[:-4]
        cmd.load(p, n)
        names.append(n)

    for n in names[1:]:
        try:
            cmd.align(f"polymer and {n}", f"polymer and {names[0]}")
        except Exception as e:
            print(f"    align failed for {n}: {e}")

    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    cmd.show("sticks", "resn LIG")
    cmd.show("sticks", "resn HEM")
    cmd.show("spheres", "elem Fe")
    cmd.set("sphere_scale", 0.3, "elem Fe")
    cmd.color("grey80", "polymer")
    cmd.color("firebrick", "resn HEM")
    cmd.util.cbay("resn LIG")
    cmd.set("cartoon_transparency", 0.6)
    cmd.dss()
    cmd.orient("resn LIG")
    cmd.save(dst)
    print(f"{os.path.basename(subdir):28} {len(names):2d} objects -> {dst}")


def main():
    one = os.environ.get("CYP_SUBMISSION")
    if one:
        targets = [one if os.path.isabs(one) else os.path.join(ROOT, one)]
    else:
        base = os.path.join(ROOT, "submissions")
        targets = [os.path.join(base, d) for d in sorted(os.listdir(base))
                   if os.path.isdir(os.path.join(base, d, "submission"))]
    print(f"building sessions for {len(targets)} submission(s)", flush=True)
    for t in targets:
        build(t)
    print("done", flush=True)


# PyMOL runs this file with __name__ != "__main__" - call unconditionally.
main()
