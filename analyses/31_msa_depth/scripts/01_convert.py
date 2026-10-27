"""Convert a variant ESMFold2 run's mmCIF output to PDB, nothing else.

Submission 23's 02_convert_validate.py already has the only mmCIF -> PDB converter in
the project, and it is the one whose output the organizers have accepted five times, so
it is imported rather than reimplemented. Its PoseBusters and zipping steps are skipped:
this is an analysis, not a submission, and all we need is a directory of PDBs that
tools/ligand_rmsd.py can read.

    python3 01_convert.py RAW_DIR OUT_DIR
"""
import os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT, "submissions", "23_esmfold2_default", "scripts"))


def main():
    raw, out = sys.argv[1], sys.argv[2]
    from importlib import import_module
    cif_to_pdb = import_module("02_convert_validate").cif_to_pdb
    os.makedirs(out, exist_ok=True)
    n = 0
    for f in sorted(os.listdir(raw)):
        if not f.endswith(".cif"):
            continue
        natoms, nhem = cif_to_pdb(os.path.join(raw, f),
                                  os.path.join(out, f[:-4] + ".pdb"))
        if nhem != 43:
            print(f"  WARNING {f}: HEM has {nhem} heavy atoms, expected 43")
        n += 1
    print(f"converted {n} structures -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
