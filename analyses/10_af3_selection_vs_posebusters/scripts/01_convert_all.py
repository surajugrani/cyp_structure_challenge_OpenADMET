"""Analysis 10, step 1: convert all 5 AF3 samples of all 20 ligands to submission-style PDB.

100 structures. Same conversion the submission uses: read the mmCIF, rename the ligand
residue to LIG (AF3 calls it LIG_C, which Schrodinger's PDB writer renders as the
3-character UNL), keep HEM as is.

Run under $SCHRODINGER/run python3.
"""
import csv, glob, os
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
AF3 = os.path.join(ROOT, "submissions", "04_af3_default", "results")
OUT = os.path.join(HERE, "work", "samples")


def main():
    os.makedirs(OUT, exist_ok=True)
    names = sorted(r["Molecule_Name"].strip() for r in csv.DictReader(
        open(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
        if r["Molecule_Name"].strip())
    n_done = 0
    for n in names:
        for cif in sorted(glob.glob(os.path.join(AF3, n.lower(),
                                                 "seed-*_sample-*", "*_model.cif"))):
            tag = os.path.basename(os.path.dirname(cif))          # seed-1_sample-K
            st = next(iter(structure.StructureReader(cif)))
            k = 0
            for a in st.atom:
                rn = a.pdbres.strip().upper()
                if rn == "UNL" or rn.startswith("LIG"):
                    a.pdbres = "LIG "
                    k += 1
            if k == 0:
                raise RuntimeError(f"{n} {tag}: no ligand residue to rename")
            st.write(os.path.join(OUT, f"{n}__{tag}.pdb"))
            n_done += 1
    print(f"converted {n_done} structures -> {OUT}")
    assert n_done == 100, f"expected 100, got {n_done}"


if __name__ == "__main__":
    main()
