"""Submission 07, step 3: convert all 100 AF3 samples to PDB and PrepWizard each.

Mirrors the path that produced submission 05: AF3 mmCIF -> PDB with the ligand renamed
to LIG (AF3 writes LIG_C, which Schrodinger's PDB writer renders as the 3-character UNL),
then PrepWizard with analysis 01's option set. PrepWizard is what took the AF3 samples
from 57/100 to 95/100 PoseBusters passes last time.

Serial and resumable: a sample whose *_pw.pdb already exists is skipped. Compute node only.
"""
import csv, glob, os, shutil, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
WORK = os.path.join(HERE, "work")
SCHRODINGER = os.environ["SCHRODINGER"]
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")


def main():
    if not (os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1"):
        raise SystemExit("REFUSING to run Schrodinger outside SLURM.")
    from schrodinger import structure
    os.makedirs(WORK, exist_ok=True)
    names = sorted(r["Molecule_Name"].strip() for r in csv.DictReader(
        open(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
        if r["Molecule_Name"].strip())

    cifs = []
    for n in names:
        found = sorted(glob.glob(os.path.join(RES, n.lower(), "seed-*_sample-*",
                                              "*_model.cif")))
        if len(found) != 5:
            print(f"  WARNING: {n} has {len(found)} samples, expected 5")
        cifs += [(n, os.path.basename(os.path.dirname(f)), f) for f in found]
    print(f"{len(cifs)} samples found (expected 100)")
    if len(cifs) < 100:
        raise SystemExit(f"FATAL: only {len(cifs)} samples - AF3 run incomplete")

    for i, (n, tag, cif) in enumerate(cifs, 1):
        base = f"{n}__{tag}"
        d = os.path.join(WORK, base)
        os.makedirs(d, exist_ok=True)
        prepped = os.path.join(d, f"{base}_pw.pdb")
        if os.path.exists(prepped):
            continue
        raw = os.path.join(d, f"{base}.pdb")
        if not os.path.exists(raw):
            st = next(iter(structure.StructureReader(cif)))
            k = 0
            for a in st.atom:
                rn = a.pdbres.strip().upper()
                if rn == "UNL" or rn.startswith("LIG"):
                    a.pdbres = "LIG "
                    k += 1
            if k == 0:
                raise RuntimeError(f"{base}: no ligand residue to rename")
            st.write(raw)
        print(f"[{i}/{len(cifs)}] {base}", flush=True)
        subprocess.run([PREPWIZARD, "-NOJOBID", "-rehtreat", "-captermini",
                        "-propka_pH", "7.4", "-epik_pH", "7.4",
                        "-f", "OPLS_2005", "-rmsd", "0.3",
                        os.path.basename(raw), os.path.basename(prepped)],
                       cwd=d, check=False, capture_output=True, text=True, timeout=3600)
        if not os.path.exists(prepped):
            print(f"    WARNING: prepwizard produced no output for {base}", flush=True)

    done = len(glob.glob(os.path.join(WORK, "*", "*_pw.pdb")))
    print(f"\nprepared {done}/{len(cifs)}")
    return 0 if done >= len(cifs) * 0.9 else 1


if __name__ == "__main__":
    sys.exit(main())
