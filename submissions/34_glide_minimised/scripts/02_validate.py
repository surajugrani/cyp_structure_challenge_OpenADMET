"""Validate and package submission 34.

The ligand coordinates changed, so this is a full re-validation: chain layout, HEM
integrity, ligand connectivity against the input SMILES, then PoseBusters in `dock`
configuration. Not a copy-only submission, so no checksum shortcut.

struct_errors and pb_fail are imported from submission 23's script rather than copied,
so the two cannot silently diverge.
"""
import csv, importlib.util, os, sys, zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "submissions", "23_esmfold2_default", "scripts",
                   "02_convert_validate.py")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_glide_minimised.zip")

spec = importlib.util.spec_from_file_location("conv23", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    tmp = os.path.join(RES, "_tmp"); os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = m.PoseBusters(config="dock")
    rows, ns, npb = [], 0, 0
    for n in sorted(smiles):
        dst = os.path.join(OUT, f"{n}.pdb")
        if not os.path.exists(dst):
            print(f"  {n}: MISSING"); continue
        se = m.struct_errors(open(dst).readlines(), n, smiles[n], tmp)
        pf = m.pb_fail(dst, smiles[n], pb)
        ns += not se; npb += not pf
        rows.append({"ligand": n, "struct_errors": "; ".join(se),
                     "pb_failed": "; ".join(pf)})
        print(f"  {n}: " + ("ok" if not (se or pf) else f"struct {se} | PB {pf}"),
              flush=True)

    with open(os.path.join(RES, "validation.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, ["ligand", "struct_errors", "pb_failed"])
        w.writeheader(); w.writerows(rows)
    print(f"\nFINAL structural {ns}/{len(rows)}   PoseBusters {npb}/{len(rows)}")

    if ns == npb == len(rows) == 20:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for r in rows:
                z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")
        print(f"wrote {ZIP}")
        return 0
    print("NOT packaged - below 20/20")
    return 1


if __name__ == "__main__":
    sys.exit(main())
