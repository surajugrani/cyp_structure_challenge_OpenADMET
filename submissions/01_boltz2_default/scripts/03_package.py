"""
Submission 01, step 3: turn Boltz-2 output into a validator-clean submission zip.

Boltz names EVERY non-polymer chain `LIG`, so its output has two LIG residues - the haem
(chain B, 43 atoms) and the blinded ligand (chain C). The organizers'
validation/structure_validation.py requires exactly one residue named LIG and rejects
anything else with "Found 2 'LIG' residues, expected 1". So chain B is renamed to HEM and
chain C is left as LIG.

Everything else the validator checks is verified here rather than discovered on upload:

  - zip contains .pdb files, NOT .cif. The .cif files in the tutorial repo's
    inputs/boltz-2_cyp3a4_predictions/ are example inputs; the validator only collects
    names ending .pdb and otherwise fails with "Zip file contains no PDB files".
  - 20 files, stems equal to the Molecule_Name values.
  - exactly one LIG residue per file.
  - at most 3 chains (protein + haem + ligand).
  - no duplicate (chain, resid) pairs, no duplicate atom names within a residue, and
    globally unique atom serials - these break OpenStructure's mapping.
  - ligand connectivity matches the given SMILES, via RDKit
    AssignBondOrdersFromTemplate, which is the same check the organizers run.

Usage:
    $SCHRODINGER/run python3 scripts/03_package.py

Outputs:
    submission/<Molecule_Name>.pdb      cleaned structures
    submission_boltz2_default.zip       the archive to upload
    results/packaging_report.csv        per-structure checks
"""

import csv
import os
import shutil
import sys
import zipfile
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
PRED = os.path.join(HERE, "results", "boltz_results_inputs", "predictions")
OUTDIR = os.path.join(HERE, "submission")
ZIPPATH = os.path.join(HERE, "submission_boltz2_default.zip")

HEME_CHAIN = "B"
LIG_CHAIN = "C"
N_HEME_ATOMS = 43


def read_smiles():
    path = os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")
    with open(path) as f:
        return {r["Molecule_Name"].strip(): r["SMILES"].strip()
                for r in csv.DictReader(f) if r["Molecule_Name"].strip()}


def rename_heme(lines):
    """Boltz calls every het chain LIG; the haem chain must become HEM."""
    out, n = [], 0
    for ln in lines:
        if ln.startswith("HETATM") and ln[21] == HEME_CHAIN and ln[17:20] == "LIG":
            ln = ln[:17] + "HEM" + ln[20:]
            n += 1
        out.append(ln)
    return out, n


def check(lines, name, smiles, tmpdir):
    """Replicate the organizers' checks locally."""
    errs = []
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]

    chains = {l[21] for l in atoms}
    if len(chains) > 3:
        errs.append(f"{len(chains)} chains (max 3)")

    lig = [l for l in atoms if l[17:20].strip() == "LIG"]
    lig_res = {(l[21], l[22:27]) for l in lig}
    if not lig:
        errs.append("no LIG residue")
    elif len(lig_res) > 1:
        errs.append(f"{len(lig_res)} LIG residues (expected 1)")

    hem = [l for l in atoms if l[17:20].strip() == "HEM"]
    if len(hem) != N_HEME_ATOMS:
        errs.append(f"HEM has {len(hem)} atoms (expected {N_HEME_ATOMS})")

    # duplicate (chain, resid) labelling different residues
    seen = {}
    for l in atoms:
        key = (l[21], l[22:27])
        rn = l[17:20].strip()
        if key in seen and seen[key] != rn:
            errs.append(f"duplicate residue {key[1].strip()} in chain {key[0]}")
        seen[key] = rn

    # duplicate atom names within one residue
    per_res = defaultdict(list)
    for l in atoms:
        per_res[(l[21], l[22:27])].append(l[12:16].strip())
    for key, names in per_res.items():
        dup = [n for n, c in Counter(names).items() if c > 1]
        if dup:
            errs.append(f"duplicate atom names in {key[0]}{key[1].strip()}: {dup[:3]}")

    # globally unique atom serials
    serials = [l[6:11] for l in atoms]
    if len(serials) != len(set(serials)):
        errs.append(f"{len(serials) - len(set(serials))} duplicate atom serials")

    # ligand connectivity vs the given SMILES
    ligpdb = os.path.join(tmpdir, f"{name}_lig.pdb")
    with open(ligpdb, "w") as f:
        f.writelines(lig + ["END\n"])
    ref = Chem.MolFromSmiles(smiles)
    probe = Chem.MolFromPDBFile(ligpdb, removeHs=True, sanitize=False)
    if ref is None:
        errs.append("cannot parse expected SMILES")
    elif probe is None:
        errs.append("RDKit cannot parse LIG residue")
    else:
        try:
            AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except ValueError:
            errs.append("ligand connectivity != expected SMILES")
    return errs, len(lig)


def main():
    smiles = read_smiles()
    os.makedirs(OUTDIR, exist_ok=True)
    tmpdir = os.path.join(HERE, "results", "_tmp_lig")
    os.makedirs(tmpdir, exist_ok=True)

    rows, ok_n = [], 0
    for name in sorted(smiles):
        src = os.path.join(PRED, name, f"{name}_model_0.pdb")
        rec = {"name": name}
        if not os.path.exists(src):
            rec.update(status="MISSING", errors="no prediction")
            rows.append(rec)
            print(f"  {name}: MISSING")
            continue
        lines = open(src).readlines()
        lines, n_renamed = rename_heme(lines)
        dst = os.path.join(OUTDIR, f"{name}.pdb")
        with open(dst, "w") as f:
            f.writelines(lines)
        errs, n_lig = check(lines, name, smiles[name], tmpdir)
        rec.update(heme_atoms_renamed=n_renamed, lig_atoms=n_lig,
                   status="ok" if not errs else "FAIL",
                   errors="; ".join(errs))
        ok_n += not errs
        rows.append(rec)
        print(f"  {name}: {'ok' if not errs else 'FAIL - ' + '; '.join(errs)}")

    shutil.rmtree(tmpdir, ignore_errors=True)
    import pandas as pd
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "packaging_report.csv"),
                              index=False)

    print(f"\n{ok_n}/{len(smiles)} pass all checks")
    if ok_n != len(smiles):
        print("NOT writing the zip - fix the failures above first.")
        return 1
    with zipfile.ZipFile(ZIPPATH, "w", zipfile.ZIP_DEFLATED) as z:
        for name in sorted(smiles):
            z.write(os.path.join(OUTDIR, f"{name}.pdb"), f"{name}.pdb")
    print(f"wrote {ZIPPATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
