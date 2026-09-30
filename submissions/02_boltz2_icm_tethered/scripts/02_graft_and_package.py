"""
Submission 02, step 2: graft ICM-minimised coordinates onto the original Boltz PDBs,
then re-run every submission check.

Grafting rather than using ICM's PDB directly: ICM adds ~4200 hydrogens and rewrites the
file, and the organizers' validator is sensitive to atom naming, duplicate atom names
within a residue, and duplicate atom serials. The original Boltz file already passed all
of that, so only the x/y/z columns of matching heavy atoms are replaced, matched on
(chain, residue number, atom name). Any atom ICM did not return is left at its Boltz
position and counted in the report.

Checks re-run afterwards are the same ones in submission 01's 03_package.py.

Usage:
    ~/micromamba/envs/cyp_posebusters_env/bin/python scripts/02_graft_and_package.py

Outputs:
    submission/<name>.pdb            grafted structures
    submission_icm_tethered.zip
    results/graft_report.csv
"""

import csv
import os
import sys
import zipfile
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, "submissions", "01_boltz2_default", "submission")
WORK = os.path.join(HERE, "work")
OUT = os.path.join(HERE, "submission")
ZIPPATH = os.path.join(HERE, "submission_icm_tethered.zip")
RESULTS = os.path.join(HERE, "results")


def key(line):
    return (line[21], line[22:27].strip(), line[12:16].strip())


def load_min_coords(path):
    """(chain, resnum, atomname) -> xyz, heavy atoms only."""
    out = {}
    for ln in open(path):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = ln[76:78].strip() or ln[12:16].strip()[0]
        if el == "H":
            continue
        out[key(ln)] = (ln[30:38], ln[38:46], ln[46:54])
    return out


def check(lines, name, smiles, tmpdir):
    errs = []
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if len({l[21] for l in atoms}) > 3:
        errs.append("more than 3 chains")
    lig = [l for l in atoms if l[17:20].strip() == "LIG"]
    if not lig:
        errs.append("no LIG residue")
    elif len({(l[21], l[22:27]) for l in lig}) > 1:
        errs.append("more than one LIG residue")
    if len([l for l in atoms if l[17:20].strip() == "HEM"]) != 43:
        errs.append("HEM atom count != 43")
    seen = {}
    for l in atoms:
        k = (l[21], l[22:27])
        if k in seen and seen[k] != l[17:20].strip():
            errs.append(f"duplicate residue {k[1].strip()} chain {k[0]}")
        seen[k] = l[17:20].strip()
    per = defaultdict(list)
    for l in atoms:
        per[(l[21], l[22:27])].append(l[12:16].strip())
    for k, names in per.items():
        if len(names) != len(set(names)):
            errs.append(f"duplicate atom names in {k[0]}{k[1].strip()}")
    ser = [l[6:11] for l in atoms]
    if len(ser) != len(set(ser)):
        errs.append("duplicate atom serials")
    p = os.path.join(tmpdir, f"{name}_lig.pdb")
    with open(p, "w") as f:
        f.writelines(lig + ["END\n"])
    ref = Chem.MolFromSmiles(smiles)
    probe = Chem.MolFromPDBFile(p, removeHs=True, sanitize=False)
    if ref is None or probe is None:
        errs.append("ligand parse failure")
    else:
        try:
            AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except ValueError:
            errs.append("ligand connectivity != SMILES")
    return errs


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    tmpdir = os.path.join(RESULTS, "_tmp")
    os.makedirs(tmpdir, exist_ok=True)

    rows, n_ok = [], 0
    for name in sorted(smiles):
        src = os.path.join(SRC, f"{name}.pdb")
        mn = os.path.join(WORK, f"{name}_min.pdb")
        rec = {"name": name}
        if not (os.path.exists(src) and os.path.exists(mn)):
            rec.update(status="MISSING")
            rows.append(rec)
            continue
        new = load_min_coords(mn)
        out_lines, grafted, missed, disp = [], 0, 0, []
        for ln in open(src):
            if ln.startswith(("ATOM", "HETATM")):
                k = key(ln)
                if k in new:
                    old = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
                    x, y, z = new[k]
                    disp.append(np.linalg.norm(old - np.array([float(x), float(y), float(z)])))
                    ln = ln[:30] + x + y + z + ln[54:]
                    grafted += 1
                else:
                    missed += 1
            out_lines.append(ln)
        with open(os.path.join(OUT, f"{name}.pdb"), "w") as f:
            f.writelines(out_lines)
        errs = check(out_lines, name, smiles[name], tmpdir)
        rec.update(grafted=grafted, not_matched=missed,
                   median_disp=float(np.median(disp)) if disp else None,
                   max_disp=float(np.max(disp)) if disp else None,
                   status="ok" if not errs else "FAIL", errors="; ".join(errs))
        n_ok += not errs
        rows.append(rec)
        print(f"  {name}: grafted {grafted}, unmatched {missed}, "
              f"med {rec['median_disp']:.2f} A, max {rec['max_disp']:.2f} A"
              f"{'' if not errs else ' | FAIL: ' + '; '.join(errs)}")

    pd.DataFrame(rows).to_csv(os.path.join(RESULTS, "graft_report.csv"), index=False)
    print(f"\n{n_ok}/{len(smiles)} pass structural checks")
    if n_ok == len(smiles):
        with zipfile.ZipFile(ZIPPATH, "w", zipfile.ZIP_DEFLATED) as z:
            for name in sorted(smiles):
                z.write(os.path.join(OUT, f"{name}.pdb"), f"{name}.pdb")
        print(f"wrote {ZIPPATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
