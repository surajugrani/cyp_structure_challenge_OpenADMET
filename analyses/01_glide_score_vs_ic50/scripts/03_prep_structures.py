"""
Step 3: Download, clean and prepare each curated structure with PrepWizard.

For each ligand in data/curated_set.csv:
  1. Download <PDB>.pdb from RCSB into structures/<PDB>/
  2. Keep one chain: protein + HEM + all copies of the ligand of interest in that
     chain. Waters, buffer/cryo molecules and other chains are removed.
     - Scoring copy = the copy closest to the heme Fe.
     - Other copies in the same chain (2V0M, 4K9T) are KEPT as part of the receptor
       (decision "option A", see README).
  3. PrepWizard: add H, fill missing side chains, cap termini, PROPKA/Epik at pH 7.4,
     H-bond optimization, restrained minimization (OPLS_2005, heavy-atom RMSD 0.3 A).
     If side-chain filling fails (e.g. Prime licence), retries without it.
  4. Sanity checks on the prepped complex (ligand present, Fe bonds, heme charge).

Usage (from the analysis folder, normally via scripts/run_prep_and_score.sh):
    $SCHRODINGER/run python3 scripts/03_prep_structures.py [--skip-prep] [--only 2V0M,4K9T]

Outputs (git-ignored):
    structures/<PDB>/<PDB>.pdb, <PDB>_clean.mae, <PDB>_prepped.mae, prep_info.json
Outputs (committed):
    results/prep_summary.csv
"""

import argparse
import json
import os
import subprocess
import time

import numpy as np
import pandas as pd
import requests
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(HERE, "data")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS_DIR = os.path.join(HERE, "results")

SCHRODINGER = os.environ["SCHRODINGER"]
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard")
# NOTE: in schrodinger/2025-3, -f takes the force field name ("OPLS_2005"), not "2005"
# Residues that are never the ligand of interest, used by the name-match fallback.
SOLVENT_AND_COFACTOR = {
    "HEM", "HOH", "WAT", "GOL", "EDO", "DMS", "PEG", "PGE", "P6G", "1PE", "2PE",
    "PG4", "TEG", "SO4", "PO4", "ACT", "MPD", "TRS", "EPE", "FMT", "ACY", "NO3",
    "CL", "NA", "MG", "K", "ZN", "CA", "IOD", "BR", "UNX", "UNK",
    # PEG/glycol fragments - PG0 made the A1ASU fallback ambiguous against UNL
    "PG0", "PG5", "PG6", "7PE", "12P", "15P", "MRD", "BU3", "IPA", "EOH", "MOH",
    "CIT", "MES", "BME", "DTT", "SCN", "AZI", "URE", "NH4", "BCT", "OXL", "TAR", "MLI",
}

PREP_OPTS = ["-NOJOBID", "-rehtreat", "-captermini", "-propka_pH", "7.4",
             "-epik_pH", "7.4", "-f", "OPLS_2005", "-rmsd", "0.3"]


def download_pdb(pdb_id, out_path, retries=3):
    """Download a structure, returning the path actually written.

    Falls back from .pdb to .cif: RCSB no longer distributes legacy PDB-format files for
    newer or large depositions, so entries like 9PLK/9COW/9BBB/9COR 404 on .pdb while the
    mmCIF exists. Also retries on 5xx, which RCSB returns transiently (seen for 5A1P,
    6BD5). Schrodinger's StructureReader reads both formats, keyed on the file extension,
    so callers must use the RETURNED path rather than assuming .pdb.
    """
    base = os.path.splitext(out_path)[0]
    for ext in (".pdb", ".cif"):
        candidate = base + ext
        if os.path.exists(candidate) and os.path.getsize(candidate) > 0:
            return candidate

    last = None
    for ext in (".pdb", ".cif"):
        candidate = base + ext
        for attempt in range(retries):
            try:
                r = requests.get(f"https://files.rcsb.org/download/{pdb_id}{ext}",
                                 timeout=120, verify=True)
                if r.status_code == 404:
                    last = f"404 for {pdb_id}{ext}"
                    break                      # no point retrying a genuine 404
                r.raise_for_status()
                with open(candidate, "w") as f:
                    f.write(r.text)
                return candidate
            except requests.RequestException as e:
                last = str(e)
                if attempt < retries - 1:
                    time.sleep(2 ** attempt)   # transient 5xx / network blip
    raise RuntimeError(f"could not download {pdb_id} as .pdb or .cif: {last}")


def clean_structure(pdb_path, chain, comp_id):
    """Return (clean structure, info dict)."""
    st = structure.StructureReader.read(pdb_path)

    fe = [a for a in st.atom if a.element == "Fe" and a.pdbres.strip() == "HEM" and a.chain == chain]
    if len(fe) != 1:
        raise RuntimeError(f"expected 1 heme Fe in chain {chain}, found {len(fe)}")
    fe_xyz = np.array(fe[0].xyz)

    # Ligand copies in this chain, keyed by residue number.
    copies = {}
    for a in st.atom:
        if a.chain == chain and a.pdbres.strip() == comp_id:
            copies.setdefault(a.resnum, []).append(a)

    matched_name = comp_id
    if not copies:
        # Fall back for 5-character CCD codes (e.g. A1CIW in recent depositions):
        # Schrodinger's pdbres field holds only 4 characters, so the reader reports these
        # as "UNL" and an exact name match can never succeed. Identify the ligand instead
        # as the single non-standard, non-cofactor, non-solvent residue in the chain.
        candidates = {}
        for a in st.atom:
            if a.chain != chain or a.getResidue().isStandardResidue():
                continue
            res = a.pdbres.strip()
            if res in SOLVENT_AND_COFACTOR:
                continue
            candidates.setdefault((res, a.resnum), []).append(a)
        if not candidates:
            raise RuntimeError(f"ligand {comp_id} not found in chain {chain}")
        names = {res for res, _ in candidates}
        if len(names) != 1:
            raise RuntimeError(
                f"ligand {comp_id} not found in chain {chain}; ambiguous fallback "
                f"candidates {sorted(names)}")
        matched_name = names.pop()
        for (_, resnum), atoms in candidates.items():
            copies[resnum] = atoms
        print(f"    note: {comp_id} not matchable by name (5-char CCD?); "
              f"using non-standard residue '{matched_name}' instead")
    dist_to_fe = {rn: min(np.linalg.norm(np.array(a.xyz) - fe_xyz) for a in atoms)
                  for rn, atoms in copies.items()}
    score_resnum = min(dist_to_fe, key=dist_to_fe.get)

    keep = []
    for a in st.atom:
        if a.chain != chain:
            continue
        res = a.pdbres.strip()
        is_protein = a.getResidue().isStandardResidue()
        if is_protein or res == "HEM" or res == matched_name:
            keep.append(a.index)
    clean = st.extract(keep, copy_props=True)
    clean.title = os.path.basename(pdb_path).split(".")[0]

    info = {
        "chain": chain,
        "comp_id": matched_name,        # what the file actually calls it
        "requested_comp_id": comp_id,   # what the RCSB report calls it
        "score_resnum": int(score_resnum),
        "score_min_dist_to_fe": round(float(dist_to_fe[score_resnum]), 2),
        "other_copies_kept_in_receptor": [int(r) for r in copies if r != score_resnum],
        "n_atoms_clean": clean.atom_total,
    }
    return clean, info


def require_compute_node(what):
    """Refuse to launch CPU-heavy Schrodinger tools outside a SLURM allocation.

    PrepWizard, grid generation and Glide must not run on the shared login node.
    Set ALLOW_LOGIN_NODE=1 to override (don't).
    """
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("ALLOW_LOGIN_NODE") == "1":
        return
    raise SystemExit(
        f"REFUSING to run {what} outside SLURM - this would load the login node.\n"
        f"Submit it instead, e.g.:  sbatch scripts/run_*.sh"
    )


def run_prepwizard(in_mae, out_mae, work_dir):
    """Run PrepWizard with side-chain filling; fall back without it on failure."""
    require_compute_node("PrepWizard")
    for extra, label in ((["-fillsidechains"], "with -fillsidechains"), ([], "without -fillsidechains")):
        cmd = [PREPWIZARD] + PREP_OPTS + extra + [os.path.basename(in_mae), os.path.basename(out_mae)]
        print("    $", " ".join(cmd))
        res = subprocess.run(cmd, cwd=work_dir)
        if res.returncode == 0 and os.path.exists(out_mae):
            return label
        print(f"    PrepWizard failed {label}")
    raise RuntimeError("PrepWizard failed")


def normalize_heme(prepped_mae):
    """Force every heme into the same protonation/charge state.

    PrepWizard is not deterministic here: across the 22 structures it produced two
    different heme states (protonated vs deprotonated propionates; mono- vs dianionic
    porphyrin). Since the heme sits directly under the ligand, an inconsistent cofactor
    charge makes Glide scores non-comparable between structures, so we enforce:

      - propionates deprotonated  -> both carboxylate O at -1  (pH 7.4)
      - porphyrin dianion         -> the two pyrrole N-H deprotonated, N at -1
      - Fe left as PrepWizard assigned it (+2)

    giving a heme formal charge of -2 in every structure.

    Returns {"propionates": n, "pyrrole_N": n} counting what this call had to change.
    """
    st = structure.StructureReader.read(prepped_mae)
    to_delete = []
    fixed = {"propionates": 0, "pyrrole_N": 0}

    for atom in st.atom:
        if atom.pdbres.strip() != "HEM":
            continue
        hydrogens = [b for b in atom.bonded_atoms if b.element == "H"]
        if not hydrogens:
            continue

        # Carboxylic acid O-H -> carboxylate
        if atom.element == "O":
            carbons = [b for b in atom.bonded_atoms if b.element == "C"]
            if len(hydrogens) == 1 and len(carbons) == 1 and \
                    sum(1 for b in carbons[0].bonded_atoms if b.element == "O") == 2:
                to_delete.append(hydrogens[0].index)
                atom.formal_charge = -1
                fixed["propionates"] += 1

        # Pyrrole N-H (only the four N coordinating Fe) -> N(-)
        elif atom.element == "N" and any(b.element == "Fe" for b in atom.bonded_atoms):
            if len(hydrogens) == 1:
                to_delete.append(hydrogens[0].index)
                atom.formal_charge = -1
                fixed["pyrrole_N"] += 1

    if to_delete:
        st.deleteAtoms(to_delete)
        st.write(prepped_mae)
    return fixed


def check_prepped(prepped_mae, info):
    st = structure.StructureReader.read(prepped_mae)
    lig = [a for a in st.atom if a.chain == info["chain"] and a.resnum == info["score_resnum"]
           and a.pdbres.strip() == info["comp_id"]]
    heme = [a for a in st.atom if a.pdbres.strip() == "HEM"]
    fe = [a for a in heme if a.element == "Fe"]
    checks = {
        "lig_atoms_prepped": len(lig),
        "lig_formal_charge": sum(a.formal_charge for a in lig),
        "heme_formal_charge": sum(a.formal_charge for a in heme),
        # Bonds to Fe (PrepWizard uses zero-order bonds): expect 4 porphyrin N + Cys S,
        # +1 if the ligand coordinates the iron.
        "fe_bonded_to": ";".join(sorted(f"{b.pdbres.strip()}:{b.element}" for b in fe[0].bonded_atoms))
        if fe else "NO_FE",
    }
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-prep", action="store_true", help="download + clean only")
    ap.add_argument("--only", help="comma-separated PDB IDs")
    args = ap.parse_args()

    cur = pd.read_csv(os.path.join(DATA_DIR, "curated_set.csv"))
    if args.only:
        cur = cur[cur["pdb_id"].isin(args.only.split(","))]
    os.makedirs(RESULTS_DIR, exist_ok=True)

    rows = []
    for _, row in cur.iterrows():
        pdb, comp, chain = row["pdb_id"], row["comp_id"], row["lig_chain"]
        wdir = os.path.join(STRUCT_DIR, pdb)
        os.makedirs(wdir, exist_ok=True)
        print(f"\n== {pdb} ({comp}, chain {chain})")
        out = {"pdb_id": pdb, "comp_id": comp}
        try:
            pdb_path = download_pdb(pdb, os.path.join(wdir, f"{pdb}.pdb"))
            clean, info = clean_structure(pdb_path, chain, comp)
            clean_mae = os.path.join(wdir, f"{pdb}_clean.mae")
            clean.write(clean_mae)
            print(f"    scoring copy {chain}:{info['score_resnum']} "
                  f"({info['score_min_dist_to_fe']} A from Fe); "
                  f"other copies kept: {info['other_copies_kept_in_receptor']}")
            out.update(info)

            if not args.skip_prep:
                prepped_mae = os.path.join(wdir, f"{pdb}_prepped.mae")
                if not os.path.exists(prepped_mae):
                    info["prep_mode"] = run_prepwizard(clean_mae, prepped_mae, wdir)
                # Idempotent: safe to re-run on an already-normalized structure.
                info.update({f"heme_fixed_{k}": v
                             for k, v in normalize_heme(prepped_mae).items()})
                checks = check_prepped(prepped_mae, info)
                info.update(checks)
                out.update(info)
                print(f"    checks: {checks}")

            with open(os.path.join(wdir, "prep_info.json"), "w") as f:
                json.dump(info, f, indent=1)
            out["status"] = "ok"
        except Exception as e:
            print(f"    ERROR: {e}")
            out["status"] = f"error: {e}"
        rows.append(out)

    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(RESULTS_DIR, "prep_summary.csv"), index=False)
    print(f"\n{(summary['status'] == 'ok').sum()}/{len(summary)} ok. Wrote results/prep_summary.csv")


if __name__ == "__main__":
    main()
