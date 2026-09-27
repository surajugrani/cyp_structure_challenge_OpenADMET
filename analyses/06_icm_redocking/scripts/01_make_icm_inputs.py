"""
Analysis 06 step 1: build ICM inputs from the PRE-PrepWizard cleaned structures.

ICM prepares the receptor itself (`_dockBatch` calls `convertObject` with hydrogen
optimisation, His optimisation and **hetero charge assignment** — the last being what
assigns charges to the heme). Feeding it the PrepWizard-prepared structure would be double
preparation and would inherit the heme quirks analyses 03/04 documented, so ICM starts
from the *cleaned* crystal structure instead:

    <PDB>.pdb          raw download
    <PDB>_clean.mae    chain + heme + ligand, waters/cryoprotectants stripped   <- HERE
    <PDB>_prepped.mae  PrepWizard output (used by Glide, NOT by ICM)

That keeps the structural starting point identical between the two pipelines — same chain,
same ligand copy (closest to the heme Fe), same het groups removed — while letting each
program do its own chemistry.

The ligand is deliberately NOT prepared by ICM: analysis 06 docks the **same LigPrep
output** that analysis 05 gave Glide, so ligand preparation is held constant and any
difference is attributable to the docking engine.

This script writes, per structure:
    <PDB>_rec.pdb        receptor: protein + heme, crystal ligand removed
    <PDB>_site_lig.sdf   crystal ligand, used by ICM only to define the binding site
    <PDB>_ligand.sdf     the LigPrep states from analysis 05 (what actually gets docked)
    <PDB>_crystal.sdf    crystal ligand again, kept as the RMSD reference

Usage (light - no docking, no licence):
    $SCHRODINGER/run python3 scripts/01_make_icm_inputs.py [--only 4D6Z,6MA7]
"""

import argparse
import os
import sys

import pandas as pd
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A01 = os.path.join(ANALYSES, "01_glide_score_vs_ic50")
A05 = os.path.join(ANALYSES, "05_full_pdb_redocking")
STRUCT_DIR = os.path.join(HERE, "structures")
RESULTS = os.path.join(HERE, "results")

sys.path.insert(0, os.path.join(A01, "scripts"))
from importlib import import_module
_p = import_module("03_prep_structures")
clean_structure = _p.clean_structure


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    args = ap.parse_args()

    s = pd.read_csv(os.path.join(A05, "data", "structure_set.csv"))
    if args.only:
        s = s[s["pdb_id"].isin(args.only.split(","))]
    os.makedirs(RESULTS, exist_ok=True)

    rows = []
    for _, row in s.iterrows():
        pdb, comp, chain = row["pdb_id"], row["comp_id"], row["chain"]
        src = os.path.join(A05, "structures", pdb)
        wdir = os.path.join(STRUCT_DIR, pdb)
        os.makedirs(wdir, exist_ok=True)
        rec = {"pdb_id": pdb, "comp_id": comp}
        try:
            raw = os.path.join(src, f"{pdb}.pdb")
            ligprep = os.path.join(src, f"{pdb}_ligprep.maegz")
            if not os.path.exists(raw):
                raise RuntimeError("raw PDB missing - run analysis 05 first")
            if not os.path.exists(ligprep):
                raise RuntimeError("LigPrep output missing - analysis 05 not done here yet")

            # Cleaned crystal structure (identical selection logic to analyses 01/05)
            clean, info = clean_structure(raw, chain, comp)
            lig_idx = [a.index for a in clean.atom
                       if a.chain == info["chain"] and a.resnum == info["score_resnum"]
                       and a.pdbres.strip() == comp]
            if not lig_idx:
                raise RuntimeError("ligand not located in cleaned structure")

            # Crystal ligand: site definition for ICM, and the RMSD reference
            lig = clean.extract(lig_idx, copy_props=True)
            lig.write(os.path.join(wdir, f"{pdb}_site_lig.sdf"))
            lig.write(os.path.join(wdir, f"{pdb}_crystal.sdf"))

            # Receptor: everything except the ligand, as PDB for ICM
            receptor = clean.copy()
            receptor.deleteAtoms(lig_idx)
            receptor.write(os.path.join(wdir, f"{pdb}_rec.pdb"))

            # Ligand to dock: the same LigPrep states Glide received (option b)
            n_states = 0
            with structure.StructureWriter(os.path.join(wdir, f"{pdb}_ligand.sdf")) as w:
                for st in structure.StructureReader(ligprep):
                    w.append(st)
                    n_states += 1

            rec.update({"chain": info["chain"], "score_resnum": info["score_resnum"],
                        "fe_dist": info["score_min_dist_to_fe"],
                        "n_lig_atoms": len(lig_idx), "n_ligprep_states": n_states,
                        "n_rec_atoms": receptor.atom_total, "status": "ok"})
            print(f"{pdb} {comp}: receptor {receptor.atom_total} atoms, "
                  f"{n_states} LigPrep state(s)")
        except Exception as e:
            rec["status"] = f"error: {e}"
            print(f"{pdb} {comp}: ERROR {e}")
        rows.append(rec)

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "icm_inputs.csv")
    if os.path.exists(out) and args.only:
        old = pd.read_csv(out)
        df = pd.concat([old[~old["pdb_id"].isin(df["pdb_id"])], df], ignore_index=True)
    df.to_csv(out, index=False)
    print(f"\n{(df['status'] == 'ok').sum()}/{len(df)} ready. Wrote {out}")


if __name__ == "__main__":
    main()
