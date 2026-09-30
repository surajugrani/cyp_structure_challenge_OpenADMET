"""Submission 04, step 3: pick the best AF3 sample per ligand and write submission PDBs.

Selection: AF3 writes 5 diffusion samples and ranks them by `ranking_score`. PXR found
selecting on iPTM instead was worth doing, so that is what happens here - specifically the
PROTEIN-LIGAND pair from `chain_pair_iptm[0][2]` (chain A vs ligand chain C) rather than
the global `iptm`, because global iptm is inflated by the very well-predicted protein-haem
interface (typically 0.97) and says little about the blinded ligand. AF3's own pick and
the best-PL-iptm pick disagree for 8 of the 20.

Conversion: AF3 emits mmCIF with the ligand residue named `LIG_C`; the organizers'
validator requires exactly one residue named `LIG`, so it is renamed. The haem is already
`HEM`.

Run under $SCHRODINGER/run python3 (compute node) - Schrodinger's reader handles mmCIF.

Outputs:
    results/iptm_selection.csv    which sample was chosen per ligand, and why
    submission/<name>.pdb
"""
import csv
import glob
import json
import os

from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "submission")


def main():
    os.makedirs(OUT, exist_ok=True)
    names = sorted(r["Molecule_Name"].strip() for r in csv.DictReader(
        open(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
        if r["Molecule_Name"].strip())

    rows = []
    for n in names:
        d = os.path.join(RES, n.lower())
        cands = []
        for f in sorted(glob.glob(d + "/seed-*_sample-*/*summary_confidences.json")):
            j = json.load(open(f))
            cp = j.get("chain_pair_iptm")
            cands.append({
                "dir": os.path.dirname(f),
                "tag": os.path.basename(os.path.dirname(f)),
                "pl_iptm": cp[0][2] if cp else -1,
                "iptm": j.get("iptm"),
                "ranking": j.get("ranking_score"),
            })
        af3 = max(cands, key=lambda c: c["ranking"])
        best = max(cands, key=lambda c: c["pl_iptm"])

        cif = glob.glob(best["dir"] + "/*_model.cif")[0]
        st = next(iter(structure.StructureReader(cif)))
        # AF3 names the ligand residue LIG_C in mmCIF. PDB residue names are 3 characters,
        # so Schrodinger's writer renders it UNL - which is what must be renamed, not
        # "LIG*". The validator wants exactly one residue called LIG.
        n_lig = 0
        for a in st.atom:
            rn = a.pdbres.strip().upper()
            if rn == "UNL" or rn.startswith("LIG"):
                a.pdbres = "LIG "
                n_lig += 1
        if n_lig == 0:
            raise RuntimeError(f"{n}: no ligand residue found to rename "
                               f"(saw {sorted({a.pdbres.strip() for a in st.atom})})")
        dst = os.path.join(OUT, f"{n}.pdb")
        st.write(dst)

        rows.append({"ligand": n, "chosen_sample": best["tag"],
                     "af3_own_pick": af3["tag"],
                     "picks_agree": af3["tag"] == best["tag"],
                     "chosen_pl_iptm": best["pl_iptm"], "af3_pl_iptm": af3["pl_iptm"],
                     "chosen_iptm": best["iptm"], "chosen_ranking": best["ranking"],
                     "lig_atoms_renamed": n_lig})
        note = ("" if rows[-1]["picks_agree"]
                else f", AF3 would pick {af3['tag']} at {af3['pl_iptm']:.2f}")
        print(f"  {n}: {best['tag']} (PL-iptm {best['pl_iptm']:.2f}{note}), "
              f"{n_lig} ligand atoms renamed to LIG")

    with open(os.path.join(RES, "iptm_selection.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    agree = sum(r["picks_agree"] for r in rows)
    print(f"\n{len(rows)} structures written; AF3's own pick matched best PL-iptm {agree}/{len(rows)}")


if __name__ == "__main__":
    main()
