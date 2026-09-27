"""
Analysis 05 step 1: define the full CYP3A4 redocking set.

Analysis 02 measured pose recovery on the 22 structures that happened to have BindingDB
affinities. Pose prediction is the actual goal of the challenge, and it does not need
affinity data — so this widens the set to every usable CYP3A4 crystal structure in the
RCSB report, which tightens the confidence interval on the recovery rate.

Selection, from `analyses/01_glide_score_vs_ic50/data/entries.csv` (122 entries):

  drop  cofactor / cryoprotectant / buffer / ion-only het groups (JUNK below)
  drop  organometallic ligands - docking force fields have no parameters for transition
        metals. Originally this excluded Ir only; Ru was added after LigPrep silently
        dropped X8S (7KS8) and X8P (7KSA), two Ru complexes, producing no output. The
        filter is now a metal list rather than a single element.
  drop  ligands with MW < 150 (metformin in 5G5J, imidazole in 9YK4): both are genuine
        ligands, but with so few heavy atoms a 2 A RMSD cutoff is nearly free and would
        inflate the success rate
  drop  entries left with no drug-like ligand (apo / heme-only structures)

Usage (light):
    $SCHRODINGER/run python3 scripts/01_build_set.py

Output: data/structure_set.csv  (pdb_id, comp_id, chain, resolution, has_affinity, ...)
"""

import os

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A01 = os.path.join(os.path.dirname(HERE), "01_glide_score_vs_ic50")
DATA = os.path.join(HERE, "data")

MIN_MW = 150.0

# Transition/heavy metals that docking force fields cannot parameterise. Checked against
# the ligand formula. (Fe is deliberately NOT here - it is the heme cofactor in the
# receptor, never the ligand of interest.)
METALS = ["Ir", "Ru", "Rh", "Pd", "Pt", "Os", "Re", "Mn", "Co", "Ni", "Cu",
          "Mo", "W", "V", "Cr", "Ti", "Au", "Ag", "Cd", "Hg"]
METAL_RE = r"\b(" + "|".join(METALS) + r")\d*\b"

# Cofactor, cryoprotectant, buffer, precipitant and ion het codes - never the ligand of
# interest for this target.
JUNK = {
    "HEM", "GOL", "EDO", "DMS", "SO4", "PO4", "CL", "NA", "MG", "K", "ZN", "CA",
    "ACT", "MPD", "TRS", "EPE", "PEG", "PGE", "P6G", "1PE", "2PE", "PG4", "TEG",
    "FMT", "ACY", "NO3", "IOD", "BR", "CIT", "MES", "BME", "DTT", "SCN", "AZI",
    "URE", "NH4", "BCT", "OXL", "TAR", "MLI", "GLY", "ALA", "SER", "GLN", "LEU",
}


def main():
    e = pd.read_csv(os.path.join(A01, "data", "entries.csv"))
    cur = pd.read_csv(os.path.join(A01, "data", "curated_set.csv"))
    os.makedirs(DATA, exist_ok=True)

    n0 = e["pdb_id"].nunique()
    c = e[~e["comp_id"].isin(JUNK)]
    metal_hits = c[c["formula"].str.contains(METAL_RE, na=False, regex=True)]
    c = c[~c["formula"].str.contains(METAL_RE, na=False, regex=True)]
    dropped_small = c[c["formula_weight"] < MIN_MW]
    c = c[c["formula_weight"] >= MIN_MW]

    print(f"{n0} entries in the RCSB report")
    if len(metal_hits):
        print(f"  dropped organometallic ligands: "
              f"{', '.join(f'{r.comp_id} in {r.pdb_id} ({r.formula})' for r in metal_hits.itertuples())}")
    print(f"  dropped small ligands (MW < {MIN_MW:g}): "
          f"{', '.join(f'{r.comp_id} in {r.pdb_id} ({r.formula_weight:.0f} Da)' for r in dropped_small.itertuples())}")
    print(f"  {c['pdb_id'].nunique()} entries retain a drug-like ligand")

    multi = c.groupby("pdb_id").size()
    if (multi > 1).any():
        print(f"  WARNING: {(multi > 1).sum()} entries have >1 candidate ligand: "
              f"{list(multi[multi > 1].index)}")

    # Ligand chain: lowest chain ID where the ligand appears (chain A where present).
    # The prep step then picks the copy closest to the heme Fe within that chain.
    def pick_chain(instances):
        chains = sorted({x.split(":")[0] for x in str(instances).split(";")})
        return chains[0]

    c = c.copy()
    c["chain"] = c["instances"].apply(pick_chain)
    c["has_affinity"] = c["pdb_id"].isin(set(cur["pdb_id"]))

    out = c[["pdb_id", "comp_id", "chain", "resolution", "method", "name",
             "formula", "formula_weight", "n_instances", "instances", "has_affinity"]]
    out = out.sort_values(["has_affinity", "resolution"], ascending=[False, True])
    out.to_csv(os.path.join(DATA, "structure_set.csv"), index=False)

    print(f"\nFinal set: {len(out)} structures "
          f"({out.has_affinity.sum()} with affinity data, {(~out.has_affinity).sum()} new)")
    print(f"  resolution {out.resolution.min():.2f}-{out.resolution.max():.2f} A "
          f"(median {out.resolution.median():.2f})")
    print(f"  ligand MW {out.formula_weight.min():.0f}-{out.formula_weight.max():.0f} Da "
          f"(median {out.formula_weight.median():.0f})")
    print(f"\nWrote {os.path.join(DATA, 'structure_set.csv')}")


if __name__ == "__main__":
    main()
