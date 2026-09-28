"""
Analysis 07: how similar are the 20 blinded challenge ligands to the CYP3A4 PDB ligands?

Analyses 05/06 measure docking performance on ~100 CYP3A4 crystal structures. That is only
a useful guide to the blind set if the two ligand sets occupy similar chemical space. This
quantifies the gap.

Four comparisons:

  1. Tanimoto similarity (Morgan r=2, 2048 bits). For each blinded ligand, the maximum
     similarity to any PDB ligand — i.e. "is there anything like this in the benchmark?"
     Conventional reference points: >0.7 very similar, 0.4-0.7 related, <0.4 distinct.

  2. Physicochemical properties: MW, cLogP, TPSA, rotatable bonds, HBD/HBA, heavy atoms,
     aromatic rings, fraction sp3. Distribution overlap, with Mann-Whitney per property.

  3. Bemis-Murcko scaffolds — how many blinded scaffolds appear in the PDB set at all.

  4. **Fe-coordinating motifs.** Analyses 01-04 found the benchmark is dominated by type-II
     inhibitors that coordinate the heme iron (16/20 of the affinity subset at 1.9-2.4 A),
     and that this is precisely what Glide's scoring mishandles. Whether the blinded
     ligands carry the same motifs (imidazole, triazole, pyridine, thiazole, ...) decides
     whether our failure mode transfers to them at all. This is the most decision-relevant
     of the four.

Usage (light - RDKit only, no licence, no docking):
    $SCHRODINGER/run python3 scripts/01_compare_chemistry.py

Outputs:
    results/similarity.csv        per-blinded-ligand nearest neighbour in the PDB set
    results/properties.csv        properties for every ligand in both sets
    results/property_stats.csv    per-property comparison with Mann-Whitney p
    results/motifs.csv            Fe-coordinating motif counts per set
    results/chemistry_comparison.png
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from scipy import stats

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A05 = os.path.join(os.path.dirname(HERE), "05_full_pdb_redocking")
RESULTS = os.path.join(HERE, "results")

# Groups that can ligate a haem iron in a type-II CYP complex. The SMARTS match the
# coordinating heteroatom's ring, not a specific substitution pattern.
FE_MOTIFS = {
    "imidazole": "c1cnc[nH0,nH]1",
    "triazole": "c1ncn[nH0,nH]1",
    "tetrazole": "c1nnn[nH0,nH]1",
    "pyridine": "n1ccccc1",
    "pyrimidine": "n1cnccc1",
    "thiazole": "c1scnc1",
    "oxazole": "c1ocnc1",
    "nitrile": "[NX1]#[CX2]",
    "primary_amine": "[NX3;H2;!$(NC=O)]",
}


def props(mol):
    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotB": Descriptors.NumRotatableBonds(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "HeavyAtoms": mol.GetNumHeavyAtoms(),
        "AromRings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "FracCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
    }


def motif_hits(mol):
    out = {}
    for name, smarts in FE_MOTIFS.items():
        patt = Chem.MolFromSmarts(smarts)
        out[name] = bool(patt) and mol.HasSubstructMatch(patt)
    out["any_Fe_motif"] = any(v for k, v in out.items())
    return out


def main():
    os.makedirs(RESULTS, exist_ok=True)

    test = pd.read_csv(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv"))
    pdb_smiles = json.load(open(os.path.join(A05, "data", "ligand_smiles.json")))

    sets = {
        "blinded": [(r.Molecule_Name, r.SMILES) for r in test.itertuples()],
        "pdb": [(k, v["smiles_stereo"]) for k, v in pdb_smiles.items()],
    }
    mols, rows = {}, []
    for label, items in sets.items():
        mols[label] = []
        for name, smi in items:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                print(f"  WARNING: unparseable SMILES for {name}")
                continue
            mols[label].append((name, m))
            rows.append({"set": label, "name": name, "smiles": smi,
                         **props(m), **motif_hits(m),
                         "scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=m)})
    allp = pd.DataFrame(rows)
    allp.to_csv(os.path.join(RESULTS, "properties.csv"), index=False)

    # ── 1. Tanimoto nearest neighbour ───────────────────────────────────────
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    fps = {k: [(n, gen.GetFingerprint(m)) for n, m in v] for k, v in mols.items()}
    sim_rows = []
    for name, fp in fps["blinded"]:
        sims = DataStructs.BulkTanimotoSimilarity(fp, [f for _, f in fps["pdb"]])
        j = int(np.argmax(sims))
        sim_rows.append({"blinded": name, "max_tanimoto": sims[j],
                         "nearest_pdb_ligand": fps["pdb"][j][0],
                         "mean_tanimoto": float(np.mean(sims))})
    sim = pd.DataFrame(sim_rows).sort_values("max_tanimoto", ascending=False)
    sim.to_csv(os.path.join(RESULTS, "similarity.csv"), index=False, float_format="%.3f")

    print("=== 1. Nearest-neighbour Tanimoto (blinded vs PDB ligands) ===")
    print(sim.to_string(index=False, float_format="%.3f"))
    mx = sim.max_tanimoto
    print(f"\n  max-similarity: median {mx.median():.3f}, range {mx.min():.3f}-{mx.max():.3f}")
    for cut, lab in ((0.7, "very similar"), (0.4, "related")):
        print(f"  >= {cut} ({lab}): {(mx >= cut).sum()}/{len(mx)}")

    # ── 2. Properties ───────────────────────────────────────────────────────
    print("\n=== 2. Physicochemical properties ===")
    pstats = []
    b = allp[allp.set == "blinded"]; p = allp[allp.set == "pdb"]
    print(f'{"property":<12}{"blinded median":>16}{"pdb median":>13}{"MWU p":>10}')
    for c in ["MW", "cLogP", "TPSA", "RotB", "HBD", "HBA", "HeavyAtoms", "AromRings", "FracCSP3"]:
        u = stats.mannwhitneyu(b[c], p[c])
        pstats.append({"property": c, "blinded_median": b[c].median(),
                       "pdb_median": p[c].median(), "mwu_p": u.pvalue})
        flag = " *" if u.pvalue < 0.05 else ""
        print(f"{c:<12}{b[c].median():>16.2f}{p[c].median():>13.2f}{u.pvalue:>10.3g}{flag}")
    pd.DataFrame(pstats).to_csv(os.path.join(RESULTS, "property_stats.csv"),
                                index=False, float_format="%.4f")

    # ── 3. Scaffolds ────────────────────────────────────────────────────────
    bs, ps = set(b.scaffold), set(p.scaffold)
    print(f"\n=== 3. Murcko scaffolds ===")
    print(f"  blinded: {len(bs)} unique | pdb: {len(ps)} unique | shared: {len(bs & ps)}")

    # ── 4. Fe-coordinating motifs ───────────────────────────────────────────
    print("\n=== 4. Fe-coordinating motifs (the decision-relevant one) ===")
    mrows = []
    for c in list(FE_MOTIFS) + ["any_Fe_motif"]:
        mrows.append({"motif": c, "blinded_n": int(b[c].sum()), "blinded_pct": 100 * b[c].mean(),
                      "pdb_n": int(p[c].sum()), "pdb_pct": 100 * p[c].mean()})
    md = pd.DataFrame(mrows)
    md.to_csv(os.path.join(RESULTS, "motifs.csv"), index=False, float_format="%.1f")
    print(md.to_string(index=False, float_format="%.1f"))

    # ── plot ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes[0, 0].hist(mx, bins=np.arange(0, 1.05, 0.05), color="tab:blue", edgecolor="k")
    axes[0, 0].axvline(0.4, ls="--", c="grey"); axes[0, 0].axvline(0.7, ls="--", c="grey")
    axes[0, 0].set_title("Max Tanimoto to any PDB ligand")
    axes[0, 0].set_xlabel("Tanimoto"); axes[0, 0].set_ylabel("blinded ligands")
    for ax, c in zip(axes.flat[1:], ["MW", "cLogP", "TPSA", "RotB", "HeavyAtoms"]):
        ax.hist([p[c], b[c]], bins=12, density=True, label=["PDB", "blinded"],
                color=["lightgrey", "tab:orange"], edgecolor="k")
        ax.set_title(c); ax.legend(fontsize=8)
    fig.suptitle("Blinded challenge ligands vs CYP3A4 PDB ligands", fontsize=13)
    fig.tight_layout()
    out = os.path.join(RESULTS, "chemistry_comparison.png")
    fig.savefig(out, dpi=150)
    print(f"\nWrote {out} and 4 CSVs to results/")


if __name__ == "__main__":
    main()
