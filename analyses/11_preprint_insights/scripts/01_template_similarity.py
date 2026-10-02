"""How similar are the 4 cryo-EM template ligands to the 20 blinded ligands?

The preprint identifies the F/G loop as what reshapes the CYP3A4 pocket per ligand, so a
template is useful insofar as its pocket conformation matches the one our ligand needs.
Chemical similarity is the available proxy for that.

Template ligands: ritonavir (RIT, 38FV), ketoconazole (KLN, 38FW), azamulin (MWY, 38FX),
vardenafil (VDN, 38FY). Apo is 38FU.
"""
import csv, json, os, urllib.request
import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
TEMPLATES = {"RIT": ("ritonavir", "38FV"), "KLN": ("ketoconazole", "38FW"),
             "MWY": ("azamulin", "38FX"), "VDN": ("vardenafil", "38FY")}


def ccd_smiles(code):
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}"
    d = json.load(urllib.request.urlopen(url, timeout=30))
    for desc in d.get("pdbx_chem_comp_descriptor", []):
        if desc.get("type") == "SMILES_CANONICAL" and desc.get("program") == "OpenEye OEToolkits":
            return desc["descriptor"]
    for desc in d.get("pdbx_chem_comp_descriptor", []):
        if desc.get("type") == "SMILES":
            return desc["descriptor"]
    return None


def main():
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    tmpl = {}
    for code, (name, pdb) in TEMPLATES.items():
        smi = ccd_smiles(code)
        m = Chem.MolFromSmiles(smi) if smi else None
        if m is None:
            print(f"  WARNING: could not fetch {code}")
            continue
        tmpl[code] = {"name": name, "pdb": pdb, "mol": m, "fp": gen.GetFingerprint(m),
                      "mw": Descriptors.MolWt(m)}
        print(f"  {code} {name:14} {pdb}  MW {Descriptors.MolWt(m):6.1f}")

    rows = list(csv.DictReader(open(os.path.join(
        ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv"))))
    codes = list(tmpl)
    print(f"\n{'blinded ligand':16}{'MW':>7}" + "".join(f"{tmpl[c]['name'][:9]:>11}" for c in codes)
          + f"{'best':>12}")
    out, best_counts = [], {}
    for r in rows:
        n, smi = r["Molecule_Name"].strip(), r["SMILES"].strip()
        if not n:
            continue
        m = Chem.MolFromSmiles(smi)
        fp = gen.GetFingerprint(m)
        sims = {c: DataStructs.TanimotoSimilarity(fp, tmpl[c]["fp"]) for c in codes}
        b = max(sims, key=sims.get)
        best_counts[b] = best_counts.get(b, 0) + 1
        out.append({"ligand": n, "MW": round(Descriptors.MolWt(m), 1),
                    **{tmpl[c]["name"]: round(sims[c], 3) for c in codes},
                    "best_template": tmpl[b]["name"], "best_pdb": tmpl[b]["pdb"],
                    "best_tanimoto": round(sims[b], 3)})
        print(f"{n:16}{Descriptors.MolWt(m):>7.1f}"
              + "".join(f"{sims[c]:>11.2f}" for c in codes)
              + f"{tmpl[b]['name'][:9]+' '+format(sims[b],'.2f'):>12}")

    allsim = [v for r in out for k, v in r.items()
              if k in [tmpl[c]["name"] for c in codes]]
    print(f"\nall 80 pairwise Tanimotos: max {max(allsim):.3f}, median {np.median(allsim):.3f}")
    print(f"pairs >= 0.4 (conventionally 'related'): {sum(s >= 0.4 for s in allsim)}/80")
    print(f"best-template counts: {best_counts}")
    with open(os.path.join(HERE, "results", "template_similarity.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    main()
