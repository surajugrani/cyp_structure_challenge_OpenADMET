"""Submission 07, step 2: AF3 inputs with the two cryo-EM templates.

Same as submission 04 (ColabFold MSA reused via unpairedMsaPath, haem as CCD HEM,
ligand from SMILES) except `templates` now carries 38FY_A and 38FU_A with the
alignment-derived index maps from step 1, instead of being empty.
"""
import csv, json, os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
MSA = os.path.join(ROOT, "submissions", "04_af3_default", "inputs",
                   "cyp3a4_colabfold.a3m")
MAPS = os.path.join(HERE, "results", "template_maps.json")
INPUTS = os.path.join(HERE, "inputs")


def main():
    os.makedirs(INPUTS, exist_ok=True)
    seq = "".join(l.strip() for l in
                  open(os.path.join(ROOT, "input", "CYP3A4_protein_sequence.fasta"))
                  if not l.startswith(">"))
    assert os.path.exists(MSA), "shared ColabFold MSA missing"
    maps = json.load(open(MAPS))
    templates = [{"mmcifPath": os.path.abspath(v["path"]),
                  "queryIndices": v["queryIndices"],
                  "templateIndices": v["templateIndices"]} for v in maps.values()]
    print(f"{len(templates)} templates: " +
          ", ".join(f"{k} ({len(v['queryIndices'])} pos)" for k, v in maps.items()))

    n = 0
    for r in csv.DictReader(open(os.path.join(
            ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv"))):
        name, smi = r["Molecule_Name"].strip(), r["SMILES"].strip()
        if not name:
            continue
        obj = {"name": name, "modelSeeds": [1],
               "sequences": [
                   {"protein": {"id": "A", "sequence": seq,
                                "unpairedMsaPath": os.path.abspath(MSA),
                                "pairedMsa": "", "templates": templates}},
                   {"ligand": {"id": "B", "ccdCodes": ["HEM"]}},
                   {"ligand": {"id": "C", "smiles": smi}}],
               "dialect": "alphafold3", "version": 3}
        json.dump(obj, open(os.path.join(INPUTS, f"{name}.json"), "w"))
        n += 1
    print(f"wrote {n} inputs to {INPUTS}")
    assert n == 20


if __name__ == "__main__":
    main()
