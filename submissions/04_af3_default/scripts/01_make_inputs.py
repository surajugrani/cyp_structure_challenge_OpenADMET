"""Submission 04, step 1: build AF3 inputs for the 20 blinded ligands.

Follows the PXR recipe (data/raw/sample_af3_json.json) rather than inventing one:
protein chain with `unpairedMsaPath`, `pairedMsa: ""`, plus the ligand.

Two CYP-specific changes:
  * the haem is added as a third chain, CCD `HEM` (PXR was protein + ligand only)
  * `templates: []` rather than null, so the run is a clean no-template baseline and
    AF3's data pipeline can be skipped entirely

MSA: all 20 ligands share one CYP3A4 sequence, so one MSA serves all of them. Rather than
run AF3's search 20 times - or even once - the ColabFold MSA that Boltz already fetched in
submission 01 is reused (verified: its first sequence is the 503 aa query, exactly).
ColabFold returns uniref and bfd/mgnify/metaeuk/smag separately; they are concatenated
here the way ColabFold itself merges them, dropping the duplicated query from the second.

Usage:
    python3 scripts/01_make_inputs.py
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
INPUTS = os.path.join(HERE, "inputs")
MSA_OUT = os.path.join(HERE, "inputs", "cyp3a4_colabfold.a3m")
BOLTZ_MSA = os.path.join(ROOT, "submissions", "01_boltz2_default", "results",
                         "boltz_results_inputs", "msa", "OCNT-2312554_unpaired_tmp_env")


def read_a3m(path):
    """Read an MMseqs2/ColabFold a3m, stripping the trailing NUL byte it writes.

    MMseqs2 terminates these files with \x00. AF3 reads that as a residue and dies with
    "Error extracting MSA or deletion features: Unknown residues in MSA" - which is what
    killed job 47834604. Boltz tolerates it, so it only shows up when reusing the file.
    """
    txt = open(path).read().replace("\x00", "")
    return [l for l in txt.split("\n") if l.strip()]


def build_msa(seq):
    uni = read_a3m(os.path.join(BOLTZ_MSA, "uniref.a3m"))
    bfd = read_a3m(os.path.join(BOLTZ_MSA, "bfd.mgnify30.metaeuk30.smag30.a3m"))
    assert uni[1].replace("-", "").upper() == seq, "uniref a3m query != CYP3A4 sequence"
    assert bfd[1].replace("-", "").upper() == seq, "bfd a3m query != CYP3A4 sequence"
    merged = uni + bfd[2:]          # drop the duplicated query header+sequence
    with open(MSA_OUT, "w") as f:
        f.write("\n".join(merged) + "\n")
    n = sum(1 for l in merged if l.startswith(">"))
    print(f"merged MSA: {n} sequences -> {MSA_OUT}")
    return MSA_OUT


def main():
    os.makedirs(INPUTS, exist_ok=True)
    seq = "".join(l.strip() for l in
                  open(os.path.join(ROOT, "input", "CYP3A4_protein_sequence.fasta"))
                  if not l.startswith(">"))
    msa = build_msa(seq)

    rows = [r for r in csv.DictReader(
        open(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
        if r["Molecule_Name"].strip()]
    for r in rows:
        name, smi = r["Molecule_Name"].strip(), r["SMILES"].strip()
        obj = {
            "name": name,
            "modelSeeds": [1],
            "sequences": [
                {"protein": {"id": "A", "sequence": seq,
                             "unpairedMsaPath": msa, "pairedMsa": "", "templates": []}},
                {"ligand": {"id": "B", "ccdCodes": ["HEM"]}},
                {"ligand": {"id": "C", "smiles": smi}},
            ],
            "dialect": "alphafold3",
            "version": 3,
        }
        json.dump(obj, open(os.path.join(INPUTS, f"{name}.json"), "w"), indent=2)
    print(f"wrote {len(rows)} AF3 inputs to {INPUTS}")
    assert len(rows) == 20


if __name__ == "__main__":
    main()
