# CYP3A4 Structure Prediction Challenge (OpenADMET)

Submissions for the structure prediction track of the OpenADMET CYP3A4 challenge. The task is to predict protein–ligand complex structures of human cytochrome P450 3A4 (UniProt [P08684](https://www.uniprot.org/uniprotkb/P08684)) with a set of blinded ligands.

> Work in progress. Methods and results will be added as the challenge proceeds.

## Inputs

| File | Description |
|---|---|
| `input/CYP3A4_protein_sequence.fasta` | Full-length CYP3A4 sequence (503 aa) |
| `input/cyp-challenge-TEST-BLINDED_structures.csv` | Blinded test ligands (`Molecule_Name`, `SMILES`) |

## Tools

- [AlphaFold3](https://github.com/google-deepmind/alphafold3): co-folding
- [Boltz-2](https://github.com/jwohlwend/boltz): co-folding
- [Schrödinger Glide](https://www.schrodinger.com/): docking
- [Molsoft ICM-Pro](https://www.molsoft.com/icm_pro.html): docking

## Related work

The approach builds on my submissions to the [OpenADMET PXR structure challenge](https://github.com/surajugrani/pxr-structure-challenge).
