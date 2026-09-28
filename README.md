# CYP3A4 Structure Prediction Challenge (OpenADMET)

Submissions for the structure prediction track of the OpenADMET CYP3A4 challenge. The task is to predict protein–ligand complex structures of human cytochrome P450 3A4 (UniProt [P08684](https://www.uniprot.org/uniprotkb/P08684)) with a set of blinded ligands.

> Work in progress. Methods and results will be added as the challenge proceeds.

## Inputs

| File | Description |
|---|---|
| `input/CYP3A4_protein_sequence.fasta` | Full-length CYP3A4 sequence (503 aa) |
| `input/cyp-challenge-TEST-BLINDED_structures.csv` | Blinded test ligands (`Molecule_Name`, `SMILES`) |
| `input/rcsb_pdb_custom_report_20260926222041.csv` | RCSB custom report: 122 CYP3A4 PDB entries with ligands and affinities |

## Analyses

Each analysis is self-contained in `analyses/`, with its own README (question, decisions, how to rerun, results). A chronological record is kept in [`PROJECT_LOG.md`](PROJECT_LOG.md).

| # | Analysis | Status |
|---|---|---|
| 01 | [Glide score of crystal poses vs IC50](analyses/01_glide_score_vs_ic50/) | Complete — no usable correlation (rho = −0.27) |
| 02 | [Redocking pose recovery](analyses/02_redocking_pose_recovery/) | Complete — 18 % top-1 success at 2 Å |
| 03 | [Heme representation & Glide's metal term](analyses/03_heme_representation/) | Complete — activation follows no chemical logic |
| 04 | [Does the wrong chemistry dock better?](analyses/04_wrong_chemistry_replication/) | Complete — no, the metal reward is unreachable |
| 05 | [Full-set redocking from SMILES (Glide)](analyses/05_full_pdb_redocking/) | 97/100 docked |
| 06 | [Redocking with ICM-Pro](analyses/06_icm_redocking/) | Running |
| 07 | [Blinded ligands vs PDB chemistry](analyses/07_test_vs_pdb_chemistry/) | Complete — sets barely overlap; benchmark likely pessimistic |

## Tools

- [AlphaFold3](https://github.com/google-deepmind/alphafold3): co-folding
- [Boltz-2](https://github.com/jwohlwend/boltz): co-folding
- [Schrödinger Glide](https://www.schrodinger.com/): docking
- [Molsoft ICM-Pro](https://www.molsoft.com/icm_pro.html): docking

## Related work

The approach builds on my submissions to the [OpenADMET PXR structure challenge](https://github.com/surajugrani/pxr-structure-challenge).
