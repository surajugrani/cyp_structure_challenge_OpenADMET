"""Submission 07, step 1: build single-chain AF3 templates from the cryo-EM structures.

Templates: 38FY (vardenafil) and 38FU (apo) - chosen on model quality, not ligand
similarity. None of the 4 cryo-EM ligands resembles the blinded set (max Tanimoto 0.218,
0/80 pairs >= 0.4), so per-ligand matching has no basis. 38FY and 38FU have the joint-best
map-model FSC (3.5 A); 38FY has the highest Q score (0.478) and 38FU the best MolProbity
(1.00) and clash score (2.24). 38FX (azamulin) is deliberately excluded: it has the worst
clash score and is the only one of the five leaving F/G loop residues 214-217 unmodelled -
and the F/G loop is the entire reason for using these templates.

Three things AF3 enforces that this script handles, each verified against the AF3 source:

 1. `get_polymer_features` raises "The structure must be filtered to a single chain" if
    the template has more than one polymer chain. Both files are C3 trimers (chains
    A/B/C), so one chain is extracted.
 2. It also raises "The structure must have a release date", read from
    `_pdbx_audit_revision_history.revision_date`. A naive chain extraction drops that
    category, so AF3's own Structure.filter()/to_mmcif() is used, which preserves it.
 3. `templateIndices` index the template's FULL sequence including unmodelled residues
    (`to_res_arrays(include_missing_residues=True)`), and the sequences do not line up
    arithmetically - 38FU chain A is 487 residues and 38FY 485, against a 503-residue
    query and a nominal Delta3-22 length of 483 (His tag, and 38FY is missing the leading
    MA). So the mapping comes from a real pairwise alignment, keeping only positions where
    the residues actually match.

Run with: module load alphafold/3.0; python 01_make_templates.py
"""
import json
import os

from alphafold3 import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
TDIR = os.path.join(HERE, "templates")
TEMPLATES = ["38FY", "38FU"]
KEEP_CHAIN = "A"


def align_map(query: str, templ: str):
    """0-based query index -> 0-based template index, from a global alignment.

    Only positions where the aligned residues are identical are kept, so a mismatch
    cannot silently place the wrong residue's coordinates onto a query position.
    """
    from Bio import pairwise2
    aln = pairwise2.align.globalms(query, templ, 2, -1, -10, -0.5,
                                   one_alignment_only=True)[0]
    qi = ti = 0
    out = {}
    for a, b in zip(aln.seqA, aln.seqB):
        if a != "-" and b != "-":
            if a == b:
                out[qi] = ti
            qi += 1
            ti += 1
        elif a != "-":
            qi += 1
        elif b != "-":
            ti += 1
    return out


def main():
    query = "".join(l.strip() for l in
                    open(os.path.join(ROOT, "input", "CYP3A4_protein_sequence.fasta"))
                    if not l.startswith(">"))
    print(f"query: {len(query)} residues")

    meta = {}
    for pid in TEMPLATES:
        src = os.path.join(TDIR, f"{pid}.cif")
        s = structure.from_mmcif(open(src).read(), fix_mse_residues=True,
                                 fix_arginines=True, include_bonds=False,
                                 include_water=False, include_other=True)
        chains = s.polymer_auth_asym_id_to_label_asym_id()
        print(f"\n{pid}: {len(chains)} polymer chains {list(chains)}, "
              f"release_date={s.release_date}")

        one = s.filter(chain_id=[c for c in s.chains
                                 if c == KEEP_CHAIN or c.startswith(KEEP_CHAIN)][:1])
        nchain = len(one.polymer_auth_asym_id_to_label_asym_id())
        if nchain != 1:
            raise RuntimeError(f"{pid}: filtered to {nchain} polymer chains, need 1")
        if one.release_date is None:
            raise RuntimeError(f"{pid}: release date lost during extraction")

        dst = os.path.join(TDIR, f"{pid}_{KEEP_CHAIN}.cif")
        with open(dst, "w") as f:
            f.write(one.to_mmcif())

        # re-read what we wrote - this is what AF3 will actually see
        chk = structure.from_mmcif(open(dst).read(), fix_mse_residues=True,
                                   fix_arginines=True, include_bonds=False,
                                   include_water=False, include_other=True)
        m = chk.polymer_auth_asym_id_to_label_asym_id()
        if len(m) != 1:
            raise RuntimeError(f"{pid}: written file has {len(m)} polymer chains")
        if chk.release_date is None:
            raise RuntimeError(f"{pid}: written file has no release date")
        label = next(iter(m.values()))
        tseq = chk.chain_single_letter_sequence()[label]
        mp = align_map(query, tseq)
        print(f"  wrote {os.path.basename(dst)}: 1 chain, {len(tseq)} residues, "
              f"release_date={chk.release_date}, {len(mp)} mapped positions "
              f"(query {min(mp)}-{max(mp)})")
        meta[pid] = {"path": dst, "n_mapped": len(mp),
                     "queryIndices": sorted(mp), "templateIndices": [mp[k] for k in sorted(mp)]}

    json.dump({k: {"path": v["path"], "queryIndices": v["queryIndices"],
                   "templateIndices": v["templateIndices"]} for k, v in meta.items()},
              open(os.path.join(HERE, "results", "template_maps.json"), "w"))
    print(f"\nwrote {HERE}/results/template_maps.json")


if __name__ == "__main__":
    main()
