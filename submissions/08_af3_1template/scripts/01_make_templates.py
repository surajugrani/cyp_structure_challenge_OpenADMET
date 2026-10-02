"""Submission 08, step 1: build ONE single-chain AF3 template from 38FW (ketoconazole).

Submission 07 used two templates, one apo (38FU) and one holo (38FY), and scored worse
than no templates at all (LDDT-PLI 0.3853 vs 0.4279). That is NOT evidence against
templates: PXR's template runs mostly beat its plain AF3 run (0.5087, 0.5104, 0.5132 and
0.5285 against 0.5055), and its single best AF3 submission used one. The distinguishing
feature of 07 is that it used TWO - and an apo and a holo pocket differ precisely in the
F/G loop region, so the pair plausibly gave AF3 contradictory restraints there.

So: one template, and a holo one.

38FW (ketoconazole) rather than 38FY (vardenafil), on binding mode rather than map
quality. 12 of the 20 blinded ligands coordinate the haem iron through an imidazole,
pyridine or triazole nitrogen at 2.06-2.32 A - classic type II. Ketoconazole is exactly
that, an imidazole nitrogen on the iron. Vardenafil is the one cryo-EM ligand the preprint
reports as NEITHER type I nor type II (no red-shifted Soret peak), so its pocket is the
least representative of our set despite the best FSC and Q score. A template supplies only
protein coordinates, so what transfers is the pocket conformation its ligand induced,
which argues for matching binding mode over map quality.

Ritonavir (38FV) is excluded at 721 Da against 240-412 Da for the blinded set; azamulin
(38FX) because it is the only structure leaving F/G loop residues 214-217 unmodelled.

Same three AF3 requirements handled as in submission 07, each verified against the AF3
source: single polymer chain (38FW is a C3 trimer), a preserved release date from
_pdbx_audit_revision_history, and templateIndices taken from a real pairwise alignment
because the sequences do not line up arithmetically.

Run with: module load alphafold/3.0; python 01_make_templates.py
"""
import json
import os

from alphafold3 import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
TDIR = os.path.join(HERE, "templates")
TEMPLATES = ["38FW"]
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
