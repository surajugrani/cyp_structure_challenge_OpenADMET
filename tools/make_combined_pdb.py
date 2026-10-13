"""Concatenate a submission's 20 PDBs into one multi-MODEL file for PyMOL.

PyMOL loads MODEL/ENDMDL records as states of a single object, so the whole submission can
be stepped through or shown at once to spot anything odd - a ligand outside the pocket, a
broken haem, a mangled chain.

    python3 tools/make_combined_pdb.py submissions/12_icm_redock_af3
    python3 tools/make_combined_pdb.py --all

Writes <submission>/visualize/all_models.pdb plus a model_index.txt mapping state number
to ligand name (PyMOL states are numbered, not named, so the index is how you tell which
state is which). Both are git-ignored - regenerable from the shipped PDBs.

Useful once loaded:
    load all_models.pdb
    split_states all_models          # one object per complex, named all_models_0001...
    set all_states, on               # or show every state superposed at once
"""
import argparse, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def combine(subdir):
    src = os.path.join(subdir, "submission")
    if not os.path.isdir(src):
        return None
    pdbs = sorted(f for f in os.listdir(src) if f.endswith(".pdb"))
    if not pdbs:
        return None
    out = os.path.join(subdir, "visualize")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "all_models.pdb")
    idx = []
    with open(path, "w") as g:
        for i, f in enumerate(pdbs, 1):
            name = f[:-4]
            idx.append(f"{i:4d}  {name}")
            g.write(f"MODEL     {i:4d}\n")
            g.write(f"REMARK   1 {name}\n")
            for l in open(os.path.join(src, f)):
                if l.startswith(("ATOM", "HETATM", "TER")):
                    g.write(l)
            g.write("ENDMDL\n")
        g.write("END\n")
    with open(os.path.join(out, "model_index.txt"), "w") as g:
        g.write("state  ligand\n")
        g.write("\n".join(idx) + "\n")
    return path, len(pdbs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("subdir", nargs="?")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    targets = []
    if a.all:
        base = os.path.join(HERE, "submissions")
        targets = [os.path.join(base, d) for d in sorted(os.listdir(base))
                   if os.path.isdir(os.path.join(base, d, "submission"))]
    elif a.subdir:
        targets = [a.subdir]
    else:
        ap.error("give a submission directory or --all")
    for t in targets:
        r = combine(t)
        if r:
            print(f"{os.path.basename(t):28} {r[1]:2d} models -> {r[0]}")
        else:
            print(f"{os.path.basename(t):28} no submission/*.pdb, skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
