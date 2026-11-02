"""Write submission 38.s structures from the PrepWizard-relaxed complexes.

01_build.py grafted the selected ICM poses into submission 29's untouched receptor and
PoseBusters failed 14 of 20 on `minimum_distance_to_protein`, with minimum ligand-protein
heavy-atom distances of 1.38-2.57 A. The relaxation in analysis 36 was not only fixing
hydrogens: the receptor's heavy atoms moved to make room for the ICM pose. Putting the
pose back into the unrelaxed receptor takes that room away.

So ship the receptor the pose was actually scored against. The cost is explicit:
submission 29's LDDT-LP of 0.8001 no longer carries over, and LDDT-PLI and LDDT-LP will
move together, so a change in the score cannot be attributed to the ligand alone.

Rather than converting the .mae wholesale, coordinates are mapped back onto submission
29's PDB by (chain, residue number, insertion code, atom name). That keeps the file
layout, numbering and chain assignment byte-identical in structure to every submission
this project has shipped, and silently drops what PrepWizard added and the challenge
does not want: hydrogens, and the ACE/NMA caps from `-captermini`.

The LIG records are handled separately. PrepWizard does not preserve the ligand's chain
and atom names, so matching it by (chain, resnum, name) silently leaves submission 29's
original ligand in place - the first run of this script did exactly that, relaxing every
receptor and transferring no pose at all. The ligand is therefore taken from the relaxed
complex in file order and checked against the PDB's element sequence before any
coordinate is written, the same guard 01_build.py uses.

The ligand shipped is the PrepWizard-relaxed ICM pose, which is the object that was
scored, not the raw ICM pose. The two differ by a median of 0.35 A.

Any receptor atom with no counterpart keeps its original coordinate and is counted in
`n_unmatched`. A nonzero count means part of the receptor was NOT relaxed and the clash
may survive there, so it is printed per ligand rather than summarised away.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORK = os.path.join(ROOT, "analyses", "36_prepped_icm_poses", "work")
SEL = os.path.join(HERE, "results", "selection.csv")
SRC = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three", "submission")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
CAPS = {"ACE", "NMA", "NME"}


def main():
    from schrodinger import structure
    picks = {r["ligand"]: int(r["pose_rank"]) for r in csv.DictReader(open(SEL))}
    rows = []
    for lig in sorted(picks):
        mae = os.path.join(WORK, lig, f"pose{picks[lig]:02d}", "complex.mae")
        if not os.path.exists(mae):
            raise SystemExit(f"{lig}: {mae} missing")
        st = next(iter(structure.StructureReader(mae)))

        coords, pose = {}, []
        for a in st.atom:
            if a.atomic_number <= 1 or a.pdbres.strip() in CAPS:
                continue
            if a.pdbres.strip() == "LIG":
                pose.append((a.element.upper(), tuple(a.xyz)))
                continue
            key = (a.chain.strip(), a.resnum, a.inscode.strip(), a.pdbname.strip())
            coords[key] = tuple(a.xyz)

        src = open(os.path.join(SRC, f"{lig}.pdb")).read().split("\n")
        lig_i = [i for i, l in enumerate(src)
                 if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "LIG"]
        if len(lig_i) != len(pose):
            raise SystemExit(f"{lig}: {len(lig_i)} LIG atoms in the PDB vs "
                             f"{len(pose)} in the relaxed complex")
        pdb_el = [(src[i][76:78].strip() or src[i][12:16].strip()[0]).upper()
                  for i in lig_i]
        if pdb_el != [e for e, _ in pose]:
            raise SystemExit(f"{lig}: element sequence differs between the PDB LIG "
                             f"records and the relaxed pose - refusing to write")
        for i, (_, c) in zip(lig_i, pose):
            coords[("LIG", i)] = c

        moved, unmatched = [], 0
        lig_set = set(lig_i)
        for i, l in enumerate(src):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            key = ("LIG", i) if i in lig_set else (
                l[21].strip(), int(l[22:26]), l[26].strip(), l[12:16].strip())
            if key not in coords:
                unmatched += 1
                continue
            x, y, z = coords[key]
            old = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            moved.append(math.dist(old, (x, y, z)))
            src[i] = f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}"
        open(os.path.join(OUT, f"{lig}.pdb"), "w").write("\n".join(src))

        moved.sort()
        rows.append({"ligand": lig, "pose_rank": picks[lig],
                     "n_atoms_updated": len(moved), "n_unmatched": unmatched,
                     "median_shift": round(moved[len(moved) // 2], 3),
                     "max_shift": round(moved[-1], 3)})
        print(f"{lig}: updated {len(moved)}, unmatched {unmatched}, "
              f"median shift {moved[len(moved)//2]:.3f} A, max {moved[-1]:.3f}",
              flush=True)

    with open(os.path.join(RES, "prepped_export.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    tot = sum(r["n_unmatched"] for r in rows)
    print(f"\n{len(rows)}/20 exported, {tot} unmatched atoms in total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
