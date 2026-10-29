"""Submission 34 - submission 29's structures with Glide-minimised ligand poses.

Refinement, not selection. Submission 29 picked, for each of the 20 blinded ligands,
whichever of the three ESMFold2 runs scored best under Glide `inplace`. That pick is kept
exactly: same receptor, same run, same ligand. The only change is the ligand's
coordinates, replaced by the pose Glide produced under `DOCKING_METHOD mininplace` -
a restrained minimisation of the co-folded pose inside the receptor's own field.

Every attempt so far has rearranged WHICH structure goes where. This is the first that
changes coordinates. The case for it is that LDDT-PLI scores interface contact quality,
and minimisation is the one operation that relieves clashes and tightens contacts.

The case against, recorded before the result is known: analysis 28 found `mininplace`
gives the WORSE gscore for 38 of 60 of these same structures, median +0.177, which says
the minimiser moves the pose somewhere Glide likes less. And post-processing a good
submission is 0 for 6 in this project.

The poses are already on disk from analysis 28 - no new docking, no new scoring.
"""
import csv, os, shutil, sys
from schrodinger import structure

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
HERE = os.path.join(ROOT, "submissions", "34_glide_minimised")
SRC = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three")
A28 = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "work")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")


def lig_coords(path):
    """Heavy-atom coordinates of the ligand entry in a Glide pose-viewer file."""
    sts = list(structure.StructureReader(path))
    if len(sts) < 2:
        raise RuntimeError(f"{path}: {len(sts)} entries, expected receptor + pose")
    pose = sts[-1]
    out = []
    for a in pose.atom:
        if a.atomic_number > 1:
            out.append((a.pdbname.strip(), a.element, tuple(a.xyz)))
    return out


def main():
    picks = {r["ligand"]: r["picked_set"]
             for r in csv.DictReader(open(os.path.join(SRC, "results", "selection.csv")))}
    rows, bad = [], []
    for lig in sorted(picks):
        sset = picks[lig]
        pv = os.path.join(A28, sset, lig, "mininplace", "score_pv.mae")
        src_pdb = os.path.join(SRC, "submission", f"{lig}.pdb")
        if not os.path.exists(pv):
            bad.append((lig, "no mininplace pose")); continue

        new = lig_coords(pv)
        lines = open(src_pdb).read().splitlines(True)
        idx = [i for i, l in enumerate(lines)
               if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "LIG"]
        if len(idx) != len(new):
            bad.append((lig, f"{len(idx)} PDB atoms vs {len(new)} pose atoms")); continue

        # Prefer matching on atom name; fall back to file order if the element
        # sequence agrees, which is how Glide writes a pose it did not re-order.
        byname = {n: xyz for n, _, xyz in new}
        names = [lines[i][12:16].strip() for i in idx]
        if len(byname) == len(new) and all(n in byname for n in names):
            match, how = [byname[n] for n in names], "name"
        else:
            pdb_el = [lines[i][76:78].strip().upper() for i in idx]
            pose_el = [e.upper() for _, e, _ in new]
            if pdb_el != pose_el:
                bad.append((lig, "atom names and element order both disagree")); continue
            match, how = [xyz for _, _, xyz in new], "order"

        sq = 0.0
        for i, (x, y, z) in zip(idx, match):
            ox, oy, oz = (float(lines[i][30:38]), float(lines[i][38:46]),
                          float(lines[i][46:54]))
            sq += (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2
            lines[i] = f"{lines[i][:30]}{x:8.3f}{y:8.3f}{z:8.3f}{lines[i][54:]}"
        rmsd = (sq / len(idx)) ** 0.5

        open(os.path.join(OUT, f"{lig}.pdb"), "w").writelines(lines)
        rows.append({"ligand": lig, "picked_set": sset, "n_atoms": len(idx),
                     "matched_by": how, "displacement_rmsd": round(rmsd, 3)})
        print(f"{lig:16} {sset:16} {len(idx):3} atoms  {how:5}  moved {rmsd:5.3f} A",
              flush=True)

    with open(os.path.join(RES, "displacement.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, ["ligand", "picked_set", "n_atoms", "matched_by",
                                "displacement_rmsd"])
        w.writeheader(); w.writerows(rows)

    if rows:
        d = sorted(r["displacement_rmsd"] for r in rows)
        n = len(d)
        med = d[n // 2] if n % 2 else (d[n // 2 - 1] + d[n // 2]) / 2
        print(f"\n{n}/20 rebuilt. ligand displacement: median {med:.3f} A, "
              f"min {d[0]:.3f}, max {d[-1]:.3f}")
    for lig, why in bad:
        print(f"FAILED {lig}: {why}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
