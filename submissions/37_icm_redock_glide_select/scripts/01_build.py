"""Submission 37: ICM redocked poses, selected by Glide `inplace` on a re-protonated
receptor, grafted into submission 29's structures.

Route. Analysis 33 redocked all 20 ligands into submission 29's receptors with ICM,
10 poses each. Analysis 34 tried to score those 200 poses with Glide against the
untouched receptor and 177 came back at the 10000.0 rejection ceiling - the receptor's
hydrogen network had been optimised around the co-folded ligand and overlapped the ICM
pose at 0.73-1.6 A. Analysis 36 rebuilt each complex, re-ran PrepWizard with `-rehtreat`
at `-rmsd 0.6` so the hydrogens relax around the ICM pose, regenerated the grid per pose,
and scored `inplace`. 167 of 200 came back usable and all 20 ligands have at least one.
This ships the lowest-gscore pose per ligand.

WHAT CHANGES AND WHAT DOES NOT. Only the ligand coordinates. The receptor shipped is
submission 29's, byte-for-byte, NOT the PrepWizard-relaxed one the pose was scored
against. The relaxation was a means of getting a usable score, not a prediction: it moved
ligand heavy atoms a median of 0.35 A and was driven by the ICM pose itself. Shipping it
would change LDDT-LP and LDDT-PLI together and leave the result unreadable. Submission 29
scored LDDT-LP 0.8001; holding the receptor fixed keeps that known and puts the whole
difference on ligand placement.

HOW FAR THIS MOVES. Median ligand RMSD to submission 29 is 5.17 A, max 8.28, and only
2 of 20 stay within 2 A. This is not a refinement of the project's best submission - it
is a near-complete replacement of its ligand placements, decided by a Glide gscore
computed against a receptor that is not the one being shipped. The Fe-to-ligand distances
also tighten, from a 1.90-6.84 A spread to 2.14-3.65 A.

Both of those are DIFFERENCES, not evidence. Nothing here says the ICM poses are closer
to the truth than submission 29's; submission 29 is the best score this project has and
these poses have never been measured against anything real. The leaderboard decides.

The ICM poses carry hydrogens and submission 29's LIG records do not, so only heavy atoms
transfer. Element sequence is checked per ligand before any coordinate is written, and a
mismatch aborts the build rather than shipping a scrambled ligand.

NOT copy-only: the ligands move, so the full gate runs - 02_validate.py.
"""
import csv
import math
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
SEL = os.path.join(ROOT, "analyses", "36_prepped_icm_poses", "results",
                   "best_pose_per_ligand.csv")
A33 = os.path.join(ROOT, "analyses", "33_icm_redock_sub29", "structures")
SRC = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three", "submission")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_icm_redock_glide_select.zip")


def sdf_mols(path):
    """Every molecule in an SDF, as [(element, (x, y, z)), ...] - hydrogens included."""
    out = []
    for block in open(path).read().split("$$$$"):
        lines = block.strip("\n").split("\n")
        if len(lines) < 4:
            continue
        try:
            n = int(lines[3][0:3])
        except ValueError:
            continue
        out.append([(l[31:34].strip(), (float(l[0:10]), float(l[10:20]),
                                        float(l[20:30]))) for l in lines[4:4 + n]])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    picks = {r["ligand"]: r for r in csv.DictReader(open(SEL))}
    if len(picks) != 20:
        raise SystemExit(f"expected 20 ligands, found {len(picks)}")

    rows = []
    for lig in sorted(picks):
        rank = int(picks[lig]["pose_rank"])
        poses = sdf_mols(os.path.join(A33, f"sub29__{lig}", "hits.sdf"))
        if not 1 <= rank <= len(poses):
            raise SystemExit(f"{lig}: pose {rank} not among {len(poses)} in hits.sdf")
        heavy = [(e, xyz) for e, xyz in poses[rank - 1] if e.upper() != "H"]

        src = open(os.path.join(SRC, f"{lig}.pdb")).read().split("\n")
        lig_i = [i for i, l in enumerate(src)
                 if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "LIG"]
        if len(lig_i) != len(heavy):
            raise SystemExit(f"{lig}: {len(lig_i)} LIG atoms in the PDB vs "
                             f"{len(heavy)} heavy atoms in the pose")
        pdb_el = [(src[i][76:78].strip() or src[i][12:16].strip()[0]).upper()
                  for i in lig_i]
        if pdb_el != [e.upper() for e, _ in heavy]:
            raise SystemExit(f"{lig}: element sequence differs between the PDB LIG "
                             f"records and the ICM pose - refusing to graft")

        moved = []
        for i, (_, (x, y, z)) in zip(lig_i, heavy):
            old = (float(src[i][30:38]), float(src[i][38:46]), float(src[i][46:54]))
            moved.append(math.dist(old, (x, y, z)))
            src[i] = f"{src[i][:30]}{x:8.3f}{y:8.3f}{z:8.3f}{src[i][54:]}"
        open(os.path.join(OUT, f"{lig}.pdb"), "w").write("\n".join(src))

        rows.append({"ligand": lig, "picked_set": picks[lig]["picked_set"],
                     "pose_rank": rank,
                     "n_usable_poses": picks[lig]["n_usable_poses"],
                     "gscore": picks[lig]["gscore"],
                     "rmsd_to_sub29_ligand": picks[lig]["rmsd_to_sub29_ligand"],
                     "n_atoms_moved": len(moved),
                     "max_atom_displacement": round(max(moved), 3)})

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for r in rows:
            z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")

    import statistics
    d = sorted(float(r["rmsd_to_sub29_ligand"]) for r in rows)
    print(f"grafted {len(rows)}/20 ligands, element sequence verified for each")
    print(f"ligand RMSD to submission 29: median {statistics.median(d):.2f} A, "
          f"min {d[0]:.2f}, max {d[-1]:.2f}")
    print(f"wrote {os.path.relpath(ZIP, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
