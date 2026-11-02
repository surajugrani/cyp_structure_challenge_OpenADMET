"""Replace submission 40's poses with the minimised ones the selection was computed from.

01_select.py ranked the three ESMFold2 runs by Glide `mininplace` gscore and staged each
winner's co-folded structure. `mininplace` restrained-minimises the pose before scoring
it, so the score belongs to the minimised pose, not the co-folded one. This grafts the
minimised ligand coordinates out of analysis 28's pose-viewer files so the object
submitted is the object that was scored.

Submission 34 did the same graft on the `inplace` picks and moved LDDT-PLI by 0.0001, so
the minimisation itself is not expected to change the score - what is being tested here
is the selection, with the inconsistency removed rather than left in.

Each ligand's source PDB is the winning run's own submission, not a single base, because
the picks span all three runs.
"""
import csv, os, shutil, sys
from schrodinger import structure

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
HERE = os.path.join(ROOT, "submissions", "40_mininplace_best_of_three")
SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
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
    picks = {r["ligand"]: r["picked_set"] for r in csv.DictReader(
        open(os.path.join(HERE, "results", "selection.csv")))}
    rows, bad = [], []
    for lig in sorted(picks):
        sset = picks[lig]
        pv = os.path.join(A28, sset, lig, "mininplace", "score_pv.mae")
        src_pdb = os.path.join(ROOT, "submissions", SETS[sset], "submission",
                               f"{lig}.pdb")
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
