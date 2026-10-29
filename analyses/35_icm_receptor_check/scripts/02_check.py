"""Were the ICM poses docked into a receptor that differs from the one we ship?

Analysis 34 scored analysis 33's ICM poses with Glide `inplace` and 177 of 200 came
back at the 10000 rejection sentinel, with van der Waals energies around 1e5 kcal/mol.
Measured against submission 29's receptor, the poses sit 1.6-2.9 A from the nearest
protein heavy atom, where an OPLS carbon-carbon contact wants ~3.4 A.

The hypothesis under test: ICM's `convertObject` does not dock into the PDB as handed
to it. It adds hydrogens, applies its own typing, and can move side chains. If the
pocket opened up during that conversion, the poses are not clashing at all in the
receptor they were produced in - they only clash in the receptor we would ship them
back into, which is a different structure.

This distinguishes the two cases:
  * contacts are short against BOTH receptors -> the poses really are strained, and
    ICM's force field simply tolerates contacts OPLS does not
  * contacts are reasonable against the ICM receptor and short against submission 29's
    -> the receptor moved, and the whole redock-into-sub29 route has a problem that no
    scoring mode fixes, because the pose and the receptor we would package it with
    never coexisted

No quality claim either way. This is about whether a pose is physically compatible with
the receptor it would be shipped in, which is a separate question from whether it is
correct.
"""
import csv, glob, math, os, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A33 = os.path.join(ROOT, "analyses", "33_icm_redock_sub29", "structures")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "submission")
ICM_HOME = "/opt/applications/icm/3.9.4"
SCRIPT = os.path.join(HERE, "scripts", "01_export.icm")
OUT = os.path.join(HERE, "results", "receptor_contacts.csv")


def heavy(path, skip_lig=True):
    out = []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        if skip_lig and l[17:20].strip() == "LIG":
            continue
        out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def poses(path):
    out = []
    for blk in open(path).read().split("$$$$"):
        ls = blk.strip("\n").split("\n")
        if len(ls) < 4:
            continue
        try:
            na = int(ls[3][:3])
        except ValueError:
            continue
        at = [(float(l[0:10]), float(l[10:20]), float(l[20:30]))
              for l in ls[4:4 + na] if l[31:34].strip() != "H"]
        if at:
            out.append(at)
    return out


def nearest(pose, rec):
    return min(math.dist(p, r) for p in pose for r in rec)


def main():
    if not (os.environ.get("SLURM_JOB_ID") or
            os.environ.get("ALLOW_LOGIN_NODE") == "1"):
        raise SystemExit("REFUSING to run ICM outside SLURM. Use sbatch.")
    rows = []
    for d in sorted(glob.glob(os.path.join(A33, "sub29__*"))):
        lig = os.path.basename(d).split("__")[1]
        ob = os.path.join(d, "proj_rec.ob")
        icm_pdb = os.path.join(d, "icm_receptor.pdb")
        if os.path.exists(ob) and not os.path.exists(icm_pdb):
            subprocess.run([os.path.join(ICM_HOME, "icm64"), "-s", SCRIPT,
                            f"s_in={ob}", f"s_out={icm_pdb}"],
                           cwd=d, capture_output=True, text=True, timeout=600)
        if not os.path.exists(icm_pdb):
            print(f"{lig}: could not export ICM receptor", flush=True)
            continue
        ps = poses(os.path.join(d, "hits.sdf"))
        shipped = heavy(os.path.join(SUB29, f"{lig}.pdb"))
        icmrec = heavy(icm_pdb)
        for i, p in enumerate(ps, 1):
            rows.append({"ligand": lig, "pose_rank": i,
                         "min_dist_icm_receptor": round(nearest(p, icmrec), 2),
                         "min_dist_sub29_receptor": round(nearest(p, shipped), 2)})
        if ps:
            a = min(r["min_dist_icm_receptor"] for r in rows if r["ligand"] == lig)
            b = min(r["min_dist_sub29_receptor"] for r in rows if r["ligand"] == lig)
            print(f"{lig:16} closest contact  ICM receptor {a:5.2f} A   "
                  f"submission 29 receptor {b:5.2f} A", flush=True)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, ["ligand", "pose_rank", "min_dist_icm_receptor",
                                "min_dist_sub29_receptor"])
        w.writeheader(); w.writerows(rows)

    if rows:
        ai = sorted(r["min_dist_icm_receptor"] for r in rows)
        bi = sorted(r["min_dist_sub29_receptor"] for r in rows)
        print(f"\n{len(rows)} poses. median closest contact: "
              f"ICM receptor {ai[len(ai)//2]:.2f} A, "
              f"submission 29 receptor {bi[len(bi)//2]:.2f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
