"""Compare constrained vs unconstrained redocking on the 10 test structures."""
import csv, glob, os, sys
import numpy as np
sys.path.insert(0, "/gpfs/home/sugrani/SelfLearning/CYP_structure_OpenADMET/.pylibs")
from spyrmsd import rmsd as srm, io as sio
from schrodinger import structure

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A05 = os.path.join(ROOT, "analyses", "05_full_pdb_redocking")
WORK = os.path.join(HERE, "work")
COORD = ["8SPD", "4I4G", "6MA7", "6DAJ", "6UNI"]
NON = ["6UNJ", "3UA1", "7LXL", "6OO9", "4D6Z"]


def rmsds(pv, ref_sdf, tmp):
    """Symmetry-corrected RMSD of each pose to the crystal ligand."""
    ref = sio.loadmol(ref_sdf); ref.strip()
    out = []
    for i, st in enumerate(list(structure.StructureReader(pv))[1:]):
        f = os.path.join(tmp, f"p{i}.sdf")
        with structure.StructureWriter(f) as w:
            w.append(st)
        try:
            m = sio.loadmol(f); m.strip()
            out.append(float(srm.symmrmsd(ref.coordinates, m.coordinates,
                                          ref.atomicnums, m.atomicnums,
                                          ref.adjacency_matrix, m.adjacency_matrix)))
        except Exception:
            pass
    return out


def main():
    tmp = os.path.join(HERE, "results", "_tmp"); os.makedirs(tmp, exist_ok=True)
    base = {r["pdb_id"]: r for r in csv.DictReader(
        open(os.path.join(A05, "results", "redock_summary.csv")))}
    rows = []
    for pdb in COORD + NON:
        arm = "coordinating" if pdb in COORD else "non_coordinating"
        pv = os.path.join(WORK, pdb, f"{pdb}_mdock_pv.mae")
        ref = os.path.join(A05, "structures", pdb, f"{pdb}_crystal_lig.sdf")
        rec = {"pdb_id": pdb, "arm": arm,
               "uncons_top1": float(base[pdb]["rmsd_top1"]),
               "uncons_best": float(base[pdb]["rmsd_best"])}
        if os.path.exists(pv) and os.path.exists(ref):
            v = rmsds(pv, ref, tmp)
            if v:
                rec.update(cons_top1=round(v[0], 3), cons_best=round(min(v), 3),
                           cons_n=len(v))
        rows.append(rec)

    with open(os.path.join(HERE, "results", "rmsd_comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}),
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    print(f"{'pdb':7}{'arm':18}{'uncons top1':>13}{'cons top1':>11}"
          f"{'uncons best':>13}{'cons best':>11}   effect")
    for r in rows:
        ct = r.get("cons_top1"); cb = r.get("cons_best")
        if ct is None:
            print(f"{r['pdb_id']:7}{r['arm']:18}{r['uncons_top1']:>13.2f}"
                  f"{'no poses':>11}{r['uncons_best']:>13.2f}{'-':>11}   constraint rejected")
            continue
        d = r["uncons_top1"] - ct
        tag = "better" if d > 0.5 else ("worse" if d < -0.5 else "~same")
        print(f"{r['pdb_id']:7}{r['arm']:18}{r['uncons_top1']:>13.2f}{ct:>11.2f}"
              f"{r['uncons_best']:>13.2f}{cb:>11.2f}   {tag} ({d:+.2f})")

    for arm in ("coordinating", "non_coordinating"):
        v = [r for r in rows if r["arm"] == arm and r.get("cons_top1") is not None]
        if not v:
            continue
        u = np.array([r["uncons_top1"] for r in v]); c = np.array([r["cons_top1"] for r in v])
        print(f"\n{arm} (n={len(v)} with poses): median top-1 "
              f"{np.median(u):.2f} -> {np.median(c):.2f}   "
              f"under 2 A: {(u<2).sum()} -> {(c<2).sum()}")


if __name__ == "__main__":
    main()
