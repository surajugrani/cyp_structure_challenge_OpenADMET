"""How far apart are the poses that submissions 08 and 09 chose?

09 is 08 with the selection key changed from mean(PL,HL) iPTM to PL-iPTM alone. Exactly 3
of 20 structures differ, yet LDDT-PLI fell 0.0258 and BiSyRMSD rose 0.17 A. Since only
~10 complexes are scored, that implies ~1.7 A of summed pose movement concentrated in
1-3 structures - a pose flip rather than a refinement.

This measures it directly: for each differing ligand, superpose the two prepped samples on
the CA atoms of the pocket (protein residues within 10 A of the ligand, mirroring how
BiSyRMSD is scored) and compute the symmetry-corrected ligand RMSD.

Also reported for context: the full pairwise RMSD matrix over ALL PoseBusters-passing
samples of those ligands, so the swap can be judged against the pool's natural spread.
"""
import csv, os, sys
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S08 = os.path.join(ROOT, "submissions", "08_af3_1template")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
RDLogger.DisableLog("rdApp.*")
POCKET_R = 10.0


def parse(path):
    ca, lig = {}, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if (l[76:78].strip() or l[12:16].strip()[0]) == "H":
            continue
        if l[17:20].strip() == "LIG":
            lig.append(l)
        elif l[12:16].strip() == "CA":
            ca[(l[21], l[22:27])] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return ca, lig


def xyz(lines):
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54]))
                     for l in lines], float)


def ligmol(lines, smi):
    m = Chem.MolFromPDBBlock("".join(lines) + "END\n", removeHs=True, sanitize=False)
    if m is None:
        return None
    try:
        return assign_bond_orders_from_template(Chem.MolFromSmiles(smi), m)
    except Exception:
        return m


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def rmsd(a, b, smi):
    (ca_a, lig_a), (ca_b, lig_b) = a, b
    la = xyz(lig_a)
    keys = [k for k, v in ca_a.items()
            if np.min(np.linalg.norm(la - np.array(v), axis=1)) <= POCKET_R
            and k in ca_b]
    if len(keys) < 20:
        return None, 0
    R, t = kabsch(np.array([ca_b[k] for k in keys], float),
                  np.array([ca_a[k] for k in keys], float))
    lb = (R @ xyz(lig_b).T).T + t
    ma, mb = ligmol(lig_a, smi), ligmol(lig_b, smi)
    if ma is not None and mb is not None and ma.GetNumAtoms() == len(la) \
            and mb.GetNumAtoms() == len(lb):
        try:
            probe = Chem.Mol(mb); cf = probe.GetConformer()
            for i, p in enumerate(lb):
                cf.SetAtomPosition(i, p.tolist())
            return float(rdMolAlign.CalcRMS(probe, ma)), len(keys)
        except Exception:
            pass
    if len(la) != len(lb):
        return None, len(keys)
    return float(np.sqrt(((la - lb) ** 2).sum(1).mean())), len(keys)


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    rows = [r for r in csv.DictReader(open(os.path.join(S08, "results", "sample_pb.csv")))
            if r["pb_pass"].lower() in ("true", "1")]
    by = defaultdict(list)
    for r in rows:
        by[r["ligand"]].append(r)

    out, pool_rows = [], []
    for lig, pool in sorted(by.items()):
        for r in pool:
            r["parsed"] = parse(r["path"])
        p08 = max(pool, key=lambda x: float(x["mean_pl_hl"]))
        p09 = max(pool, key=lambda x: float(x["pl_iptm"]))
        d, nk = (0.0, 0) if p08["sample"] == p09["sample"] else \
            rmsd(p08["parsed"], p09["parsed"], smiles[lig])
        vals = []
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                v, _ = rmsd(pool[i]["parsed"], pool[j]["parsed"], smiles[lig])
                if v is not None:
                    vals.append(v)
                    pool_rows.append({"ligand": lig, "a": pool[i]["sample"],
                                      "b": pool[j]["sample"], "rmsd": round(v, 3)})
        out.append({"ligand": lig, "differs": p08["sample"] != p09["sample"],
                    "s08_pick": p08["sample"], "s09_pick": p09["sample"],
                    "swap_rmsd": None if d is None else round(d, 3),
                    "n_pocket_ca": nk, "n_pool": len(pool),
                    "pool_mean_rmsd": round(float(np.mean(vals)), 3) if vals else None,
                    "pool_max_rmsd": round(float(np.max(vals)), 3) if vals else None})
        tag = "  <<< DIFFERS" if out[-1]["differs"] else ""
        print(f"{lig}  swap {str(out[-1]['swap_rmsd']):>6} A   "
              f"pool mean {str(out[-1]['pool_mean_rmsd']):>6} / "
              f"max {str(out[-1]['pool_max_rmsd']):>6} A  (n={len(pool)}){tag}",
              flush=True)
        for r in pool:
            r.pop("parsed", None)

    R = os.path.join(HERE, "results")
    for fn, rr in (("swap_rmsd.csv", out), ("pool_pairwise_rmsd.csv", pool_rows)):
        with open(os.path.join(R, fn), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rr[0])); w.writeheader(); w.writerows(rr)
    diff = [r for r in out if r["differs"]]
    print(f"\n{len(diff)} ligands differ between 08 and 09")
    for r in diff:
        print(f"  {r['ligand']}: {r['swap_rmsd']} A  "
              f"(pool mean {r['pool_mean_rmsd']}, max {r['pool_max_rmsd']})")
    allp = [r["rmsd"] for r in pool_rows]
    print(f"\nall {len(allp)} within-ligand pairs: median {np.median(allp):.2f} A, "
          f"mean {np.mean(allp):.2f} A, max {np.max(allp):.2f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
