"""Same pose-divergence measurement as 01, but on the NO-TEMPLATE run (submission 04).

01 showed that in the 38FW-template run the 5 AF3 samples of a ligand differ by a median
1.81 A, mean 2.01 A, max 6.42 A - genuinely different binding modes, not variations. This
repeats it on the no-template predictions so the two runs can be compared, and so the
pool spread of the run that produced the best submissions (05/06) is known.
"""
import csv, os, sys
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S04 = os.path.join(ROOT, "submissions", "04_af3_default")
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
    mol = Chem.MolFromPDBBlock("".join(lines) + "END\n", removeHs=True, sanitize=False)
    if mol is None:
        return None
    try:
        return assign_bond_orders_from_template(Chem.MolFromSmiles(smi), mol)
    except Exception:
        return mol


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
        return None
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
            return float(rdMolAlign.CalcRMS(probe, ma))
        except Exception:
            pass
    if len(la) != len(lb):
        return None
    return float(np.sqrt(((la - lb) ** 2).sum(1).mean()))


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    by = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(
            S04, "glide_rescoring", "results", "pb_after_prepwizard.csv"))):
        if r["after_pass"] != "True":
            continue
        lig, tag = r["ligand"], r["sample"]
        p = os.path.join(S04, "glide_rescoring", "work",
                         f"{lig}__{tag}", f"{lig}__{tag}_pw.pdb")
        if os.path.exists(p):
            by[lig].append({"sample": tag, "path": p})

    out, pairs = [], []
    for lig, pool in sorted(by.items()):
        for r in pool:
            r["parsed"] = parse(r["path"])
        vals = []
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                v = rmsd(pool[i]["parsed"], pool[j]["parsed"], smiles[lig])
                if v is not None:
                    vals.append(v)
                    pairs.append({"ligand": lig, "a": pool[i]["sample"],
                                  "b": pool[j]["sample"], "rmsd": round(v, 3)})
        out.append({"ligand": lig, "n_pool": len(pool),
                    "pool_mean_rmsd": round(float(np.mean(vals)), 3) if vals else None,
                    "pool_median_rmsd": round(float(np.median(vals)), 3) if vals else None,
                    "pool_min_rmsd": round(float(np.min(vals)), 3) if vals else None,
                    "pool_max_rmsd": round(float(np.max(vals)), 3) if vals else None})
        print(f"{lig}  n={len(pool)}  mean {out[-1]['pool_mean_rmsd']}  "
              f"min {out[-1]['pool_min_rmsd']}  max {out[-1]['pool_max_rmsd']}", flush=True)
        for r in pool:
            r.pop("parsed", None)

    R = os.path.join(HERE, "results")
    for fn, rr in (("notemplate_pool_summary.csv", out),
                   ("notemplate_pairwise_rmsd.csv", pairs)):
        with open(os.path.join(R, fn), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rr[0])); w.writeheader(); w.writerows(rr)
    allp = [r["rmsd"] for r in pairs]
    print(f"\nall {len(allp)} within-ligand pairs: median {np.median(allp):.2f} A, "
          f"mean {np.mean(allp):.2f} A, max {np.max(allp):.2f} A")
    return 0


if __name__ == "__main__":
    sys.exit(main())
