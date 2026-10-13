"""Cheap read: how far would each candidate selection rule move us from submission 06?

No zip is built. For every candidate rule, report (a) how many of the 20 structures would
change versus submission 06's picks - best PL-iPTM among PoseBusters-passing samples of
the no-template run, our best leaderboard result at 0.4291 - and (b) the pocket-superposed,
symmetry-corrected ligand RMSD between the candidate pose and 06's pose.

Rules compared:
  06_pliptm      baseline, best PL-iPTM          (no-template pool)
  best_gscore    best Glide gscore               (no-template pool)
  centroid_04    lowest mean RMSD to the rest    (no-template pool)
  centroid_pool  lowest mean RMSD to the rest    (no-template + 38FW-template pools)

A rule that changes many structures is a large bet; the RMSD column says whether the
changes are refinements or reorientations. Analysis 13 showed one 3.5 A swap cost 0.026
LDDT-PLI, so these numbers set the risk before a submission slot is spent.
"""
import csv, glob, json, os, sys
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S04 = os.path.join(ROOT, "submissions", "04_af3_default")
S08 = os.path.join(ROOT, "submissions", "08_af3_1template")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
RDLogger.DisableLog("rdApp.*")
POCKET_R = 10.0


def parse(p):
    ca, lig = {}, []
    for l in open(p):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if (l[76:78].strip() or l[12:16].strip()[0]) == "H":
            continue
        if l[17:20].strip() == "LIG":
            lig.append(l)
        elif l[12:16].strip() == "CA":
            ca[(l[21], l[22:27])] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return ca, lig


def xyz(ls):
    return np.array([(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in ls], float)


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def ligmol(ls, smi):
    m = Chem.MolFromPDBBlock("".join(ls) + "END\n", removeHs=True, sanitize=False)
    if m is None:
        return None
    try:
        return assign_bond_orders_from_template(Chem.MolFromSmiles(smi), m)
    except Exception:
        return m


def rmsd(a, b, smi):
    (ca_a, la_l), (ca_b, lb_l) = a, b
    la = xyz(la_l)
    keys = [k for k, v in ca_a.items()
            if np.min(np.linalg.norm(la - np.array(v), axis=1)) <= POCKET_R and k in ca_b]
    if len(keys) < 20:
        return None
    R, t = kabsch(np.array([ca_b[k] for k in keys], float),
                  np.array([ca_a[k] for k in keys], float))
    lb = (R @ xyz(lb_l).T).T + t
    ma, mb = ligmol(la_l, smi), ligmol(lb_l, smi)
    if ma is not None and mb is not None and ma.GetNumAtoms() == len(la) \
            and mb.GetNumAtoms() == len(lb):
        try:
            pr = Chem.Mol(mb); cf = pr.GetConformer()
            for i, p in enumerate(lb):
                cf.SetAtomPosition(i, p.tolist())
            return float(rdMolAlign.CalcRMS(pr, ma))
        except Exception:
            pass
    if len(la) != len(lb):
        return None
    return float(np.sqrt(((la - lb) ** 2).sum(1).mean()))


def pl_of(resdir, lig, tag):
    g = glob.glob(os.path.join(resdir, lig.lower(), tag, "*summary_confidences.json"))
    return json.load(open(g[0]))["chain_pair_iptm"][0][2] if g else None


def main():
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    gsc = {}
    for r in csv.DictReader(open(os.path.join(
            S04, "glide_rescoring", "results", "glide_scores.csv"))):
        if r["status"] == "ok" and r["gscore"] and abs(float(r["gscore"])) <= 100:
            gsc[(r["ligand"], r["sample"])] = float(r["gscore"])

    pool04, pool08 = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(os.path.join(
            S04, "glide_rescoring", "results", "pb_after_prepwizard.csv"))):
        if r["after_pass"] != "True":
            continue
        p = os.path.join(S04, "glide_rescoring", "work",
                         f"{r['ligand']}__{r['sample']}", f"{r['ligand']}__{r['sample']}_pw.pdb")
        if os.path.exists(p):
            pool04[r["ligand"]].append({"run": "04", "sample": r["sample"], "path": p,
                                        "pl": pl_of(os.path.join(S04, "results"),
                                                    r["ligand"], r["sample"]),
                                        "g": gsc.get((r["ligand"], r["sample"]))})
    for r in csv.DictReader(open(os.path.join(S08, "results", "sample_pb.csv"))):
        if r["pb_pass"].lower() in ("true", "1"):
            pool08[r["ligand"]].append({"run": "08", "sample": r["sample"],
                                        "path": r["path"], "pl": float(r["pl_iptm"]),
                                        "g": None})

    def centroid(cands, smi):
        if len(cands) < 2:
            return cands[0]
        M = np.full((len(cands), len(cands)), np.nan)
        for i in range(len(cands)):
            for j in range(len(cands)):
                if i != j:
                    v = rmsd(cands[i]["parsed"], cands[j]["parsed"], smi)
                    if v is not None:
                        M[i, j] = v
        mean = [np.nanmean(M[i]) if not np.all(np.isnan(M[i])) else np.inf
                for i in range(len(cands))]
        return cands[int(np.argmin(mean))]

    out = []
    for lig in sorted(smiles):
        p4 = pool04.get(lig, [])
        pall = p4 + pool08.get(lig, [])
        if not p4:
            continue
        for c in pall:
            c["parsed"] = parse(c["path"])
        base = max([c for c in p4 if c["pl"] is not None], key=lambda x: x["pl"])
        cands = {
            "best_gscore": min([c for c in p4 if c["g"] is not None],
                               key=lambda x: x["g"], default=None),
            "centroid_04": centroid(p4, smiles[lig]),
            "centroid_pool": centroid(pall, smiles[lig]),
        }
        row = {"ligand": lig, "n_04": len(p4), "n_pool": len(pall),
               "base_06": base["sample"]}
        for name, c in cands.items():
            if c is None:
                row[name] = ""; row[name + "_rmsd"] = ""
                continue
            same = c["run"] == "04" and c["sample"] == base["sample"]
            d = 0.0 if same else rmsd(base["parsed"], c["parsed"], smiles[lig])
            row[name] = f"{c['run']}/{c['sample']}"
            row[name + "_rmsd"] = "" if d is None else round(d, 2)
        out.append(row)
        for c in pall:
            c.pop("parsed", None)

    with open(os.path.join(HERE, "results", "candidate_rules.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

    print(f"{len(out)} ligands\n")
    print(f"{'rule':16} {'differs from 06':>16} {'median shift':>13} {'max shift':>10} {'>2A':>5}")
    for name in ("best_gscore", "centroid_04", "centroid_pool"):
        ds = [r[name + "_rmsd"] for r in out
              if r[name] and r[name + "_rmsd"] != "" and r[name + "_rmsd"] > 0]
        n = len(ds)
        print(f"{name:16} {n:>13}/{len(out)} "
              f"{(np.median(ds) if ds else 0):13.2f} {(max(ds) if ds else 0):10.2f} "
              f"{sum(d > 2 for d in ds):5d}")
    npool = sum(r["centroid_pool"].startswith("08") for r in out if r["centroid_pool"])
    print(f"\ncentroid_pool picks a template-run (08) pose for {npool}/{len(out)} ligands")
    print("per-ligand detail in results/candidate_rules.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
