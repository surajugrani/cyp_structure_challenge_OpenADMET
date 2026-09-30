"""Analysis 10, step 2: PoseBusters every AF3 sample, then ask which confidence metric
best predicts a PoseBusters pass.

The question this answers: AF3 emits 5 samples per ligand, and submission 04 picked the
one with the highest protein-ligand iPTM - which gave only 13/20 PoseBusters passes. If
some other reported metric tracks physical validity better, selecting on that would raise
the pass rate for free, with no extra compute.

Metrics tested, all from AF3's own summary_confidences.json:
  pl_iptm          chain_pair_iptm[0][2], protein vs ligand chain  (what submission 04 used)
  iptm             global, dominated by the well-predicted protein-haem interface
  ptm, ranking_score
  chain_iptm_lig   the ligand chain's own iptm
  has_clash        AF3's own clash flag - the most direct candidate
  fraction_disordered

For each metric the report gives: how often the metric-selected sample passes PoseBusters
(the number that matters), and the pass-rate difference between samples scoring high vs
low on it.
"""
import csv, glob, json, os, sys
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
AF3 = os.path.join(ROOT, "submissions", "04_af3_default", "results")
SAMP = os.path.join(HERE, "work", "samples")
RES = os.path.join(HERE, "results")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")

METRICS = ["pl_iptm", "iptm", "ptm", "ranking_score", "chain_iptm_lig",
           "has_clash", "fraction_disordered"]
# True  = a HIGHER value should mean a better (more likely PoseBusters-passing) pose
HIGHER_BETTER = {"pl_iptm": True, "iptm": True, "ptm": True, "ranking_score": True,
                 "chain_iptm_lig": True, "has_clash": False, "fraction_disordered": False}


def confidences(name, tag):
    f = glob.glob(os.path.join(AF3, name.lower(), tag, "*summary_confidences.json"))[0]
    j = json.load(open(f))
    cp = j.get("chain_pair_iptm") or [[None] * 3] * 3
    ci = j.get("chain_iptm") or [None, None, None]
    return {"pl_iptm": cp[0][2], "iptm": j.get("iptm"), "ptm": j.get("ptm"),
            "ranking_score": j.get("ranking_score"), "chain_iptm_lig": ci[2],
            "has_clash": j.get("has_clash"),
            "fraction_disordered": j.get("fraction_disordered")}


def pb_failures(path, smi, pb, tmp, name):
    pred = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    frags = Chem.SplitMolByPDBResidues(pred)
    lig = frags.pop("LIG", None)
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    if lig is None or prot is None:
        return ["parse_failure"]
    tpl = Chem.MolFromSmiles(smi)
    try:
        lig = assign_bond_orders_from_template(tpl, lig)
    except Exception:
        pass
    r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot)
    bc = [c for c in r.columns if r[c].dtype == bool]
    row = r[bc].iloc[0]
    return sorted(row.index[~row].tolist())


def main():
    os.makedirs(RES, exist_ok=True)
    tmp = os.path.join(RES, "_tmp"); os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = PoseBusters(config="dock")

    rows = []
    for i, f in enumerate(sorted(glob.glob(os.path.join(SAMP, "*.pdb"))), 1):
        base = os.path.basename(f)[:-4]
        name, tag = base.split("__")
        failed = pb_failures(f, smiles[name], pb, tmp, name)
        rec = {"ligand": name, "sample": tag, "pb_pass": not failed,
               "pb_failed": "; ".join(failed), **confidences(name, tag)}
        rows.append(rec)
        if i % 10 == 0:
            print(f"  {i}/100 scored", flush=True)

    with open(os.path.join(RES, "sample_level.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    npass = sum(r["pb_pass"] for r in rows)
    print(f"\n=== {npass}/{len(rows)} samples pass PoseBusters ===")
    byl = defaultdict(list)
    for r in rows:
        byl[r["ligand"]].append(r)
    solvable = sum(any(x["pb_pass"] for x in v) for v in byl.values())
    print(f"ligands with at least one passing sample: {solvable}/20  <- selection ceiling")

    print(f"\n{'metric':22}{'picked sample passes':>22}{'pass rate hi vs lo':>22}")
    out = []
    for m in METRICS:
        hi = HIGHER_BETTER[m]
        picked = 0
        for lig, v in byl.items():
            vals = [x for x in v if x[m] is not None]
            if not vals:
                continue
            best = (max if hi else min)(vals, key=lambda x: x[m])
            picked += best["pb_pass"]
        vv = [(r[m], r["pb_pass"]) for r in rows if r[m] is not None]
        if vv:
            med = float(np.median([a for a, _ in vv]))
            above = [p for a, p in vv if a > med]
            below = [p for a, p in vv if a <= med]
            ra = np.mean(above) if above else float("nan")
            rb = np.mean(below) if below else float("nan")
            good, bad = (ra, rb) if hi else (rb, ra)
            spread = f"{good:.2f} vs {bad:.2f}"
        else:
            spread = "n/a"
        out.append({"metric": m, "picked_passes": picked, "spread": spread})
        print(f"{m:22}{f'{picked}/20':>22}{spread:>22}")

    with open(os.path.join(RES, "metric_comparison.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
    print(f"\nWrote {RES}/sample_level.csv and metric_comparison.csv")


if __name__ == "__main__":
    main()
