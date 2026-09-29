"""
Analysis 09, step 2a: can ICM's neural-network RTCNN score pick better poses than ICM's
own total Score?

No docking or rescoring is run. RTCNN (and mfScore) were already written into analysis
06's `<pdb>_hits.sdf` by `_dockScan`; the analysis 06 parser silently dropped RTCNN because
it filtered on a tag named `RTCNNscore` while the actual SD tag is `RTCNN`. This script
re-reads those files and re-ranks the poses that already exist.

Set: every structure in the full redocking set with RotB <= 8 that ICM docked (n = 33).
Poses are matched back to analysis 06's `icm_poses.csv` (which carries the RMSD to the
crystal pose) on the Score value, since the SDF pose order is not the CSV `pose` index.

Usage (light - pandas/RDKit only, no licence, no docking):
    $SCHRODINGER/run python3 scripts/02_rtcnn_rescore.py

Outputs:
    results/rtcnn_rescore_poses.csv     every pose with Score, RTCNN, mfScore, rmsd
    results/rtcnn_rescore_summary.csv   per-structure top-1 pick of each score
"""

import json
import os
import re

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANALYSES = os.path.dirname(HERE)
A05 = os.path.join(ANALYSES, "05_full_pdb_redocking")
A06 = os.path.join(ANALYSES, "06_icm_redocking")
RESULTS = os.path.join(HERE, "results")

MAX_ROTB = 8
CUTOFFS = (2.0, 2.5)
SCORES = ["Score", "RTCNN", "mfScore"]  # all "lower is better"


def sdf_tags(path):
    """Minimal SD-tag reader - avoids a Schrodinger StructureReader round-trip."""
    out = []
    for block in open(path).read().split("$$$$"):
        if "<Score>" not in block:
            continue
        out.append({m.group(1): m.group(2).strip()
                    for m in re.finditer(r"> *<([^>]+)>\s*\n([^\n]*)", block)})
    return out


def main():
    os.makedirs(RESULTS, exist_ok=True)

    smi = json.load(open(os.path.join(A05, "data", "ligand_smiles.json")))
    rotb = {}
    for k, v in smi.items():
        m = Chem.MolFromSmiles(v["smiles_stereo"])
        if m:
            rotb[k] = Descriptors.NumRotatableBonds(m)

    summ = pd.read_csv(os.path.join(A06, "results", "icm_summary.csv"))
    summ["RotB"] = summ.comp_id.map(rotb)
    keep = summ[(summ.RotB <= MAX_ROTB) & (summ.status == "ok")]

    poses = pd.read_csv(os.path.join(A06, "results", "icm_poses.csv"))

    allp, rows = [], []
    for r in keep.itertuples():
        sdf = os.path.join(A06, "structures", r.pdb_id, f"{r.pdb_id}_hits.sdf")
        if not os.path.exists(sdf):
            # A few structure dirs were rebuilt after analysis 06 and lost their pose
            # SDFs. RTCNN cannot be recovered for these without re-docking.
            print(f"  {r.pdb_id}: no hits.sdf (pose file lost), skipped")
            continue
        tags = pd.DataFrame([{"Score": float(t["Score"]), "RTCNN": float(t["RTCNN"]),
                              "mf": float(t["mfScore"])} for t in sdf_tags(sdf)])
        p = poses[poses.pdb_id == r.pdb_id].copy()
        # Join on Score; the SDF order is not the CSV pose index. The key must be the
        # "%.4f" STRING - icm_poses.csv was written with float_format="%.4f" (round-half-
        # away) while Series.round(4) uses banker's rounding, so e.g. -9.12645 becomes
        # -9.1265 in the CSV but -9.1264 from round(). That silently lost one pose per
        # affected structure and dropped the whole structure from the analysis.
        p["_k"] = p.Score.map(lambda v: f"{v:.4f}")
        tags["_k"] = tags.Score.map(lambda v: f"{v:.4f}")
        m = p.merge(tags[["_k", "RTCNN", "mf"]].drop_duplicates("_k"), on="_k", how="left")
        m = m.rename(columns={"mf": "mfScore_sdf"})
        if m.RTCNN.isna().any():
            print(f"  {r.pdb_id}: {int(m.RTCNN.isna().sum())} poses unmatched, skipped")
            continue
        m["RotB"] = r.RotB
        allp.append(m)

        rec = {"pdb_id": r.pdb_id, "comp_id": r.comp_id, "RotB": r.RotB,
               "n_poses": len(m), "oracle": m.rmsd.min()}
        for s, col in (("Score", "Score"), ("RTCNN", "RTCNN"), ("mfScore", "mfScore_sdf")):
            rec[f"{s}_top1"] = m.loc[m[col].idxmin(), "rmsd"]
            rec[f"{s}_rank_of_best"] = int(m[col].rank().iloc[m.rmsd.values.argmin()])
        rows.append(rec)

    pd.concat(allp).to_csv(os.path.join(RESULTS, "rtcnn_rescore_poses.csv"),
                           index=False, float_format="%.4f")
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(RESULTS, "rtcnn_rescore_summary.csv"),
             index=False, float_format="%.4f")

    print(d[["pdb_id", "comp_id", "RotB", "n_poses", "Score_top1", "RTCNN_top1",
             "mfScore_top1", "oracle"]].to_string(index=False, float_format="%.2f"))

    n = len(d)
    print(f"\n=== top-1 success, n = {n} ===")
    for cut in CUTOFFS:
        print(f"  < {cut} A")
        for s in SCORES:
            k = (d[f"{s}_top1"] < cut).sum()
            print(f"    {s:8} {k:2d}/{n} ({100*k/n:4.1f}%)")
        k = (d.oracle < cut).sum()
        print(f"    {'oracle':8} {k:2d}/{n} ({100*k/n:4.1f}%)")

    # who does RTCNN save, and who does it break?
    for cut in (2.0,):
        sw = d[(d.Score_top1 >= cut) & (d.RTCNN_top1 < cut)]
        br = d[(d.Score_top1 < cut) & (d.RTCNN_top1 >= cut)]
        print(f"\nRTCNN vs Score at {cut} A:  saves {len(sw)}  breaks {len(br)}")
        for lab, sub in (("SAVED", sw), ("BROKEN", br)):
            for r in sub.itertuples():
                print(f"  {lab:6} {r.pdb_id} ({r.comp_id}) RotB={r.RotB}  "
                      f"Score {r.Score_top1:.2f} -> RTCNN {r.RTCNN_top1:.2f} "
                      f"(oracle {r.oracle:.2f})")


if __name__ == "__main__":
    main()
