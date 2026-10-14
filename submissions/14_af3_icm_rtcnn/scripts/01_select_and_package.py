"""Submission 14: AF3 no-template samples selected by ICM's RTCNN score.

Test #3 of the agreed three, and the last selector with measured skill behind it.

WHY RTCNN. Analysis 09 ranked every available pose-scoring term by how well it tracks
RMSD *within* one ligand's pool of poses - the quantity a selector actually needs:

    ICM RTCNN    +0.334        <- by far the best
    Glide lipo   +0.115
    Glide gscore +0.086
    Glide emodel +0.057

Submission 13 tested the Glide end of that table and lost badly: 0.3730 against
submission 06's 0.4291, a drop of 0.056. If the ranking in analysis 09 is real and
transfers from crystal receptors to co-folded ones, RTCNN should do much better than
that - possibly better than PL-iPTM itself.

WHAT CHANGES. 16 of 20 structures differ from submission 06. The Fe chemistry barely
moves: of the 15 ligands carrying a nitrogen donor, RTCNN selection coordinates 11
against 06's 12, median Fe-ligand distance 2.34 A against 2.25 A. So, as with Glide, the
swaps are largely lateral - different sample, same binding mode.

Notably RTCNN and analysis 09's preferred RTCNN+Hbond model pick identically on all 20
ligands here, so the Hbond term contributes nothing on this data and the simpler rule is
used.

Scores come from analysis 20: each of the 100 PrepWizard-minimised AF3 samples scored in
place with `_dockScan -n` on its own receptor conformation. 94 samples are both scored
and PoseBusters-passing; 1 failed scoring outright.

CAVEAT on `-n`: it minimises before scoring, so the scored pose sits ~0.3-1.4 A from the
input (median 0.86 A, measured on the crystal set in analysis 08). The shipped coordinates
are the ORIGINAL PrepWizard pose, not ICM's minimised one - ICM is used only to rank.
"""
import csv, glob, json, os, sys, zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S04 = os.path.join(ROOT, "submissions", "04_af3_default")
AF3 = os.path.join(S04, "results")
WORK = os.path.join(S04, "glide_rescoring", "work")
PBCSV = os.path.join(S04, "glide_rescoring", "results", "pb_after_prepwizard.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_af3_icm_rtcnn.zip")
A20 = os.path.join(ROOT, "analyses", "20_icm_score_af3", "results")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")


def icm_scores():
    """RTCNN per (ligand, sample) from analysis 20; lower is better."""
    out = {}
    for r in csv.DictReader(open(os.path.join(A20, "icm_sample_scores.csv"))):
        if r["status"] != "ok" or not r.get("RTCNN"):
            continue
        try:
            out[(r["ligand"], r["sample"])] = float(r["RTCNN"])
        except ValueError:
            pass
    return out


def conf(name, tag):
    f = glob.glob(os.path.join(AF3, name.lower(), tag, "*summary_confidences.json"))[0]
    cp = json.load(open(f))["chain_pair_iptm"]
    return cp[0][2], cp[1][2]


def struct_errors(lines, name, smi, tmp):
    e = []
    at = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if len({l[21] for l in at}) > 3:
        e.append("more than 3 chains")
    lig = [l for l in at if l[17:20].strip() == "LIG"]
    if not lig:
        e.append("no LIG residue")
    elif len({(l[21], l[22:27]) for l in lig}) > 1:
        e.append("more than one LIG residue")
    nhem = len([l for l in at if l[17:20].strip() == "HEM"
                and (l[76:78].strip() or l[12:16].strip()[0]) != "H"])
    if nhem != 43:
        e.append(f"HEM has {nhem} heavy atoms")
    seen = {}
    for l in at:
        k = (l[21], l[22:27])
        if k in seen and seen[k] != l[17:20].strip():
            e.append("duplicate residue number")
        seen[k] = l[17:20].strip()
    per = defaultdict(list)
    for l in at:
        per[(l[21], l[22:27])].append(l[12:16].strip())
    if any(len(v) != len(set(v)) for v in per.values()):
        e.append("duplicate atom names in a residue")
    s = [l[6:11] for l in at]
    if len(s) != len(set(s)):
        e.append("duplicate atom serials")
    p = os.path.join(tmp, f"{name}_lig.pdb")
    open(p, "w").writelines(lig + ["END\n"])
    ref = Chem.MolFromSmiles(smi)
    probe = Chem.MolFromPDBFile(p, removeHs=True, sanitize=False)
    if ref is None or probe is None:
        e.append("ligand parse failure")
    else:
        try:
            AllChem.AssignBondOrdersFromTemplate(ref, probe)
        except ValueError:
            e.append("ligand connectivity != SMILES")
    return e


def pb_fail(path, smi, pb):
    pred = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    if pred is None:
        return ["unparseable"]
    frags = Chem.SplitMolByPDBResidues(pred)
    lig = frags.pop("LIG", None)
    prot = None
    for fr in frags.values():
        prot = fr if prot is None else Chem.CombineMols(prot, fr)
    if lig is None or prot is None:
        return ["no LIG or protein"]
    tpl = Chem.MolFromSmiles(smi)
    try:
        lig = assign_bond_orders_from_template(tpl, lig)
    except Exception:
        pass
    r = pb.bust(mol_pred=lig, mol_true=tpl, mol_cond=prot)
    bc = [c for c in r.columns if r[c].dtype == bool]
    return sorted(r[bc].iloc[0].index[~r[bc].iloc[0]].tolist())


def main():
    for d in (OUT, RES):
        os.makedirs(d, exist_ok=True)
    tmp = os.path.join(RES, "_tmp"); os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    gsc = icm_scores()
    passing = defaultdict(list)
    for r in csv.DictReader(open(PBCSV)):
        if r["after_pass"] != "True":
            continue
        if (r["ligand"], r["sample"]) not in gsc:
            continue
        p = os.path.join(WORK, f"{r['ligand']}__{r['sample']}",
                         f"{r['ligand']}__{r['sample']}_pw.pdb")
        if os.path.exists(p):
            pl, hl = conf(r["ligand"], r["sample"])
            passing[r["ligand"]].append({"sample": r["sample"], "path": p,
                                         "pl_iptm": pl, "hl_iptm": hl,
                                         "mean_pl_hl": (pl + hl) / 2.0,
                                         "rtcnn": gsc[(r["ligand"], r["sample"])]})
    print(f"{sum(len(v) for v in passing.values())} PB-passing prepped samples, "
          f"{len(passing)}/20 ligands covered")

    sel = []
    for n in sorted(smiles):
        pool = passing.get(n, [])
        if not pool:
            raise RuntimeError(f"{n}: no PoseBusters-passing prepped sample")
        best = min(pool, key=lambda x: x["rtcnn"])   # most negative wins
        alt_pl = max(pool, key=lambda x: x["pl_iptm"])["sample"]
        alt_mean = max(pool, key=lambda x: x["mean_pl_hl"])["sample"]
        with open(best["path"]) as f, open(os.path.join(OUT, f"{n}.pdb"), "w") as g:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    if (l[76:78].strip() or l[12:16].strip()[0]) == "H":
                        continue
                g.write(l)
        sel.append({"ligand": n, "n_passing": len(pool), "chosen": best["sample"],
                    "rtcnn": best["rtcnn"],
                    "hl_iptm": best["hl_iptm"], "pl_iptm": best["pl_iptm"],
                    "mean_pl_hl": round(best["mean_pl_hl"], 4),
                    "pl_only_would_pick": alt_pl, "mean_would_pick": alt_mean,
                    "differs_from_06_pl": best["sample"] != alt_pl,
                    "differs_from_05_mean": best["sample"] != alt_mean})
        marks = []
        if sel[-1]["differs_from_06_pl"]:
            marks.append("!=06")
        if sel[-1]["differs_from_05_mean"]:
            marks.append("!=05")
        print(f"  {n}: {best['sample']} ({len(pool)} passing, "
              f"RTCNN {best['rtcnn']:.2f}, PL {best['pl_iptm']:.2f})  "
              f"{' '.join(marks)}")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sel[0])); w.writeheader(); w.writerows(sel)
    print(f"\ndiffers from 06 (PL-iPTM) in "
          f"{sum(r['differs_from_06_pl'] for r in sel)}/20; "
          f"from 05 (mean) in {sum(r['differs_from_05_mean'] for r in sel)}/20")

    pb = PoseBusters(config="dock")
    ns = npb = 0
    for n in sorted(smiles):
        p = os.path.join(OUT, f"{n}.pdb")
        se = struct_errors(open(p).readlines(), n, smiles[n], tmp)
        pf = pb_fail(p, smiles[n], pb)
        ns += not se
        npb += not pf
        if se or pf:
            print(f"  FINAL {n}: struct {se} | PB {pf}")
    print(f"\nFINAL structural {ns}/20   PoseBusters {npb}/20")
    if ns == 20 and npb == 20:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for n in sorted(smiles):
                z.write(os.path.join(OUT, f"{n}.pdb"), f"{n}.pdb")
        print(f"wrote {ZIP}")
    else:
        print("zip not written - failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
