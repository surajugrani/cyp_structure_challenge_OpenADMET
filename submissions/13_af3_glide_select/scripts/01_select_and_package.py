"""Submission 13: AF3 no-template samples selected by best Glide docking score.

The first attempt to select with a physics score rather than a confidence metric. Same
100 PrepWizard-minimised AF3 samples as submissions 05, 06 and 10; only the key changes:
the most negative Glide `gscore` among a ligand's PoseBusters-passing samples.

WHAT CHANGES. 16 of 20 structures differ from submission 06 (our best, 0.4291) - the
largest perturbation of any attempt so far. Analysis 15 measured the pose movement:
median 2.23 A, max 6.65 A, with 8 structures moving more than 2 A.

WHAT DOES NOT CHANGE. Analysis 16 measured the Fe chemistry of both rules and they are
near-identical. Of the 12 ligands carrying a free aromatic nitrogen and therefore able to
coordinate the haem iron, both PL-iPTM and Glide selection land on a coordinating pose
for 11. Median Fe-ligand distance 2.25 A (06) against 2.32 A (Glide). So the 16 swaps are
largely lateral - a different sample with the same binding mode - and the real movement is
concentrated in ligands that cannot coordinate at all: OCNT-2395328 (4.35 -> 2.86 A) and
OCNT-2395342 (4.50 -> 3.48 A), where Glide pulls the ligand toward an iron it has no
chemistry to contact.

THE EVIDENCE IS AGAINST THIS WORKING:
  * analysis 10 - Glide scoring selected PoseBusters-valid poses BELOW random (10/20 vs
    11.4/20 expected) and mildly preferred the poses PoseBusters rejects.
  * analysis 14 - Glide scores do not track potency either; every cross-validated model
    on (evdw, ecoul) has a negative Q2, and gscore's correlation has the wrong sign.
  * analysis 09 - Glide's within-ligand pose-ranking power is near zero (rho +0.086 for
    gscore, +0.115 for lipo), far below ICM RTCNN's +0.334.

So this is high variance with no measured skill behind it. Submitted to close the question
rather than because it is expected to win. Submission 09 is the cautionary precedent: one
ligand moving 3.5 A cost 0.0258 LDDT-PLI, and this run moves 8 ligands by more than 2 A.

Nothing is recomputed - AF3 inference, PrepWizard and the Glide rescoring all exist from
submission 04. Glide's 10000.0 "no score" sentinel is filtered by magnitude, since it can
appear with status "ok".
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
ZIP = os.path.join(HERE, "submission_af3_glide_select.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")


def glide_scores():
    """gscore per (ligand, sample); the 10000.0 sentinel can carry status ok."""
    out = {}
    for r in csv.DictReader(open(os.path.join(
            S04, "glide_rescoring", "results", "glide_scores.csv"))):
        if r["status"] != "ok" or not r["gscore"]:
            continue
        v = float(r["gscore"])
        if abs(v) > 100:
            continue
        out[(r["ligand"], r["sample"])] = v
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

    gsc = glide_scores()
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
                                         "gscore": gsc[(r["ligand"], r["sample"])]})
    print(f"{sum(len(v) for v in passing.values())} PB-passing prepped samples, "
          f"{len(passing)}/20 ligands covered")

    sel = []
    for n in sorted(smiles):
        pool = passing.get(n, [])
        if not pool:
            raise RuntimeError(f"{n}: no PoseBusters-passing prepped sample")
        best = min(pool, key=lambda x: x["gscore"])   # most negative wins
        alt_pl = max(pool, key=lambda x: x["pl_iptm"])["sample"]
        alt_mean = max(pool, key=lambda x: x["mean_pl_hl"])["sample"]
        with open(best["path"]) as f, open(os.path.join(OUT, f"{n}.pdb"), "w") as g:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    if (l[76:78].strip() or l[12:16].strip()[0]) == "H":
                        continue
                g.write(l)
        sel.append({"ligand": n, "n_passing": len(pool), "chosen": best["sample"],
                    "gscore": best["gscore"],
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
              f"gscore {best['gscore']:.3f}, PL {best['pl_iptm']:.2f})  "
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
