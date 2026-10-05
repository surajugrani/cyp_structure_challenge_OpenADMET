"""Submission 09: AF3 + single 38FW holo template, selected on PL-iPTM ALONE.

Submission 08 picked, among a ligand's PoseBusters-passing PrepWizard-minimised samples,
the one with the highest mean(chain_pair_iptm[0][2], chain_pair_iptm[1][2]) - the
protein-ligand and haem-ligand interface confidences averaged. This run changes exactly
one thing: the key is chain_pair_iptm[0][2] on its own.

It is the same A/B that 05 -> 06 ran on the no-template predictions, where dropping the
haem term moved 2 of 20 structures and the score by +0.0012 LDDT-PLI. Here it moves
3 of 20: OCNT-2313300, OCNT-2395328 and OCNT-2395415.

Nothing is recomputed. AF3 inference (job 47869043), PrepWizard on all 100 samples
(job 47869044) and PoseBusters on the 80 that survived prep (job 47869045) all live in
submission 08; this script re-reads 08's results/sample_pb.csv, re-selects, strips
hydrogens, re-runs the organizers' structural checks plus PoseBusters on the final 20,
and writes the zip.

CAVEAT on interpreting the result: the leaderboard scores ~10 of the 20 complexes, and in
2 of the 3 differing cases the two candidate samples are tied in PL-iPTM to two decimals.
If none of the three fall in the scored subset the metrics will come back identical to
08's, which would show nothing either way - not that the two rules are equivalent.
"""
import csv, os, sys, zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S08 = os.path.join(ROOT, "submissions", "08_af3_1template")
PBCSV = os.path.join(S08, "results", "sample_pb.csv")
OUT = os.path.join(HERE, "submission")
ZIP = os.path.join(HERE, "submission_af3_1template_pliptm.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")

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
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    tmp = os.path.join(HERE, "results", "_tmp"); os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = PoseBusters(config="dock")

    # 1. reuse submission 08's per-sample PoseBusters table
    byl = defaultdict(list)
    for r in csv.DictReader(open(PBCSV)):
        r["pb_pass"] = r["pb_pass"].lower() in ("true", "1")
        r["pl_iptm"] = float(r["pl_iptm"]); r["hl_iptm"] = float(r["hl_iptm"])
        r["mean_pl_hl"] = float(r["mean_pl_hl"])
        byl[r["ligand"]].append(r)
    print(f"read {sum(len(v) for v in byl.values())} prepped samples from 08")

    # 2. select on pl_iptm alone, strip hydrogens, write
    sel = []
    for n in sorted(smiles):
        v = byl.get(n, [])
        passing = [x for x in v if x["pb_pass"]]
        pool = passing or v
        if not pool:
            raise RuntimeError(f"{n}: no prepped samples at all")
        best = max(pool, key=lambda x: x["pl_iptm"])
        alt08 = max(pool, key=lambda x: x["mean_pl_hl"])
        with open(best["path"]) as f, open(os.path.join(OUT, f"{n}.pdb"), "w") as g:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    if (l[76:78].strip() or l[12:16].strip()[0]) == "H":
                        continue
                g.write(l)
        sel.append({"ligand": n, "chosen": best["sample"], "from_passing": bool(passing),
                    "n_passing": len(passing), "pl_iptm": best["pl_iptm"],
                    "hl_iptm": best["hl_iptm"],
                    "mean_pl_hl": round(best["mean_pl_hl"], 4),
                    "s08_chose": alt08["sample"],
                    "differs_from_08": best["sample"] != alt08["sample"]})
        flag = "   <-- DIFFERS from 08" if best["sample"] != alt08["sample"] else ""
        print(f"  {n}: {best['sample']} ({len(passing)} passing, "
              f"PL {best['pl_iptm']:.2f}){flag}")
    with open(os.path.join(HERE, "results", "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sel[0])); w.writeheader(); w.writerows(sel)
    nd = sum(r["differs_from_08"] for r in sel)
    print(f"\n{nd}/20 structures differ from submission 08")

    # 3. re-validate the final set
    ns = npb = 0
    for n in sorted(smiles):
        p = os.path.join(OUT, f"{n}.pdb")
        lines = open(p).readlines()
        se = struct_errors(lines, n, smiles[n], tmp)
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
