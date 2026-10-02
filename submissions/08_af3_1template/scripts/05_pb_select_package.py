"""Submission 07, step 4: PoseBusters every prepped sample, select, package, validate.

Same rule that produced submission 05, so the two are comparable and the only variable is
the templates: among a ligand's PoseBusters-PASSING prepped samples, take the highest
mean(chain_pair_iptm[0][2], chain_pair_iptm[1][2]) - protein-ligand and haem-ligand
interface confidence averaged. Hydrogens are stripped (PrepWizard adds them; LDDT-PLI and
BiSyRMSD are heavy-atom metrics) and the organizers' structural checks are re-run before
the zip is written.
"""
import csv, glob, json, os, sys, zipfile
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
AF3 = os.path.join(HERE, "results")
WORK = os.path.join(HERE, "work")
OUT = os.path.join(HERE, "submission")
ZIP = os.path.join(HERE, "submission_af3_1template.zip")
sys.path.insert(0, os.path.join(ROOT, "submissions", "01_boltz2_default",
                                "scripts", "organizer_code"))
from posebusters_utils import assign_bond_orders_from_template  # noqa: E402
from posebusters import PoseBusters  # noqa: E402
RDLogger.DisableLog("rdApp.*")


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
    os.makedirs(OUT, exist_ok=True)
    tmp = os.path.join(HERE, "results", "_tmp"); os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    pb = PoseBusters(config="dock")

    # 1. PoseBusters every prepped sample
    rows = []
    files = sorted(glob.glob(os.path.join(WORK, "*", "*_pw.pdb")))
    print(f"PoseBusters on {len(files)} prepped samples")
    for i, f in enumerate(files, 1):
        base = os.path.basename(f)[:-7]
        name, tag = base.split("__")
        failed = pb_fail(f, smiles[name], pb)
        pl, hl = conf(name, tag)
        rows.append({"ligand": name, "sample": tag, "path": f,
                     "pb_pass": not failed, "pb_failed": "; ".join(failed),
                     "pl_iptm": pl, "hl_iptm": hl, "mean_pl_hl": (pl + hl) / 2})
        if i % 20 == 0:
            print(f"  {i}/{len(files)}", flush=True)
    with open(os.path.join(HERE, "results", "sample_pb.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    byl = defaultdict(list)
    for r in rows:
        byl[r["ligand"]].append(r)
    npass = sum(r["pb_pass"] for r in rows)
    cover = sum(any(x["pb_pass"] for x in v) for v in byl.values())
    print(f"\n{npass}/{len(rows)} samples pass; {cover}/20 ligands have a passing sample")

    # 2. select, strip hydrogens, write
    sel = []
    for n in sorted(smiles):
        v = byl.get(n, [])
        passing = [x for x in v if x["pb_pass"]]
        pool = passing or v          # fall back to all if none pass
        if not pool:
            raise RuntimeError(f"{n}: no prepped samples at all")
        best = max(pool, key=lambda x: x["mean_pl_hl"])
        with open(best["path"]) as f, open(os.path.join(OUT, f"{n}.pdb"), "w") as g:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    if (l[76:78].strip() or l[12:16].strip()[0]) == "H":
                        continue
                g.write(l)
        sel.append({"ligand": n, "chosen": best["sample"], "from_passing": bool(passing),
                    "n_passing": len(passing), "pl_iptm": best["pl_iptm"],
                    "hl_iptm": best["hl_iptm"],
                    "mean_pl_hl": round(best["mean_pl_hl"], 4)})
        print(f"  {n}: {best['sample']} ({len(passing)} passing, "
              f"mean {best['mean_pl_hl']:.3f})"
              + ("" if passing else "   <-- NO passing sample, used best available"))
    with open(os.path.join(HERE, "results", "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sel[0])); w.writeheader(); w.writerows(sel)

    # 3. re-validate the final set
    ns = npb2 = 0
    for n in sorted(smiles):
        p = os.path.join(OUT, f"{n}.pdb")
        lines = open(p).readlines()
        se = struct_errors(lines, n, smiles[n], tmp)
        pf = pb_fail(p, smiles[n], pb)
        ns += not se
        npb2 += not pf
        if se or pf:
            print(f"  FINAL {n}: struct {se} | PB {pf}")
    print(f"\nFINAL structural {ns}/20   PoseBusters {npb2}/20")
    if ns == 20 and npb2 == 20:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for n in sorted(smiles):
                z.write(os.path.join(OUT, f"{n}.pdb"), f"{n}.pdb")
        print(f"wrote {ZIP}")
    else:
        print("zip not written - failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
