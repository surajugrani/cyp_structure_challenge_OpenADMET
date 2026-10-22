"""Submission 24: per-ligand best Glide score across THREE co-folding models.

The one selection rule that has ever worked here is cross-model Glide. Submission 19
applied it to two models - AF3 (submission 06) and Boltz-2 (submission 18) - and beat
both parents, 0.4331 against 0.4291 and 0.4174. Analysis 26 now supplies the third
model's scores, so this is that same rule over {AF3, Boltz-2, ESMFold2}.

Taking the best score over {ESMFold2, submission 19} is identical to taking it over
{AF3, Boltz-2, ESMFold2}, because 19 is already the per-ligand min of the first two.
This script therefore starts from 19 and swaps in ESMFold2 where ESMFold2 wins, which
keeps the provenance obvious and reuses two already-validated structure sets.

Glide prefers ESMFold2 for 5 of 20:

    OCNT-2395415   -5.772 (Boltz)  ->  -5.826
    OCNT-2395436   -6.464 (Boltz)  ->  -7.354
    OCNT-2395464   -5.268 (Boltz)  ->  -6.065
    OCNT-2395503   -5.372 (AF3)    ->  -6.013
    OCNT-2395789   -6.581 (Boltz)  ->  -6.969     AF3 has no structure for this ligand

THE EVIDENCE IS MIXED AND THE POINT IS TO SETTLE IT.

Against: submission 23, pure ESMFold2, scored 0.4527 - our best, and well above
submission 19's 0.4331. Glide nonetheless picks ESMFold2 for only 5 of 20, so this
submission discards the better model on three quarters of the ligands. If Glide really
could rank structures across models it should have preferred ESMFold2 far more often,
and the fact that it does not is itself evidence against the rule.

For: submission 19 is the only selection rule that ever beat its parents, and it did so
on exactly this logic. Two models is a thin basis for believing it; three is a real
test. The user's purpose here is explicit - if this scores poorly, cross-model Glide
selection is not a capability we have, and we stop building on it.

Note the scores are not strictly commensurable: each is computed on that model's own
receptor and its own grid. Submission 19 relied on the same assumption and won, so the
assumption is being tested too.

Inputs are both already PrepWizard-minimised and validated 20/20, so nothing is
re-prepared; the structures are copied byte-for-byte and only re-validated.
"""
import csv
import os
import shutil
import sys
import zipfile
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S19 = os.path.join(ROOT, "submissions", "19_hybrid_glide", "submission")
S23 = os.path.join(ROOT, "submissions", "23_esmfold2_default", "submission")
A26 = os.path.join(ROOT, "analyses", "26_glide_on_esm", "results",
                   "glide_scores.csv")
SEL19 = os.path.join(ROOT, "submissions", "19_hybrid_glide", "results",
                     "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_hybrid_three_way.zip")
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
    nhem = len([l for l in at if l[17:20].strip() == "HEM"])
    if nhem != 43:
        e.append(f"HEM has {nhem} heavy atoms")
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
    os.makedirs(RES, exist_ok=True)
    tmp = os.path.join(RES, "_tmp")
    os.makedirs(tmp, exist_ok=True)
    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip()
              for r in csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}

    esm = {r["ligand"]: float(r["gscore"])
           for r in csv.DictReader(open(A26)) if r["status"] == "ok"}
    # submission 19 kept both parents' scores and which one it took
    prev = {r["ligand"]: r for r in csv.DictReader(open(SEL19))}

    rows = []
    for n in sorted(smiles):
        p = prev[n]
        cand = {}
        if p["gscore_06"]:
            cand["af3_06"] = float(p["gscore_06"])
        if p["gscore_18"]:
            cand["boltz_18"] = float(p["gscore_18"])
        if n in esm:
            cand["esmfold2_23"] = esm[n]
        if not cand:
            raise RuntimeError(f"{n}: no Glide score from any model")
        win = min(cand, key=cand.get)
        src = S23 if win == "esmfold2_23" else S19
        shutil.copy(os.path.join(src, f"{n}.pdb"), os.path.join(OUT, f"{n}.pdb"))
        rows.append({"ligand": n, "source": win,
                     "gscore_af3_06": p["gscore_06"],
                     "gscore_boltz_18": p["gscore_18"],
                     "gscore_esmfold2_23": esm.get(n, ""),
                     "gscore_best": round(cand[win], 3),
                     "sub19_took": p["source"],
                     "changed_vs_19": "yes" if win == "esmfold2_23" else "no"})
        if win == "esmfold2_23":
            print(f"  {n}: ESMFold2 wins at {cand[win]:.3f}   <-- SWAPPED")

    nsw = sum(r["changed_vs_19"] == "yes" for r in rows)
    print(f"\n{nsw}/20 taken from ESMFold2, {20 - nsw}/20 unchanged from submission 19")
    by = defaultdict(int)
    for r in rows:
        by[r["source"]] += 1
    print("  per model: " + ", ".join(f"{k} {v}" for k, v in sorted(by.items())))

    pb = PoseBusters(config="dock")
    ns = npb = 0
    for r in rows:
        p = os.path.join(OUT, f"{r['ligand']}.pdb")
        se = struct_errors(open(p).readlines(), r["ligand"],
                           smiles[r["ligand"]], tmp)
        pf = pb_fail(p, smiles[r["ligand"]], pb)
        r["struct_errors"] = "; ".join(se)
        r["pb_failed"] = "; ".join(pf)
        ns += not se
        npb += not pf
        if se or pf:
            print(f"  FINAL {r['ligand']}: struct {se} | PB {pf}")
    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nFINAL structural {ns}/20   PoseBusters {npb}/20")
    if ns == 20 and npb == 20:
        with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
            for r in rows:
                z.write(os.path.join(OUT, f"{r['ligand']}.pdb"),
                        f"{r['ligand']}.pdb")
        print(f"wrote {ZIP}")
    else:
        print("zip not written - failures above")
    return 0


if __name__ == "__main__":
    sys.exit(main())
