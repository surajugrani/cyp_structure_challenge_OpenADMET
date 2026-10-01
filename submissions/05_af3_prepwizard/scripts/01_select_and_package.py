"""Submission 05: AF3 + PrepWizard, selected on mean(protein-ligand, haem-ligand) iPTM.

Pipeline, all of it already computed upstream:
  * AF3 default co-folding, 5 samples per ligand   (submission 04)
  * PrepWizard on all 100 samples                  (submission 04 glide_rescoring, 9 h)
  * PoseBusters on all 100 prepped samples         -> 95/100 pass, all 20 ligands covered

Selection: among a ligand's PoseBusters-PASSING prepped samples, take the highest
mean(chain_pair_iptm[0][2], chain_pair_iptm[1][2]) - the protein-ligand and haem-ligand
interface confidences averaged. The first scores placement in the pocket, the second
scores the coordination geometry, which matters for a type-II CYP binder; global iptm is
not used because it averages in the near-perfect protein-haem interface (~0.97) and says
little about the blinded ligand.

Why PoseBusters is the filter rather than a confidence metric: analysis 10 measured every
AF3 metric against PoseBusters outcome and none predicted it (best 14/20 vs 11.4/20
random; has_clash is identically 0.0 across all 100). Glide scoring was tested too and
selected BELOW random (10/20), mildly preferring the poses PoseBusters rejects.

Also writes the alternative selections (pl_iptm only, haem-ligand only, AF3 ranking) to
the report so later attempts can be compared without recomputing anything.
"""
import csv, glob, json, os, sys, zipfile
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
AF3 = os.path.join(ROOT, "submissions", "04_af3_default", "results")
PWORK = os.path.join(ROOT, "submissions", "04_af3_default", "glide_rescoring", "work")
PBCSV = os.path.join(ROOT, "submissions", "04_af3_default", "glide_rescoring",
                     "results", "pb_after_prepwizard.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_af3_prepwizard.zip")


def conf(name, tag):
    f = glob.glob(os.path.join(AF3, name.lower(), tag, "*summary_confidences.json"))[0]
    j = json.load(open(f))
    cp = j.get("chain_pair_iptm")
    return {"pl_iptm": cp[0][2], "hl_iptm": cp[1][2],
            "mean_pl_hl": (cp[0][2] + cp[1][2]) / 2.0,
            "iptm": j.get("iptm"), "ranking_score": j.get("ranking_score")}


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    names = sorted(r["Molecule_Name"].strip() for r in csv.DictReader(
        open(os.path.join(ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
        if r["Molecule_Name"].strip())

    passing = defaultdict(list)
    for r in csv.DictReader(open(PBCSV)):
        if r["after_pass"] == "True":
            passing[r["ligand"]].append(r["sample"])

    rows = []
    for n in names:
        cands = passing.get(n, [])
        if not cands:
            raise RuntimeError(f"{n}: no PoseBusters-passing prepped sample")
        scored = [(t, conf(n, t)) for t in cands]
        tag, c = max(scored, key=lambda x: x[1]["mean_pl_hl"])

        src = os.path.join(PWORK, f"{n}__{tag}", f"{n}__{tag}_pw.pdb")
        if not os.path.exists(src):
            raise RuntimeError(f"{n}: prepped file missing: {src}")
        # Strip hydrogens: PrepWizard adds them (HEM goes 43 -> 73 atoms), while
        # submission 03 and the organizers' own sample predictions are heavy-atom only.
        # LDDT-PLI and BiSyRMSD are heavy-atom metrics, so nothing is lost.
        kept = 0
        with open(src) as f, open(os.path.join(OUT, f"{n}.pdb"), "w") as g:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    el = l[76:78].strip() or l[12:16].strip()[0]
                    if el == "H":
                        continue
                    kept += 1
                g.write(l)

        alt = {k: max(scored, key=lambda x: x[1][k])[0]
               for k in ("pl_iptm", "hl_iptm", "ranking_score")}
        rows.append({"ligand": n, "n_pb_passing": len(cands), "chosen": tag,
                     "pl_iptm": c["pl_iptm"], "hl_iptm": c["hl_iptm"],
                     "mean_pl_hl": round(c["mean_pl_hl"], 4),
                     "alt_best_pl_iptm": alt["pl_iptm"],
                     "alt_best_hl_iptm": alt["hl_iptm"],
                     "alt_best_ranking": alt["ranking_score"]})
        print(f"  {n}: {tag} from {len(cands)} passing  "
              f"(PL {c['pl_iptm']:.2f}, HL {c['hl_iptm']:.2f}, mean {c['mean_pl_hl']:.3f})")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)}/20 selected. Mean of chosen mean_pl_hl: "
          f"{sum(r['mean_pl_hl'] for r in rows)/len(rows):.3f}")
    same = sum(r["chosen"] == r["alt_best_pl_iptm"] for r in rows)
    print(f"selection agrees with pl_iptm-only in {same}/20")
    return 0


if __name__ == "__main__":
    sys.exit(main())
