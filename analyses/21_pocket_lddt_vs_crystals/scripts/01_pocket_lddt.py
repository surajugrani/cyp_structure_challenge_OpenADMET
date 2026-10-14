"""Analysis 21 - how native-like is each predicted pocket, against real CYP3A4 crystals?

The question: the leaderboard says Boltz-2 (submission 03) scores LDDT-LP 0.7962 while
AF3 (submission 06) scores 0.6985 - a 0.098 gap on the protein pocket - and we have no
per-structure breakdown to explain it.

We cannot reproduce the organizers' LDDT-LP, because that needs the blinded reference
structures. But CYP3A4 is the same protein in every PDB entry, so a crystal pocket is a
legitimate yardstick for "does this look like a real CYP3A4 conformation".

WHY ALL THE CRYSTALS, NOT ONE. The CYP3A4 pocket is famously plastic; any single
reference bakes in that crystal's particular conformation. Each model is therefore scored
against EVERY filtered crystal, and two numbers are reported:

  max     the best-matching crystal - "does this pocket resemble ANY real conformation?"
          This is the meaningful one, since the blinded references are themselves some
          conformation from this family.
  median  how typical the pocket is across the ensemble.

The best-matching PDB id is recorded too: if Boltz and AF3 systematically match different
crystals, that localises the gap to a conformational state rather than to general quality.

REFERENCE FILTER: ligand-bound structures at <=2.5 A resolution (38 of 100). Apo pockets
are a different shape and poor maps add noise.

POCKET DEFINITION: a FIXED residue-number set - residues with a heavy atom within 10 A of
the haem iron in the highest-resolution reference. Fixed rather than per-model, so both
methods are scored on identical atoms and the comparison is not confounded by each
model's own ligand placement.

lDDT comes from OpenStructure (the organizers' own tool), in the pxr-eval environment.
"""
import csv, json, os, sys
from collections import defaultdict

import numpy as np
from ost import io
from ost.mol.alg.lddt import lDDTScorer

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
A05 = os.path.join(ROOT, "analyses", "05_full_pdb_redocking")
RESULTS = os.path.join(HERE, "results")
MODELS = {"boltz_03": os.path.join(ROOT, "submissions", "03_boltz2_prepwizard",
                                   "submission"),
          "af3_06": os.path.join(ROOT, "submissions", "06_af3_pliptm", "submission")}
MAX_RES = 2.5
POCKET_R = 10.0


def references():
    out = []
    for r in csv.DictReader(open(os.path.join(A05, "results", "redock_summary.csv"))):
        try:
            res = float(r["resolution"])
        except (TypeError, ValueError):
            continue
        if res > MAX_RES:
            continue
        p = os.path.join(A05, "structures", r["pdb_id"], f"{r['pdb_id']}.pdb")
        if os.path.exists(p):
            out.append((r["pdb_id"], res, p))
    return sorted(out, key=lambda x: x[1])


def pocket_resnums(path):
    """Residue numbers with a heavy atom within POCKET_R of the haem iron."""
    fe, atoms = None, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[:2]).upper()
        if el == "H":
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if el == "FE" and fe is None:
            fe = xyz
        if l.startswith("ATOM"):
            atoms.append((int(l[22:26]), xyz))
    if fe is None:
        return set()
    fe = np.array(fe)
    keep = set()
    for num, xyz in atoms:
        if np.linalg.norm(np.array(xyz) - fe) <= POCKET_R:
            keep.add(num)
    return keep


def load_pocket(path, resnums):
    ent = io.LoadPDB(path)
    sel = ent.Select("peptide=true and rnum=" + ",".join(str(n) for n in sorted(resnums)))
    return sel


def main():
    os.makedirs(RESULTS, exist_ok=True)
    refs = references()
    if not refs:
        raise SystemExit("no references passed the filter")
    pocket = pocket_resnums(refs[0][2])
    print(f"{len(refs)} reference crystals <= {MAX_RES} A; best is {refs[0][0]} "
          f"({refs[0][1]} A)")
    print(f"pocket = {len(pocket)} residues within {POCKET_R} A of Fe in {refs[0][0]}\n",
          flush=True)
    json.dump(sorted(pocket), open(os.path.join(RESULTS, "pocket_residues.json"), "w"))

    scorers = {}
    for pdb, res, path in refs:
        try:
            s = load_pocket(path, pocket)
            if s.residue_count < 10:
                continue
            scorers[pdb] = (lDDTScorer(s), res, s.residue_count)
        except Exception as e:
            print(f"  reference {pdb} unusable: {e}", flush=True)
    print(f"{len(scorers)} usable references\n", flush=True)

    rows, per_model = [], defaultdict(list)
    for mname, mdir in MODELS.items():
        for f in sorted(os.listdir(mdir)):
            if not f.endswith(".pdb"):
                continue
            lig = f[:-4]
            try:
                mdl = load_pocket(os.path.join(mdir, f), pocket)
            except Exception as e:
                print(f"{mname} {lig}: unreadable ({e})", flush=True)
                continue
            scores = {}
            for pdb, (sc, res, nres) in scorers.items():
                try:
                    v, _ = sc.lDDT(mdl, chain_mapping={"A": "A"}) \
                        if False else sc.lDDT(mdl)
                    if v is not None:
                        scores[pdb] = float(v)
                except Exception:
                    continue
            if not scores:
                print(f"{mname} {lig}: no score against any reference", flush=True)
                continue
            best = max(scores, key=scores.get)
            vals = sorted(scores.values())
            rec = {"model": mname, "ligand": lig,
                   "n_refs": len(scores),
                   "lddt_max": round(scores[best], 4),
                   "best_ref": best,
                   "best_ref_resolution": scorers[best][1],
                   "lddt_median": round(float(np.median(vals)), 4),
                   "lddt_min": round(vals[0], 4)}
            rows.append(rec)
            per_model[mname].append(rec)
            print(f"{mname:9} {lig}  max {rec['lddt_max']:.3f} ({best})  "
                  f"median {rec['lddt_median']:.3f}", flush=True)

    with open(os.path.join(RESULTS, "pocket_lddt.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print("\n" + "=" * 62)
    for m, rs in per_model.items():
        mx = [r["lddt_max"] for r in rs]
        md = [r["lddt_median"] for r in rs]
        print(f"{m:9} n={len(rs)}  max: mean {np.mean(mx):.4f} median {np.median(mx):.4f}"
              f"   median-vs-ensemble: mean {np.mean(md):.4f}")
    if len(per_model) == 2:
        a, b = list(per_model)
        pa = {r["ligand"]: r["lddt_max"] for r in per_model[a]}
        pb = {r["ligand"]: r["lddt_max"] for r in per_model[b]}
        common = sorted(set(pa) & set(pb))
        d = [pb[l] - pa[l] for l in common]
        print(f"\nper-ligand {b} - {a} (max lDDT): mean {np.mean(d):+.4f}, "
              f"{sum(1 for x in d if x > 0)}/{len(d)} favour {b}")
        for l in sorted(common, key=lambda x: -(pb[x] - pa[x]))[:5]:
            print(f"  {l}  {a} {pa[l]:.3f}  {b} {pb[l]:.3f}  diff {pb[l]-pa[l]:+.3f}")
        print("  ...")
        for l in sorted(common, key=lambda x: (pb[x] - pa[x]))[:5]:
            print(f"  {l}  {a} {pa[l]:.3f}  {b} {pb[l]:.3f}  diff {pb[l]-pa[l]:+.3f}")
        from collections import Counter
        for m in (a, b):
            c = Counter(r["best_ref"] for r in per_model[m])
            print(f"\n{m} best-matching crystals: {c.most_common(5)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
