"""Fit pIC50 from the Glide vdW and electrostatic terms, honestly cross-validated.

Features: evdw and ecoul of one pose per ligand. Two selection rules, as in 01/02:
  rule 1 - the sample with the best (most negative) gscore
  rule 2 - the sample submission 06 ships, i.e. best PL-iPTM among PoseBusters-passing

Models, all fitted with numpy least squares (no sklearn in this env):
  linear      pIC50 ~ a*evdw + b*ecoul + c
  quadratic   adds evdw^2, ecoul^2, evdw*ecoul            (6 parameters)
  mean-only   the null model, predicts the training mean   (baseline)

THE STATISTICS MATTER MORE THAN THE FIT HERE. n = 15. A 3-parameter linear model on 15
points will always produce a non-zero in-sample R^2, and the 6-parameter quadratic will
produce a large one, purely from fitting noise. So the reported number is leave-one-out
cross-validated Q^2:

    Q^2 = 1 - sum((y_i - yhat_i_from_a_fit_that_excluded_i)^2) / sum((y_i - ybar)^2)

Q^2 <= 0 means the model predicts a held-out ligand no better than simply guessing the
mean pIC50 - i.e. it has learned nothing generalisable. Unlike R^2, Q^2 can go negative,
and for a model fitted to noise it usually does.

PRIOR: analyses 01 and 08 found no single-term correlation with these pIC50 values, and
02 here found evdw alone at rho = -0.03. A 2-feature combination of terms that are
individually flat, on 15 points over ~2 log units, is very unlikely to generalise.
"""
import csv, glob, json, os, sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S04 = os.path.join(ROOT, "submissions", "04_af3_default")
GS = os.path.join(S04, "glide_rescoring", "results", "glide_scores.csv")
PBC = os.path.join(S04, "glide_rescoring", "results", "pb_after_prepwizard.csv")
AF3 = os.path.join(S04, "results")
IC50 = os.path.join(ROOT, "input", "cyp-challenge_cyp-inhibition-data.csv")


def pl_iptm(name, tag):
    g = glob.glob(os.path.join(AF3, name.lower(), tag, "*summary_confidences.json"))
    return json.load(open(g[0]))["chain_pair_iptm"][0][2] if g else None


def design(X, kind):
    ones = np.ones((len(X), 1))
    if kind == "mean":
        return ones
    if kind == "linear":
        return np.hstack([ones, X])
    v, c = X[:, :1], X[:, 1:2]
    return np.hstack([ones, X, v * v, c * c, v * c])


def fit_predict(Xtr, ytr, Xte, kind):
    A = design(Xtr, kind)
    beta, *_ = np.linalg.lstsq(A, ytr, rcond=None)
    return design(Xte, kind) @ beta, beta


def evaluate(X, y, kind):
    n = len(y)
    insample, _ = fit_predict(X, y, X, kind)
    r2 = 1 - ((y - insample) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    loo = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        pred, _ = fit_predict(X[m], y[m], X[i:i + 1], kind)
        loo[i] = float(pred[0])
    q2 = 1 - ((y - loo) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    rmse = float(np.sqrt(((y - loo) ** 2).mean()))
    return r2, q2, rmse, insample, loo


def main():
    pic50 = {}
    for r in csv.DictReader(open(IC50)):
        n, v = r.get("Molecule_Name", "").strip(), r.get("CYP3A4_pIC50_direct_inhibition", "").strip()
        if n and v:
            try:
                pic50[n] = float(v)
            except ValueError:
                pass
    pbpass = {(r["ligand"], r["sample"]) for r in csv.DictReader(open(PBC))
              if r["after_pass"] == "True"}

    pool = defaultdict(list)
    for r in csv.DictReader(open(GS)):
        if r["status"] != "ok" or not r["gscore"] or abs(float(r["gscore"])) > 100:
            continue
        if r["ligand"] not in pic50:
            continue
        pool[r["ligand"]].append({
            "sample": r["sample"], "gscore": float(r["gscore"]),
            "evdw": float(r["evdw"]), "ecoul": float(r["ecoul"]),
            "pb": (r["ligand"], r["sample"]) in pbpass,
            "pl": pl_iptm(r["ligand"], r["sample"])})

    rules = {}
    for lig in sorted(pool):
        p = pool[lig]
        pbp = [x for x in p if x["pb"]] or p
        rules.setdefault("rule1_best_gscore", {})[lig] = min(p, key=lambda x: x["gscore"])
        withpl = [x for x in pbp if x["pl"] is not None]
        rules.setdefault("rule2_best_pliptm", {})[lig] = max(withpl, key=lambda x: x["pl"])

    rows = []
    print(f"n = {len(pool)} ligands with pIC50 and Glide-scored AF3 samples\n")
    print(f"{'rule':20} {'model':10} {'params':>6} {'R2 in-sample':>13} {'Q2 LOO':>9} {'RMSE LOO':>9}")
    for rname, sel in rules.items():
        ligs = sorted(sel)
        X = np.array([[sel[l]["evdw"], sel[l]["ecoul"]] for l in ligs], float)
        y = np.array([pic50[l] for l in ligs], float)
        for kind, npar in (("mean", 1), ("linear", 3), ("quadratic", 6)):
            r2, q2, rmse, ins, loo = evaluate(X, y, kind)
            print(f"{rname:20} {kind:10} {npar:6d} {r2:13.3f} {q2:9.3f} {rmse:9.3f}")
            if kind != "mean":
                for l, a, b, p_, lo in zip(ligs, X[:, 0], X[:, 1], ins, loo):
                    rows.append({"rule": rname, "model": kind, "ligand": l,
                                 "evdw": a, "ecoul": b, "pIC50": pic50[l],
                                 "pred_insample": round(float(p_), 4),
                                 "pred_loo": round(float(lo), 4)})
        _, beta = fit_predict(X, y, X, "linear")
        print(f"{'':20} linear coefficients: intercept {beta[0]:+.3f}  "
              f"evdw {beta[1]:+.4f}  ecoul {beta[2]:+.4f}")
        print(f"{'':20} pIC50 sd {y.std(ddof=1):.3f}\n")

    with open(os.path.join(HERE, "results", "model_fit.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("per-ligand predictions written to results/model_fit.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
