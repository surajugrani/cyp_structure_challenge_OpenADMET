"""Same two selection rules as 01, but on a single GlideScore TERM instead of gscore.

Term is given on the command line (default evdw). Rule 1 selects the sample with the best
(most negative) value of that term; rule 2 keeps submission 06's choice - the best PL-iPTM
PoseBusters-passing sample - and reports that sample's term value.

Original 01 docstring follows.

Glide score of the AF3 poses vs experimental pIC50, under two selection rules.

(1) best Glide score per ligand - the most favourable gscore among that ligand's
    PrepWizard-minimised AF3 samples;
(2) best PL-iPTM per ligand - the sample submission 06 would ship (06 is the best
    leaderboard result so far at 0.4291 LDDT-PLI), scored by Glide.

Both plot the same y (CYP3A4_pIC50_direct_inhibition) against the same x (gscore); only
which sample supplies the score changes. Rule 2 selects among PoseBusters-PASSING samples,
matching what submission 06 actually did; rule 1 is reported over all successfully scored
samples and, for comparison, over the PB-passing subset.

PRIOR: expect nothing. Analysis 01 found Glide crystal-pose score vs pIC50 rho = -0.27,
analysis 08 the same for ICM, and ~30 correlations have now been tested against these same
15 pIC50 values with one nominal hit traced to multiple testing. n = 15 over a ~2 log
range is the regime those nulls came from. The point of this run is the COMPARISON between
selection rules, not a hit.
"""
import csv, glob, json, math, os, sys
from collections import defaultdict

import numpy as np


def _rank(a):
    """Average ranks, so ties are handled the way Spearman requires."""
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    r[order] = np.arange(1, len(a) + 1, dtype=float)
    _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    for i in np.where(cnt > 1)[0]:
        m = inv == i
        r[m] = r[m].mean()
    return r


def _pearson(x, y):
    """r and a two-sided p from the t approximation; scipy is not in this env."""
    n = len(x)
    r = float(np.corrcoef(x, y)[0, 1])
    if n < 3 or abs(r) >= 1.0:
        return r, float("nan")
    t = abs(r) * np.sqrt((n - 2) / (1 - r * r))
    df = n - 2
    # two-sided tail of Student's t by continued-fraction incomplete beta
    xb = df / (df + t * t)
    a, b = df / 2.0, 0.5
    lbeta = (math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))
    if xb < (a + 1) / (a + b + 2):
        cf, flip = _betacf(a, b, xb), False
        ib = np.exp(a * np.log(xb) + b * np.log1p(-xb) - lbeta) * cf / a
    else:
        cf, flip = _betacf(b, a, 1 - xb), True
        ib = 1 - np.exp(b * np.log1p(-xb) + a * np.log(xb) - lbeta) * cf / b
    return r, float(min(1.0, max(0.0, ib)))


def _betacf(a, b, x, itmax=200, eps=3e-12):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-30 else 1e-30)
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > 1e-30 else 1e-30)
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > 1e-30 else 1e-30)
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def spearman(x, y):
    return _pearson(_rank(x), _rank(y))


def pearson(x, y):
    return _pearson(x, y)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
S04 = os.path.join(ROOT, "submissions", "04_af3_default")
TERM = sys.argv[1] if len(sys.argv) > 1 else "evdw"
GS = os.path.join(S04, "glide_rescoring", "results", "glide_scores.csv")
PBC = os.path.join(S04, "glide_rescoring", "results", "pb_after_prepwizard.csv")
AF3 = os.path.join(S04, "results")
IC50 = os.path.join(ROOT, "input", "cyp-challenge_cyp-inhibition-data.csv")
OUTCSV = "term_vs_ic50_{}.csv"


def pl_iptm(name, tag):
    g = glob.glob(os.path.join(AF3, name.lower(), tag, "*summary_confidences.json"))
    if not g:
        return None
    return json.load(open(g[0]))["chain_pair_iptm"][0][2]


def main():
    pic50 = {}
    for r in csv.DictReader(open(IC50)):
        v = r.get("CYP3A4_pIC50_direct_inhibition", "").strip()
        n = r.get("Molecule_Name", "").strip()
        if n and v:
            try:
                pic50[n] = float(v)
            except ValueError:
                pass

    pbpass = {(r["ligand"], r["sample"]) for r in csv.DictReader(open(PBC))
              if r["after_pass"] == "True"}

    rows = defaultdict(list)
    for r in csv.DictReader(open(GS)):
        if r["status"] != "ok" or not r["gscore"] or not r.get(TERM):
            continue
        # Glide writes 10000.0 as its "no score" sentinel; one sample of
        # OCNT-2395789 carries it with status ok, so filter on the value too.
        if abs(float(r["gscore"])) > 100 or abs(float(r[TERM])) > 1000:
            continue
        lig, tag = r["ligand"], r["sample"]
        if lig not in pic50:
            continue
        rows[lig].append({"sample": tag, "gscore": float(r[TERM]),
                          "pb": (lig, tag) in pbpass, "pl": pl_iptm(lig, tag)})

    print(f"TERM = {TERM}")
    print(f"{len(pic50)} ligands with CYP3A4 pIC50; "
          f"{len(rows)} of them have Glide-scored AF3 samples\n")

    out = []
    for lig in sorted(rows):
        pool = rows[lig]
        pbp = [x for x in pool if x["pb"]]
        best_g_all = min(pool, key=lambda x: x["gscore"])
        best_g_pb = min(pbp, key=lambda x: x["gscore"]) if pbp else None
        withpl = [x for x in (pbp or pool) if x["pl"] is not None]
        best_pl = max(withpl, key=lambda x: x["pl"]) if withpl else None
        out.append({
            "ligand": lig, "pIC50": pic50[lig], "n_samples": len(pool),
            "n_pb_passing": len(pbp),
            "rule1_best_gscore_sample": best_g_all["sample"],
            "rule1_best_gscore": best_g_all["gscore"],
            "rule1_pb_sample": best_g_pb["sample"] if best_g_pb else "",
            "rule1_pb_gscore": best_g_pb["gscore"] if best_g_pb else "",
            "rule2_best_pliptm_sample": best_pl["sample"] if best_pl else "",
            "rule2_pl_iptm": best_pl["pl"] if best_pl else "",
            "rule2_gscore": best_pl["gscore"] if best_pl else "",
        })

    with open(os.path.join(HERE, "results", OUTCSV.format(TERM)), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

    y = np.array([r["pIC50"] for r in out])
    series = [
        ("rule 1  best %s (all scored samples)" % TERM,
         np.array([r["rule1_best_gscore"] for r in out])),
        ("rule 1b best %s (PB-passing only)" % TERM,
         np.array([r["rule1_pb_gscore"] for r in out], dtype=float)),
        ("rule 2  best PL-iPTM sample, its %s" % TERM,
         np.array([r["rule2_gscore"] for r in out], dtype=float)),
    ]
    print(f"{'rule':48} {'n':>3} {'Spearman':>9} {'p':>7} {'Pearson':>8} {'p':>7}")
    for lab, x in series:
        m = ~np.isnan(x)
        rs, ps = spearman(x[m], y[m])
        rp, pp = pearson(x[m], y[m])
        print(f"{lab:48} {m.sum():3d} {rs:9.3f} {ps:7.3f} {rp:8.3f} {pp:7.3f}")

    same = sum(r["rule1_best_gscore_sample"] == r["rule2_best_pliptm_sample"]
               for r in out)
    print(f"\nrules 1 and 2 pick the same sample for {same}/{len(out)} ligands")
    print("\nper-ligand values written to results/glide_vs_ic50.csv")
    for r in out:
        print(f"  {r['ligand']}  pIC50 {r['pIC50']:.2f}  "
              f"rule1 {r['rule1_best_gscore']:7.3f}  rule2 "
              f"{r['rule2_gscore'] if r['rule2_gscore'] == '' else format(r['rule2_gscore'], '7.3f')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
