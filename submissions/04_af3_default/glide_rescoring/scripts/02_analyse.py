"""Does Glide scoring of AF3 samples select better poses, and does it track AF3's metrics?

Joins the Glide terms for each of the 100 AF3 samples to that sample's PoseBusters result
and AF3 confidences (analysis 10).

Two questions:
 1. If we pick the best-Glide-scoring sample per ligand, how many pass PoseBusters?
    Benchmarks: random 11.4/20, best AF3 metric 14/20, PoseBusters-filtered 19/20.
 2. Do Glide scores correlate with the AF3 metrics? Analysis 10 showed those metrics carry
    no information about pose validity, so if Glide agrees with them it is measuring the
    same non-signal.
"""
import csv, os, itertools
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
G = os.path.join(HERE, "results", "glide_scores.csv")
S = os.path.join(ROOT, "analyses", "10_af3_selection_vs_posebusters",
                 "results", "sample_level.csv")

GTERMS = ["gscore", "emodel", "evdw", "ecoul", "lipo", "einternal", "energy"]
AF3M = ["pl_iptm", "iptm", "ranking_score", "chain_iptm_lig"]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    g = {(r["ligand"], r["sample"]): r for r in csv.DictReader(open(G))
         if r["status"] == "ok"}
    rows = []
    for r in csv.DictReader(open(S)):
        k = (r["ligand"], r["sample"])
        if k not in g:
            continue
        rec = {"ligand": r["ligand"], "sample": r["sample"],
               "pb_pass": r["pb_pass"] == "True"}
        for m in AF3M:
            rec[m] = f(r[m])
        for t in GTERMS:
            rec[t] = f(g[k][t])
        rows.append(rec)
    print(f"{len(rows)} samples with both Glide scores and PoseBusters results\n")

    byl = {}
    for r in rows:
        byl.setdefault(r["ligand"], []).append(r)

    print("=== 1. Selecting by each Glide term: how many of 20 pass PoseBusters? ===")
    print(f"{'rule':24}{'PB passes':>12}")
    for t in GTERMS:
        n = 0
        for lig, v in byl.items():
            vv = [x for x in v if x[t] is not None]
            if vv:
                n += min(vv, key=lambda x: x[t])["pb_pass"]   # lower = better
        print(f"  best {t:18}{f'{n}/{len(byl)}':>10}")
    print(f"  {'(random baseline)':22}{'11.4/20':>10}")
    print(f"  {'(best AF3 metric)':22}{'14/20':>10}")
    print(f"  {'(PoseBusters filter)':22}{'19/20':>10}")

    print("\n=== 2. Glide term vs PoseBusters pass (pooled, point-biserial) ===")
    for t in GTERMS:
        v = [(x[t], x["pb_pass"]) for x in rows if x[t] is not None]
        a = np.array([p for p, _ in v]); b = np.array([q for _, q in v], dtype=float)
        r_, p_ = stats.pointbiserialr(b, a)
        print(f"  {t:12} r={r_:+.3f} p={p_:.3f}   mean(pass)={a[b==1].mean():8.2f} "
              f"mean(fail)={a[b==0].mean():8.2f}")

    print("\n=== 3. Glide vs AF3 metrics (Spearman, pooled over samples) ===")
    print(f"{'':12}" + "".join(f"{m:>16}" for m in AF3M))
    for t in GTERMS:
        line = f"  {t:10}"
        for m in AF3M:
            v = [(x[t], x[m]) for x in rows if x[t] is not None and x[m] is not None]
            if len(v) > 5:
                rho = stats.spearmanr([a for a, _ in v], [b for _, b in v]).statistic
                line += f"{rho:>+16.3f}"
            else:
                line += f"{'n/a':>16}"
        print(line)

    with open(os.path.join(HERE, "results", "joined.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\nWrote {HERE}/results/joined.csv")


if __name__ == "__main__":
    main()
