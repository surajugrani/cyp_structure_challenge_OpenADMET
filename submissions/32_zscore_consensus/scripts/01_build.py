"""Submission 32: z-score consensus of Glide `inplace` gscore and ICM RTCNN.

Same three-ESMFold2-run pool as submissions 28 and 29, a third rule over it.

    28  best ICM RTCNN                 0.4411 / 3.7722 / 0.7965
    29  best Glide inplace gscore      0.4638 / 3.7277 / 0.8001   <- project best
    23  best single parent             0.4527 / 3.8551 / 0.7985

The case for consensus is that the two scorers are nearly decorrelated: analysis 29
found they pick the same winning structure for only 10 of 20 ligands, Pearson between
them +0.373. Decorrelated errors are the precondition under which averaging two rankers
beats either alone, and RTCNN is not a junk input - analysis 09 measured its
within-ligand Spearman against true RMSD on crystals at +0.334, roughly four times
Glide gscore's +0.086.

The case against, stated plainly because it is the stronger one on this pool: the two
have already met prospectively here, and Glide won by 0.0227. Every structure consensus
changes relative to submission 29 is a case of RTCNN overruling Glide, so this imports
precisely the disagreements that submission 28 lost on.

Method. gscore and RTCNN are different units on different scales, so a raw sum is
meaningless. Each is z-scored WITHIN a ligand's three candidates - subtract that
ligand's mean, divide by its population standard deviation - which puts both on a
common scale before they are added. Both are "more negative is better", so the minimum
of the summed z wins. With three candidates the standard deviation is estimated from
three points and is therefore noisy; that is a real weakness of the method here and the
reason rank-averaging was tried first. Rank-averaging turned out to have no resolution
at n=3 (ranks 0/1/2 tie constantly) and collapsed onto Glide for 18 of 20, so it would
have tested nothing.

Measured: this changes 7 of 20 against submission 29, and for 3 of those - 
OCNT-2315379, OCNT-2395457 and OCNT-2535422 - it selects a structure that neither
scorer picks on its own. Those three are the only places in this project where a
combined rule has produced a genuinely new choice rather than siding with one input.

Chosen in the knowledge that the evidence leans against it. The alternative considered
was Glide-primary with RTCNN breaking ties below a 0.25 gscore margin, which changes
only 3 of 20 - all with Glide margins of 0.00-0.12, i.e. where Glide is indifferent.
That is the safer rule but it moves ~1-2 structures inside the ~10 that are scored,
which is unreadable. Between a safe unreadable test and a risky readable one, this is
the readable one.

COPY-ONLY: every structure is copied unchanged from a submission that already passed
the full gate 20/20, verified by sha256. PoseBusters is not re-run - a byte-identical
copy cannot return a different verdict.
"""
import collections
import csv
import os
import statistics
import sys
import zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from copy_check import copy_verified, report  # noqa: E402

SETS = {"23_default_msa": "23_esmfold2_default",
        "26_max_accuracy": "26_esmfold2_maxacc",
        "27_default_nomsa": "27_esmfold2_nomsa"}
GLIDE = os.path.join(ROOT, "analyses", "28_inplace_vs_mininplace", "results",
                     "scores_both_modes.csv")
RTCNN = os.path.join(ROOT, "analyses", "29_icm_on_esm", "results",
                     "icm_esm_scores.csv")
SUB29 = os.path.join(ROOT, "submissions", "29_glide_inplace_best_of_three",
                     "results", "selection.csv")
OUT = os.path.join(HERE, "submission")
RES = os.path.join(HERE, "results")
ZIP = os.path.join(HERE, "submission_zscore_consensus.zip")


def zscore(d):
    """Within-ligand z-scores. Population sd, so a flat set maps to all zeros."""
    vals = list(d.values())
    mean = statistics.mean(vals)
    sd = statistics.pstdev(vals)
    return {k: (0.0 if sd == 0 else (v - mean) / sd) for k, v in d.items()}


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RES, exist_ok=True)

    g = collections.defaultdict(dict)
    for r in csv.DictReader(open(GLIDE)):
        if r["status"] == "ok" and r["mode"] == "inplace" and r["gscore"]:
            g[r["ligand"]][r["set"]] = float(r["gscore"])
    t = collections.defaultdict(dict)
    for r in csv.DictReader(open(RTCNN)):
        if r["status"] == "ok" and r.get("RTCNN"):
            t[r["ligand"]][r["set"]] = float(r["RTCNN"])

    prev = {r["ligand"]: r["picked_set"] for r in csv.DictReader(open(SUB29))}

    ligands = sorted(set(g) & set(t))
    if len(ligands) != 20:
        raise SystemExit(f"expected 20 ligands in both scorers, found {len(ligands)}")

    rows, copies = [], []
    for n in ligands:
        sets = sorted(set(g[n]) & set(t[n]))
        if len(sets) != len(SETS):
            raise SystemExit(f"{n}: scored in {len(sets)} sets, expected {len(SETS)}")
        zg = zscore({s: g[n][s] for s in sets})
        zr = zscore({s: t[n][s] for s in sets})
        tot = {s: zg[s] + zr[s] for s in sets}
        win = min(tot, key=tot.get)
        gwin = min(sets, key=lambda s: g[n][s])
        rwin = min(sets, key=lambda s: t[n][s])
        src = os.path.join(ROOT, "submissions", SETS[win], "submission", f"{n}.pdb")
        copies.append(copy_verified(src, os.path.join(OUT, f"{n}.pdb"),
                                    ligand=n, source=win))
        rows.append({"ligand": n, "picked_set": win,
                     "glide_pick": gwin, "rtcnn_pick": rwin,
                     "z_total": round(tot[win], 3),
                     "z_glide": round(zg[win], 3), "z_rtcnn": round(zr[win], 3),
                     "scorers_agree": "yes" if gwin == rwin else "no",
                     "new_choice": "yes" if win not in (gwin, rwin) else "no",
                     "sub29_pick": prev.get(n, ""),
                     "changed_vs_sub29": "yes" if win != prev.get(n) else "no"})

    tally = collections.Counter(r["picked_set"] for r in rows)
    print("picked per set: " + ", ".join(f"{k} {v}" for k, v in sorted(tally.items())))
    print(f"scorers agreed outright on "
          f"{sum(1 for r in rows if r['scorers_agree'] == 'yes')}/20")
    print(f"differs from submission 29 (Glide) on "
          f"{sum(1 for r in rows if r['changed_vs_sub29'] == 'yes')}/20")
    print(f"differs from submission 23 alone on "
          f"{sum(1 for r in rows if r['picked_set'] != '23_default_msa')}/20")
    print(f"chose a set NEITHER scorer picked alone: "
          f"{sum(1 for r in rows if r['new_choice'] == 'yes')}/20")
    if not report(copies):
        raise SystemExit("copy verification failed - zip not written")

    with open(os.path.join(RES, "selection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for r in rows:
            z.write(os.path.join(OUT, f"{r['ligand']}.pdb"), f"{r['ligand']}.pdb")
    print(f"wrote {os.path.relpath(ZIP, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
