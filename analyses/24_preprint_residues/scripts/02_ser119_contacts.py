"""Count Ser119 polar contacts in each submission's 20 complexes.

Why Ser119: the cryo-EM preprint reports it as the SECOND most universal active-site
contact across all CYP3A4 structures analysed - 90.9% of structures in which it was
modelled, behind only Ala305 at 92.6% and ahead of Phe304 (86.8%), Thr309 (80.2%),
Ile301 (76.0%), Phe108 (74.4%). All four of their cryo-EM ligand complexes contact it.
For azamulin it is the specific anchor: the pleuromutilin hydroxyl is coordinated by
S119 while the ring flips around the axis the hydroxyl defines.

This is therefore a benchmark grounded in experimental statistics rather than in a
docking score. If a submission's ligands contact Ser119 at a rate near 90%, its poses
behave like real CYP3A4 complexes on this axis.

Polar contact = any ligand N/O/S within 3.5 A of Ser119's OG (the hydroxyl) or, more
loosely, of any Ser119 polar atom (N, O, OG). Heavy-atom distances; no hydrogens needed.
Also reports the nearest ligand atom to OG regardless of element, to catch near-misses.
"""
import csv, math, os, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
SUBS = {s: os.path.join(ROOT, "submissions", s, "submission") for s in (
    "06_af3_pliptm", "16_boltz2_pliptm", "18_boltz2_glide",
    "19_hybrid_glide", "20_union_glide")}
CUT = 3.5


def parse(path):
    ser, lig = {}, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        name = l[12:16].strip()
        if l[17:20].strip() == "SER" and int(l[22:26]) == 119 and l[21] == "A":
            ser[name] = xyz
        elif l[17:20].strip() == "LIG":
            lig.append((name, el, xyz))
    return ser, lig


rows = []
for sub, d in SUBS.items():
    if not os.path.isdir(d):
        continue
    for f in sorted(os.listdir(d)):
        if not f.endswith(".pdb"):
            continue
        ser, lig = parse(os.path.join(d, f))
        if "OG" not in ser or not lig:
            continue
        og = ser["OG"]
        polar_ser = [v for k, v in ser.items() if k in ("OG", "O", "N")]
        d_og_any = min(math.dist(og, x) for _, _, x in lig)
        d_og_polar = min([math.dist(og, x) for _, e, x in lig if e in "NOS"],
                         default=float("inf"))
        d_anypolar = min(math.dist(s, x) for s in polar_ser
                         for _, e, x in lig if e in "NOS") if any(
                             e in "NOS" for _, e, _ in lig) else float("inf")
        rows.append({"submission": sub, "ligand": f[:-4],
                     "og_to_nearest_lig_atom": round(d_og_any, 2),
                     "og_to_nearest_lig_polar": None if math.isinf(d_og_polar)
                     else round(d_og_polar, 2),
                     "ser119polar_to_ligpolar": None if math.isinf(d_anypolar)
                     else round(d_anypolar, 2),
                     "hbond_og": d_og_polar <= CUT,
                     "contact_4A": d_og_any <= 4.0})

with open(os.path.join(ROOT, "analyses", "24_preprint_residues", "results",
                       "ser119_contacts.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

by = defaultdict(list)
for r in rows:
    by[r["submission"]].append(r)
print(f"Ser119 OG to nearest ligand polar atom, cutoff {CUT} A\n")
print(f"{'submission':20} {'H-bond to OG':>13} {'any contact <4A':>16} {'median OG dist':>15}")
import statistics as st
for s in SUBS:
    if s not in by:
        continue
    v = by[s]
    print(f"{s:20} {sum(r['hbond_og'] for r in v):10}/{len(v)} "
          f"{sum(r['contact_4A'] for r in v):13}/{len(v)} "
          f"{st.median(r['og_to_nearest_lig_atom'] for r in v):15.2f}")
print(f"\n(preprint: Ser119 contacted in 90.9% of CYP3A4 structures, 2nd of 68 "
      f"active-site residues)\n")
print(f"{'ligand':16} " + " ".join(f"{s.split('_')[0]:>6}" for s in SUBS if s in by))
ligs = sorted({r["ligand"] for r in rows})
for lg in ligs:
    cells = []
    for s in SUBS:
        if s not in by:
            continue
        m = [r for r in by[s] if r["ligand"] == lg]
        cells.append(f"{m[0]['og_to_nearest_lig_polar'] if m and m[0]['og_to_nearest_lig_polar'] is not None else '-':>6}")
    print(f"{lg:16} " + " ".join(str(c) for c in cells))
