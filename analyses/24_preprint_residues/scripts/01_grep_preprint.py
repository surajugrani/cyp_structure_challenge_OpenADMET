"""Search the cryo-EM preprint for active-site residue mentions.

Question that prompted this: many ligands appear to make a polar contact with Ser119 —
does the paper discuss it?

A hand-rolled zlib/regex extractor returned 23 MB of noise and zero matches even for
"ketoconazole", because the text is CID-font encoded. pypdf (in cyp_posebusters_env)
decodes it properly.
"""
import re, sys
from pypdf import PdfReader

PDF = ("input/Structural Plasticity and Ligand Promiscuity of CYP3A4 "
       "Revealed by Cryo-EM.pdf")
r = PdfReader(PDF)
pages = [p.extract_text() or "" for p in r.pages]
full = "\n".join(pages)
flat = re.sub(r"\s+", " ", full)
print(f"{len(r.pages)} pages, {len(flat)} characters extracted")
print(f"sanity - 'ketoconazole' appears {len(re.findall('ketoconazole', flat, re.I))} times\n")

TARGETS = ["Ser ?-?119", "S119", "Arg ?-?105", "R105", "Arg ?-?212",
           "Phe ?-?304", "Ala ?-?305", "Thr ?-?309", "Glu ?-?374",
           "Phe ?-?215", "Phe ?-?108", "Phe ?-?213", "Leu ?-?211",
           "Ile ?-?369", "Ser ?-?119", "F-G loop", "substrate recognition"]
seen = set()
for t in TARGETS:
    if t in seen:
        continue
    seen.add(t)
    hits = list(re.finditer(t, flat, re.I))
    print(f"=== {t}: {len(hits)} match(es)")
    for m in hits[:4]:
        a, b = max(0, m.start() - 300), min(len(flat), m.end() + 300)
        print(f"    ...{flat[a:b]}...\n")

# every "<ThreeLetter><number>" style residue mention, ranked
res = re.findall(r"\b(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|"
                 r"Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s?-?\s?(\d{2,3})\b", flat)
from collections import Counter
c = Counter(f"{a}{b}" for a, b in res)
print("most-mentioned residues in the paper:")
for k, v in c.most_common(25):
    print(f"   {k:8} {v}")
