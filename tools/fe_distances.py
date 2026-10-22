"""Minimum haem-iron to ligand heavy-atom distance, read once and cached.

Every Fe-distance figure needs the same number: for one complex, the shortest distance
from the HEM iron to any heavy atom of the ligand. Before this module that routine was
written out separately for each source - shipped submission PDBs, pre-PrepWizard PDBs,
ESMFold2 mmCIF - and a one-off comparison meant writing a fourth copy in /tmp. One
parser now serves all of them.

Two file formats appear in this project and both are handled by `fe_min`:
  *.pdb   fixed-column; ligand residue is LIG, iron is element FE
  *.cif   mmCIF from a co-folding model; ligand is label_comp_id LIG in chain C

Results are cached in figures/.fe_cache.json keyed on (path, mtime, size), so
regenerating all the figures re-reads only structures that actually changed. Delete the
cache file to force a full recompute.
"""
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "figures", ".fe_cache.json")
_cache = None
_dirty = False


def _load():
    global _cache
    if _cache is None:
        try:
            _cache = json.load(open(CACHE))
        except (OSError, ValueError):
            _cache = {}
    return _cache


def save():
    """Write the cache back. Call once after a batch of reads."""
    global _dirty
    if _dirty and _cache is not None:
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        json.dump(_cache, open(CACHE, "w"))
        _dirty = False


def _from_pdb(path):
    fe, lig = None, []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        # Element column is authoritative; fall back to the atom name's first
        # character for writers that leave columns 77-78 blank.
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if el == "FE":
            fe = xyz
        elif l[17:20].strip() == "LIG":
            lig.append(xyz)
    return fe, lig


def _from_cif(path):
    lines = open(path).read().splitlines()
    cols = [l.strip().split(".", 1)[1] for l in lines if l.startswith("_atom_site.")]
    idx = {c: i for i, c in enumerate(cols)}
    need = max(idx[c] for c in ("type_symbol", "label_comp_id",
                                "Cartn_x", "Cartn_y", "Cartn_z"))
    fe, lig = None, []
    for l in lines:
        if not l.startswith(("ATOM", "HETATM")):
            continue
        r = l.split()
        if len(r) <= need:
            continue
        el = r[idx["type_symbol"]].upper()
        if el == "H":
            continue
        xyz = (float(r[idx["Cartn_x"]]), float(r[idx["Cartn_y"]]),
               float(r[idx["Cartn_z"]]))
        if el == "FE":
            fe = xyz
        elif r[idx["label_comp_id"]].upper() == "LIG":
            lig.append(xyz)
    return fe, lig


def fe_min(path):
    """Shortest Fe-to-ligand heavy-atom distance in Angstrom, or None if either is
    missing (a receptor-only file, or a complex whose ligand is not named LIG)."""
    global _dirty
    st = os.stat(path)
    key = f"{os.path.relpath(path, ROOT)}|{int(st.st_mtime)}|{st.st_size}"
    c = _load()
    if key in c:
        return c[key]
    fe, lig = _from_cif(path) if path.endswith(".cif") else _from_pdb(path)
    v = None if (fe is None or not lig) else round(
        min(math.dist(fe, x) for x in lig), 2)
    c[key] = v
    _dirty = True
    return v


def series(d, ext=None, strip=None):
    """{ligand: distance} for every structure in directory `d`.

    ext     restrict to this extension (default: .pdb and .cif both)
    strip   suffix to drop from the filename to recover the ligand name, for
            directories holding e.g. OCNT-1234567_pw.pdb
    """
    out = {}
    for f in sorted(os.listdir(d)):
        if ext and not f.endswith(ext):
            continue
        if not ext and not f.endswith((".pdb", ".cif")):
            continue
        name = f.rsplit(".", 1)[0]
        if strip and name.endswith(strip):
            name = name[: -len(strip)]
        elif strip:
            continue
        v = fe_min(os.path.join(d, f))
        if v is not None:
            out[name] = v
    save()
    return out
