"""
Step 2: Curate the IC50 dataset: one structure and one pIC50 per ligand.

Rules (see ../README.md for the reasoning):
  1. IC50 only (Kd/Ki dropped: different quantities, too few to use).
  2. Drop qualified values ('<', '>'): they are bounds, not measurements.
  3. Drop entries whose ligand of interest contains iridium (no docking parameters).
  4. Keep values with reference_sequence_identity >= 95 %.
  5. pIC50 = -log10(IC50 [M]); per ligand (comp_id) take the MEDIAN pIC50.
     Values are pooled across entries that contain the same ligand, de-duplicated
     so the same BindingDB value attached to two entries is counted once.
  6. One structure per ligand: best (lowest) resolution.
  7. Flags (not auto-removed, for manual review):
       spread_flag   : max - min pIC50 > 1.5 log units
       shared_value  : identical IC50 set assigned to >1 different ligand
                       (e.g. stereoisomers given the same BindingDB value)
     Among ligands sharing the same value, only the best-resolution one is kept in
     the main set (`in_main_set`); the rest are kept in the table for reference.

Usage (from the analysis folder):
    $SCHRODINGER/run python3 scripts/02_curate.py

Outputs:
    data/ic50_values_used.csv   every IC50 value that went into a median
    data/curated_set.csv        one row per ligand: structure, chain, pIC50, flags
"""

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(HERE, "data")

MIN_SEQ_IDENTITY = 95
SPREAD_CUTOFF = 1.5  # log units


def main():
    aff = pd.read_csv(os.path.join(DATA_DIR, "affinities_all.csv"))
    ligs = pd.read_csv(os.path.join(DATA_DIR, "entries.csv"))

    # ── Filters ──────────────────────────────────────────────────────────────
    n0 = len(aff)
    aff = aff[aff["type"] == "IC50"]
    print(f"IC50 only: {n0} -> {len(aff)} values")

    n0 = len(aff)
    aff = aff[aff["symbol"].isna()]
    print(f"Drop qualified (<, >): {n0} -> {len(aff)}")

    n0 = len(aff)
    aff = aff[aff["reference_sequence_identity"] >= MIN_SEQ_IDENTITY]
    print(f"Seq identity >= {MIN_SEQ_IDENTITY}%: {n0} -> {len(aff)}")

    ir_entries = set(ligs.loc[ligs["formula"].str.contains(r"\bIr\b", na=False), "pdb_id"])
    n0 = len(aff)
    aff = aff[~aff["pdb_id"].isin(ir_entries)]
    print(f"Drop Ir-containing entries ({len(ir_entries)} in report): {n0} -> {len(aff)}")

    assert (aff["unit"] == "nM").all(), "unexpected unit"
    aff = aff.assign(pIC50=9 - np.log10(aff["value"]))

    # ── Pool values per ligand, counting each distinct value once per entry ──
    # The same BindingDB record is attached to every entry containing the ligand,
    # so take the multiset of values from the entry that lists the most values.
    per_entry = aff.groupby(["comp_id", "pdb_id"])["value"].apply(lambda v: sorted(v))
    values_used = []
    for comp, grp in per_entry.groupby(level=0):
        src_pdb, vals = max(grp.items(), key=lambda kv: len(kv[1]))
        for v in vals:
            values_used.append({"comp_id": comp, "source_pdb": src_pdb[1], "IC50_nM": v,
                                "pIC50": 9 - np.log10(v)})
    values_used = pd.DataFrame(values_used)
    values_used.to_csv(os.path.join(DATA_DIR, "ic50_values_used.csv"), index=False)

    stats = values_used.groupby("comp_id").agg(
        n_values=("pIC50", "size"),
        pIC50_median=("pIC50", "median"),
        pIC50_min=("pIC50", "min"),
        pIC50_max=("pIC50", "max"),
        IC50_values_nM=("IC50_nM", lambda v: ";".join(f"{x:g}" for x in sorted(v))),
    )
    stats["pIC50_spread"] = stats["pIC50_max"] - stats["pIC50_min"]
    stats["spread_flag"] = stats["pIC50_spread"] > SPREAD_CUTOFF

    # ── One structure per ligand: best resolution ────────────────────────────
    cand = ligs[ligs["comp_id"].isin(stats.index) & ligs["pdb_id"].isin(aff["pdb_id"])]
    cand = cand.sort_values(["comp_id", "resolution", "pdb_id"])
    all_pdbs = cand.groupby("comp_id")["pdb_id"].apply(lambda s: ";".join(s))
    best = cand.drop_duplicates("comp_id").set_index("comp_id")

    cur = best[["pdb_id", "resolution", "name", "formula", "formula_weight",
                "n_instances", "instances"]].join(stats)
    cur["all_pdbs_with_ligand"] = all_pdbs
    # Ligand copy to use: lowest chain ID (A if present), lowest resnum; protein chain =
    # same chain. Count copies in that chain (e.g. 2V0M has 2 stacked ketoconazoles).
    def pick_copy(instances):
        copies = sorted((c, int(s)) for c, s in (x.split(":") for x in instances.split(";")))
        chain, resnum = copies[0]
        return pd.Series({"lig_chain": chain, "lig_resnum": resnum,
                          "lig_copies_in_chain": sum(c == chain for c, _ in copies)})

    cur = cur.join(cur["instances"].apply(pick_copy))

    # ── Flag ligands that share an identical IC50 value set ─────────────────
    cur["shared_value"] = cur.duplicated("IC50_values_nM", keep=False)
    cur = cur.sort_values(["IC50_values_nM", "resolution"])
    cur["in_main_set"] = ~cur.duplicated("IC50_values_nM", keep="first")

    cur = cur.reset_index().sort_values("pIC50_median", ascending=False)
    cols = ["comp_id", "pdb_id", "resolution", "lig_chain", "lig_resnum",
            "lig_copies_in_chain", "name",
            "formula", "formula_weight", "n_values", "IC50_values_nM", "pIC50_median",
            "pIC50_spread", "spread_flag", "shared_value", "in_main_set",
            "all_pdbs_with_ligand"]
    cur[cols].to_csv(os.path.join(DATA_DIR, "curated_set.csv"), index=False, float_format="%.3f")

    print(f"\n{len(cur)} unique ligands; {cur['in_main_set'].sum()} in main set")
    print(cur[["comp_id", "pdb_id", "resolution", "lig_chain", "lig_copies_in_chain",
               "n_values", "pIC50_median",
               "pIC50_spread", "spread_flag", "shared_value", "in_main_set"]]
          .to_string(index=False, float_format="%.2f"))
    print("\nWrote data/ic50_values_used.csv, data/curated_set.csv")


if __name__ == "__main__":
    main()
