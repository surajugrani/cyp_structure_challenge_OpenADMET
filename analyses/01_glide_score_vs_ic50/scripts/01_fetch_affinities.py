"""
Step 1: Fetch binding affinities and ligand info from the RCSB PDB for every entry
in the RCSB custom report.

The custom report CSV is hard to parse reliably (multi-row entries, ligands with no
affinity, affinity rows not tied to a ligand). Instead we take only the entry IDs from
it and query the RCSB Data API, which links each affinity value to a specific ligand
(comp_id) and its source database.

Usage (from the analysis folder):
    $SCHRODINGER/run python3 scripts/01_fetch_affinities.py

Outputs:
    data/rcsb_raw.json        raw API response (for provenance)
    data/affinities_all.csv   one row per affinity value (all types: IC50/Kd/Ki/...)
    data/entries.csv          one row per entry x non-polymer ligand (+ resolution, chains)
"""

import csv
import json
import os

import pandas as pd
import requests

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
REPORT_CSV = os.path.join(ROOT, "input", "rcsb_pdb_custom_report_20260926222041.csv")
DATA_DIR = os.path.join(HERE, "data")

GRAPHQL_URL = "https://data.rcsb.org/graphql"
BATCH_SIZE = 50

QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_entry_info { resolution_combined experimental_method }
    rcsb_binding_affinity {
      comp_id type value unit symbol provenance_code reference_sequence_identity link
    }
    nonpolymer_entities {
      pdbx_entity_nonpoly { comp_id name }
      nonpolymer_comp { chem_comp { formula formula_weight } }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers { auth_asym_id auth_seq_id }
      }
    }
  }
}
"""


def read_entry_ids(report_csv):
    """Entry IDs are the non-empty values in the first column after the 2 header rows."""
    with open(report_csv) as f:
        rows = list(csv.reader(f))[2:]
    ids = [r[0].strip() for r in rows if r and r[0].strip()]
    return list(dict.fromkeys(ids))  # de-duplicate, keep order


def query_rcsb(ids):
    entries = []
    for i in range(0, len(ids), BATCH_SIZE):
        batch = ids[i:i + BATCH_SIZE]
        r = requests.post(GRAPHQL_URL, json={"query": QUERY, "variables": {"ids": batch}},
                          timeout=120, verify=True)
        r.raise_for_status()
        payload = r.json()
        if "errors" in payload:
            raise RuntimeError(payload["errors"])
        entries.extend(payload["data"]["entries"])
        print(f"  fetched {len(entries)}/{len(ids)}")
    return entries


def flatten(entries):
    aff_rows, lig_rows = [], []
    for e in entries:
        pdb = e["rcsb_id"]
        info = e.get("rcsb_entry_info") or {}
        res = (info.get("resolution_combined") or [None])[0]
        method = info.get("experimental_method")

        for a in e.get("rcsb_binding_affinity") or []:
            aff_rows.append({"pdb_id": pdb, "resolution": res, **a})

        for ent in e.get("nonpolymer_entities") or []:
            comp = ent["pdbx_entity_nonpoly"]["comp_id"]
            chem = (ent.get("nonpolymer_comp") or {}).get("chem_comp") or {}
            instances = [
                (inst["rcsb_nonpolymer_entity_instance_container_identifiers"]["auth_asym_id"],
                 inst["rcsb_nonpolymer_entity_instance_container_identifiers"]["auth_seq_id"])
                for inst in ent.get("nonpolymer_entity_instances") or []
            ]
            lig_rows.append({
                "pdb_id": pdb,
                "resolution": res,
                "method": method,
                "comp_id": comp,
                "name": ent["pdbx_entity_nonpoly"]["name"],
                "formula": chem.get("formula"),
                "formula_weight": chem.get("formula_weight"),
                "n_instances": len(instances),
                # chain:resnum for each copy, e.g. "A:601;B:601"
                "instances": ";".join(f"{c}:{s}" for c, s in instances),
            })
    return pd.DataFrame(aff_rows), pd.DataFrame(lig_rows)


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    ids = read_entry_ids(REPORT_CSV)
    print(f"{len(ids)} entry IDs read from {os.path.basename(REPORT_CSV)}")

    entries = query_rcsb(ids)
    with open(os.path.join(DATA_DIR, "rcsb_raw.json"), "w") as f:
        json.dump(entries, f, indent=1)

    aff, ligs = flatten(entries)
    aff.to_csv(os.path.join(DATA_DIR, "affinities_all.csv"), index=False)
    ligs.to_csv(os.path.join(DATA_DIR, "entries.csv"), index=False)

    print(f"{len(aff)} affinity values across {aff['pdb_id'].nunique()} entries")
    print(aff["type"].value_counts().to_string())
    print("Wrote data/rcsb_raw.json, data/affinities_all.csv, data/entries.csv")


if __name__ == "__main__":
    main()
