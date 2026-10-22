"""Fold the 20 complexes again through ESMFold2's `fold_max_accuracy` endpoint.

Submission 23 used `fold_all_atom`, which is the ordinary endpoint: it takes a
FoldingConfig of quality knobs (num_sampling_steps, num_loops, lm_dropout, MSA
subsampling) and we ran it at defaults with our ColabFold MSA. It scored 0.4527
LDDT-PLI, our best.

`fold_max_accuracy` is a DIFFERENT endpoint on the same model family, described in the
SDK as "fold a molecular complex at the highest-accuracy settings". Three differences
matter and none of them are our choice:

  * It exposes no quality knobs at all - FoldMaxAccuracyConfig carries only the three
    include_* output flags. The server decides the settings.
  * It REFUSES a supplied MSA: validate_fold_max_accuracy_input raises on any protein
    chain with msa set, because "fold_max_accuracy generates its own MSA". So this run
    is necessarily MSA-free from our side, and partly answers whether our ColabFold MSA
    was pulling any weight - though not cleanly, since the endpoint changes more than
    just the MSA.
  * It pins its own model (ESMFOLD2_CUTOFF_2025) rather than the esmfold2-2026-05 we
    requested in submission 23, so the weights may differ too.

It is also a batch endpoint reached through ForgeBatchClient, not the inference client,
and its own default timeout is 2 hours per complex - it is expected to be much slower
than the ~2 min/complex of submission 23.

Everything else is held identical to submission 23: same sequence, same chain layout
(A protein / B HEM / C ligand), same SMILES source, mmCIF out.

Runs on the login node: this is an HTTPS client, the folding happens on Biohub's GPUs.

Usage:
    source ~/.cyp_biohub_token
    python 06_fold_max_accuracy.py --only OCNT-2395436
    python 06_fold_max_accuracy.py
"""
import argparse, csv, json, os, sys, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
RAW = os.path.join(HERE, "max_accuracy", "raw")
RES = os.path.join(HERE, "max_accuracy", "results")

def load_seq():
    """Reuse submission 23's sequence rather than copying it, so the two runs cannot
    silently diverge. The module name starts with a digit, so import_module is needed."""
    import importlib.util
    p = os.path.join(HERE, "scripts", "01_fold.py")
    spec = importlib.util.spec_from_file_location("fold23", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.SEQ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True); os.makedirs(RES, exist_ok=True)

    tok = os.environ.get("BIOHUB_TOKEN")
    if not tok:
        raise SystemExit("BIOHUB_TOKEN not set - run: source ~/.cyp_biohub_token")

    from esm.sdk import ForgeBatchClient
    from esm.sdk.api import FoldMaxAccuracyConfig
    from esm.utils.structure import input_builder as ib

    seq = load_seq()
    client = ForgeBatchClient(url="https://biohub.ai", token=tok)
    cfg = FoldMaxAccuracyConfig(include_pae=True, include_pair_chains_iptm=True)

    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    names = [a.only] if a.only else sorted(smiles)

    rows = []
    for i, n in enumerate(names, 1):
        out = os.path.join(RAW, f"{n}.cif")
        if os.path.exists(out):
            print(f"[{i}/{len(names)}] {n}: already folded", flush=True)
            continue
        print(f"[{i}/{len(names)}] {n} ...", flush=True)
        t0 = time.time()
        # msa deliberately omitted - the endpoint rejects a supplied one.
        inp = ib.StructurePredictionInput(sequences=[
            ib.ProteinInput(id="A", sequence=seq),
            ib.LigandInput(id="B", ccd=["HEM"]),
            ib.LigandInput(id="C", smiles=smiles[n]),
        ])
        r = client.fold_max_accuracy.run(all_atom_input=inp, config=cfg)
        if isinstance(r, list):
            print(f"    NOTE: got a list of {len(r)} results; using the first")
            r = r[0]
        if hasattr(r, "error") or type(r).__name__ == "ESMProteinError":
            print(f"    FAILED: {r}", flush=True)
            rows.append({"ligand": n, "seconds": round(time.time() - t0, 1),
                         "error": str(r)[:300]})
            continue
        with open(out, "w") as fh:
            fh.write(r.complex.to_mmcif())
        rec = {"ligand": n, "seconds": round(time.time() - t0, 1)}
        for k in ("ptm", "iptm"):
            v = getattr(r, k, None)
            if v is not None:
                rec[k] = float(v)
        pl = getattr(r, "plddt", None)
        if pl is not None:
            rec["plddt_mean"] = round(float(pl.mean()), 4)
        pc = getattr(r, "pair_chains_iptm", None)
        if pc is not None:
            rec["pair_chains_iptm"] = json.dumps(
                pc.tolist() if hasattr(pc, "tolist") else pc)
        rows.append(rec)
        print(f"    {rec}", flush=True)

        p = os.path.join(RES, "fold_metrics.csv")
        new = not os.path.exists(p)
        with open(p, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["ligand", "seconds", "ptm", "iptm",
                                               "plddt_mean", "pair_chains_iptm",
                                               "error"], extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerow(rows[-1])
    print(f"\ndone: {len(os.listdir(RAW))} cif files in {RAW}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
