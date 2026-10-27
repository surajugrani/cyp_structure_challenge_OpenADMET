"""Fold CYP3A4 + haem + blinded ligand with ESMFold2 (remote Biohub API).

ESMFold2 is the first method here that is NOT run locally - the folding happens on
Biohub's GPUs and this machine is only an HTTPS client. There is therefore no SLURM job
for the folding, and it must run somewhere with outbound internet (login node: yes;
compute nodes: usually not).

Chain layout matches every other submission so the organizers' validator is happy:
    A = protein (503 aa)
    B = haem, as CCD "HEM"
    C = the blinded ligand, as SMILES

Config: defaults apart from two flags that cost nothing and we need downstream -
include_pair_chains_iptm (the PL-iPTM analogue that produced submission 16's +0.0166)
and include_pae.

ESMFold2 returns mmCIF via result.complex.to_mmcif(); the organizers reject a zip
without .pdb files, so conversion happens later in the pipeline.

Usage:
    source ~/.cyp_biohub_token
    python 01_fold.py --only OCNT-2395436        # single trial
    python 01_fold.py                            # all 20
"""
import argparse, csv, json, os, sys, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(HERE))
RAW = os.path.join(HERE, "raw")
RES = os.path.join(HERE, "results")
SEQ = ("MALIPDLAMETWLLLAVSLVLLYLYGTHSHGLFKKLGIPGPTPLPFLGNILSYHKGFCMFDMECHKKYGKVWGFYDGQ"
       "QPVLAITDPDMIKTVLVKECYSVFTNRRPFGPVGFMKSAISIAEDEEWKRLRSLLSPTFTSGKLKEMVPIIAQYGDVLV"
       "RNLRREAETGKPVTLKDVFGAYSMDVITSTSFGVNIDSLNNPQDPFVENTKKLLRFDFLDPFFLSITVFPFLIPILEVL"
       "NICVFPREVTNFLRKSVKRMKESRLEDTQKHRVDFLQLMIDSQNSKETESHKALSDLELVAQSIIFIFAGYETTSSVLS"
       "FIMYELATHPDVQQKLQEEIDAVLPNKAPPTYDTVLQMEYLDMVVNETLRLFPIAMRLERVCKKDVEINGMFIPKGVVV"
       "MIPSYALHRDPKYWTEPEKFLPERFSKKNKDNIDPYIYTPFGSGPRNCIGMRFALMNMKLALIRVLQNFSFKPCKETQI"
       "PLKLSLGGLLQPEKPVVLKVESRDGTVSGA")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    ap.add_argument("--model", default="esmfold2-2026-05")
    ap.add_argument("--no-msa", action="store_true",
                    help="fold from sequence alone (the 21 s trial used this)")
    ap.add_argument("--sampling-steps", type=int, default=None,
                    help="FoldingConfig.num_sampling_steps, the diffusion ODE solver "
                         "steps. The API default is 100; higher is slower and meant "
                         "to be more accurate. Left unset, the API default stands.")
    ap.add_argument("--lm-dropout", type=float, default=None,
                    help="FoldingConfig.lm_dropout, dropout on the LM pair embeddings. "
                         "The API default is 0.3, i.e. dropout is ON at inference and "
                         "repeated runs differ. 0.0 turns it off.")
    ap.add_argument("--msa-depth", type=int, default=None,
                    help="FoldingConfig.msa_max_depth, the number of MSA rows randomly "
                         "subsampled each loop. The API default is 1024 out of our "
                         "14301. Raising it narrows the spread between runs, lowering "
                         "it widens them. Keeps subsampling ON, unlike --full-msa.")
    ap.add_argument("--full-msa", action="store_true",
                    help="set msa_max_depth=None, so the whole 14301-sequence MSA is "
                         "used every loop instead of a fresh random 1024-row subsample.")
    ap.add_argument("--out", help="write .cif and metrics here instead of this "
                                  "submission's raw/ and results/, so a variant run "
                                  "does not collide with the folded set")
    a = ap.parse_args()
    raw, res = (os.path.join(a.out, "raw"), os.path.join(a.out, "results")) \
        if a.out else (RAW, RES)
    os.makedirs(raw, exist_ok=True); os.makedirs(res, exist_ok=True)

    tok = os.environ.get("BIOHUB_TOKEN")
    if not tok:
        raise SystemExit("BIOHUB_TOKEN not set - run: source ~/.cyp_biohub_token")

    from esm.sdk import esmfold2_client
    from esm.sdk.api import FoldingConfig
    from esm.utils.structure import input_builder as ib

    client = esmfold2_client(model=a.model, url="https://biohub.ai", token=tok)
    kw = {}
    if a.lm_dropout is not None:
        kw["lm_dropout"] = a.lm_dropout
        print(f"lm_dropout = {a.lm_dropout} (API default is 0.3)", flush=True)
    if a.msa_depth is not None:
        kw["msa_max_depth"] = a.msa_depth
        print(f"msa_max_depth = {a.msa_depth} (API default is 1024)", flush=True)
    if a.full_msa:
        kw["msa_max_depth"] = None
        print("msa_max_depth = None, full MSA every loop (API default subsamples 1024)",
              flush=True)
    if a.sampling_steps is not None:
        kw["num_sampling_steps"] = a.sampling_steps
        print(f"num_sampling_steps = {a.sampling_steps} (API default is 100)",
              flush=True)
    cfg = FoldingConfig(include_pae=True, include_pair_chains_iptm=True, **kw)

    # ColabFold MSA (api.colabfold.com / MMseqs2), 14301 sequences, shared by all 20
    # ligands since the protein is identical. MMseqs2 writes trailing NUL bytes that
    # break strict parsers - AF3 died on them with "Unknown residues in MSA" - so the
    # file is cleaned into the scratch copy before loading.
    msa = None
    if not a.no_msa:
        from esm.utils.msa.msa import MSA
        src = os.path.join(ROOT, "submissions", "04_af3_default", "inputs",
                           "cyp3a4_colabfold.a3m")
        clean = os.path.join(RES, "cyp3a4_colabfold_clean.a3m")
        if not os.path.exists(clean):
            with open(src, errors="replace") as f, open(clean, "w") as g:
                for line in f:
                    g.write(line.replace("\x00", ""))
        msa = MSA.from_a3m(clean)
        print(f"MSA loaded: {len(msa.sequences) if hasattr(msa, 'sequences') else '?'} "
              f"sequences from {os.path.basename(src)}", flush=True)

    smiles = {r["Molecule_Name"].strip(): r["SMILES"].strip() for r in
              csv.DictReader(open(os.path.join(
                  ROOT, "input", "cyp-challenge-TEST-BLINDED_structures.csv")))
              if r["Molecule_Name"].strip()}
    names = [a.only] if a.only else sorted(smiles)

    rows = []
    for i, n in enumerate(names, 1):
        out = os.path.join(raw, f"{n}.cif")
        if os.path.exists(out):
            print(f"[{i}/{len(names)}] {n}: already folded")
            continue
        print(f"[{i}/{len(names)}] {n} ...", flush=True)
        t0 = time.time()
        inp = ib.StructurePredictionInput(sequences=[
            ib.ProteinInput(id="A", sequence=SEQ, msa=msa),
            ib.LigandInput(id="B", ccd=["HEM"]),
            ib.LigandInput(id="C", smiles=smiles[n]),
        ])
        r = client.fold_all_atom(inp, config=cfg)
        if isinstance(r, list):
            print(f"    NOTE: got a list of {len(r)} results; using the first")
            r = r[0]
        if hasattr(r, "error") or type(r).__name__ == "ESMProteinError":
            raise SystemExit(f"{n}: API error: {r}")
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

    if rows:
        p = os.path.join(res, "fold_metrics.csv")
        new = not os.path.exists(p)
        with open(p, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
