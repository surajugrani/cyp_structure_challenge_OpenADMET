#!/usr/bin/env python3
"""Block `sbatch` unless the script justifies its resource request.

Why this exists: resource padding kept recurring (40G for a 2G job, 16 CPUs for 10
concurrent tasks, 24h for a 3h run) despite a written rule to derive every number. The
rule failed because sbatch headers get written from a copied template before the
workload is thought about, and nobody revisits them. A hook cannot be forgotten.

The script must contain a line starting with `# RESOURCES:` that states where the
numbers come from. The hook also applies hard sanity limits and prints the project's
measured baselines so the derivation can be done on the spot.

Escape hatch: prefix the command with ALLOW_PADDED=1 for a deliberate exception.
"""
import json
import os
import re
import sys

MEM_CAP_GB = 24          # above this, justify explicitly in the RESOURCES line
CPU_CAP = 16
TIME_CAP_H = 12

BASELINES = """  measured baselines for this project (sacct MaxRSS / Elapsed):
    PrepWizard, 1 complex .............. <2 GB, ~4.5 min      (so 10 concurrent -> 16G)
    Glide prepwiz+grid+score, 4 conc. ... 7.5 GB total, 1h44m for 100 samples
    Glide mininplace only ............... ~1.5 GB per concurrent run
    ICM dock or in-place score .......... 16G ample; 100 samples ~1 h serial
    PoseBusters / RDKit / CSV work ...... 2-4 GB, 1 CPU
    Boltz-2 GPU, 20 ligands x 5 ......... 32G, 1 GPU, ~3 h
  rule: cpus-per-task = number of CONCURRENT processes, nothing more."""


def fail(msg):
    print(msg, file=sys.stderr)
    sys.exit(2)


def main():
    data = json.load(sys.stdin)
    cmd = (data.get("tool_input") or {}).get("command", "")
    if not re.search(r"\bsbatch\b", cmd) or "ALLOW_PADDED=1" in cmd:
        sys.exit(0)

    # Only inspect the script the scheduler is actually given. An earlier version
    # scanned the whole command for *.sh and picked up unrelated paths - e.g. a printf
    # writing another script's name into a file - and blocked on the wrong script.
    script = None
    for seg in re.split(r"(?:\|\||&&|\||;|\n)", cmd):
        m = re.search(r"\bsbatch\b(.*)", seg)
        if not m:
            continue
        for tok in re.findall(r"[\w./@-]+\.(?:sh|sbatch)", m.group(1)):
            p = tok if os.path.isabs(tok) else os.path.join(
                os.environ.get("CLAUDE_PROJECT_DIR", "."), tok)
            if os.path.exists(p):
                script = p
    if script is None:
        sys.exit(0)        # inline/heredoc submission - nothing to inspect

    text = open(script, errors="replace").read()
    has_just = re.search(r"^\s*#\s*RESOURCES:", text, re.M | re.I)

    def directive(name):
        m = re.search(rf"^#SBATCH\s+(?:--{name}|-{name[:1]})[=\s]+(\S+)", text, re.M)
        return m.group(1) if m else None

    mem = directive("mem")
    cpus = directive("cpus-per-task")
    tm = re.search(r"^#SBATCH\s+(?:-t|--time)[=\s]+(\S+)", text, re.M)

    problems = []
    if mem:
        g = float(re.sub(r"[^\d.]", "", mem)) / (1024 if mem.rstrip().upper().endswith("M") else 1)
        if g > MEM_CAP_GB:
            problems.append(f"--mem={mem} exceeds the {MEM_CAP_GB}G soft cap")
    if cpus and cpus.isdigit() and int(cpus) > CPU_CAP:
        problems.append(f"--cpus-per-task={cpus} exceeds the {CPU_CAP} soft cap")
    if tm:
        parts = tm.group(1).split("-")
        hh = int(parts[0]) * 24 if len(parts) > 1 else 0
        hh += int(parts[-1].split(":")[0])
        if hh > TIME_CAP_H:
            problems.append(f"--time={tm.group(1)} exceeds the {TIME_CAP_H}h soft cap")

    if has_just and not problems:
        sys.exit(0)

    msg = [f"BLOCKED: {os.path.relpath(script, os.environ.get('CLAUDE_PROJECT_DIR','.'))}"]
    if not has_just:
        msg.append("  missing a '# RESOURCES:' line justifying cpus / mem / time.")
    for p in problems:
        msg.append(f"  {p}")
    msg += [
        f"  requested: cpus={cpus} mem={mem} time={tm.group(1) if tm else None}",
        "",
        BASELINES,
        "",
        "  Add a line like:",
        "  # RESOURCES: 6 structures to prep, 6 concurrent -> 6 CPUs; PrepWizard <2G each",
        "  #            -> 12G; ~4.5 min each so ~10 min -> 1h.",
        "",
        "  Deliberate override: prefix the command with ALLOW_PADDED=1.",
    ]
    fail("\n".join(msg))


if __name__ == "__main__":
    main()
