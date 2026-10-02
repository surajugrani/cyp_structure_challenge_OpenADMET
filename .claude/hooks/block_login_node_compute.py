#!/usr/bin/env python3
"""PreToolUse(Bash) guard: refuse heavy compute on the HPC login node.

Repeatedly during this project, Schrodinger/ICM/Boltz/PoseBusters work was started
directly on the login node instead of through SLURM. This hook makes that structurally
impossible rather than relying on remembering.

A Bash call is DENIED when a command segment invokes one of HEAVY_TOOLS, unless any of:
  * the command goes through sbatch or srun
  * SLURM_JOB_ID is set (we are already inside an allocation)
  * the command sets ALLOW_LOGIN_NODE=1 (deliberate, documented override)

Matching is on the COMMAND POSITION of each `;`/`&&`/`||`/`|`-separated segment, after
stripping leading VAR=value assignments and any directory prefix. So `grep glide out.csv`
and `ls *boltz*` are fine, while `glide -WAIT x.in` and `$SCHRODINGER/run python3 s.py`
are refused.

Exit 2 + stderr = deny, per the PreToolUse hook contract.
"""

import json
import os
import re
import sys

HEAVY_TOOLS = {
    "prepwizard", "glide", "ligprep", "generate_glide_grids",
    "icm64", "_dockScan", "_dockBatch", "boltz", "posebusters", "bust",
}
ESCAPES = ("sbatch", "srun")

# Interpreters that are heavy by virtue of the environment they come from: running
# `~/micromamba/envs/cyp_posebusters_env/bin/python` IS a PoseBusters run even though the
# binary is called "python". This is the case that slipped through on 2026-09-30.
HEAVY_ENVS = ("cyp_posebusters_env", "boltz")


def strip_heredocs(cmd):
    """Remove heredoc bodies before scanning.

    A `cat > f <<EOF ... EOF` body is DATA, not commands. Without this, writing
    documentation that mentions a tool name at the start of a line -- e.g. a markdown
    code block showing a generate_glide_grids invocation -- is misread as running it.
    """
    out = []
    lines = cmd.split("\n")
    i = 0
    while i < len(lines):
        out.append(lines[i])
        m = re.search(r"""<<-?\s*['"]?([A-Za-z_][A-Za-z0-9_]*)['"]?""", lines[i])
        if m:
            term = m.group(1)
            i += 1
            while i < len(lines) and lines[i].strip() != term:
                i += 1
        i += 1
    return "\n".join(out)


def segments(cmd):
    cmd = strip_heredocs(cmd)
    return [s for s in re.split(r"\|\||&&|[;|&\n]", cmd) if s.strip()]


def first_token(seg):
    toks = seg.strip().split()
    for t in toks:
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", t):   # env assignment prefix
            continue
        return t
    return ""


def offending(cmd):
    for seg in segments(cmd):
        tok = first_token(seg)
        if not tok:
            continue
        base = tok.rsplit("/", 1)[-1]
        if base in HEAVY_TOOLS:
            return base
        # `$SCHRODINGER/run python3 ...` and friends
        if base == "run" and "SCHRODINGER" in tok:
            return "$SCHRODINGER/run"
        # interpreter taken from a heavy micromamba env
        for env in HEAVY_ENVS:
            if f"/envs/{env}/" in tok:
                return f"{env} interpreter"
    return None


def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        return 0                      # never block on a malformed payload
    cmd = (data.get("tool_input") or {}).get("command", "")
    if not cmd:
        return 0

    if os.environ.get("SLURM_JOB_ID"):
        return 0
    if re.search(r"\bALLOW_LOGIN_NODE=1\b", cmd):
        return 0
    if any(re.search(rf"(^|[\s;|&(]){e}\b", cmd) for e in ESCAPES):
        return 0

    tool = offending(cmd)
    if tool:
        print(
            f"BLOCKED: '{tool}' would run on the login node.\n"
            f"Submit it with sbatch (or wrap in srun) instead.\n"
            f"Deliberate override: prefix the command with ALLOW_LOGIN_NODE=1.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
