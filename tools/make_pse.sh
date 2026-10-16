#!/bin/bash
# Build the PyMOL session for one submission. Safe to call at the end of any
# validate/packaging script: it is a no-op if no zip was written.
#
#     bash tools/make_pse.sh submissions/16_boltz2_pliptm
#
# Runs in-process (PyMOL is cheap here, ~1 min) rather than submitting a second job, so
# the session always exists by the time the packaging job finishes.
set -uo pipefail
SUB="${1:?usage: make_pse.sh <submission dir>}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if ! compgen -G "$SUB/*.zip" >/dev/null; then
  echo "make_pse: no zip in $SUB, skipping session"
  exit 0
fi
module load pymol/3.1.0 2>/dev/null || true
CYP_ROOT="$ROOT" CYP_SUBMISSION="$SUB" pymol -cq tools/make_pymol_session.py
