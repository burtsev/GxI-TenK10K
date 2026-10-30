#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Requires the original full-summary archive and official GRCh38/GENCODE44
# sources. The source-fetch command below verifies checksums and reuses cache.
# API key must be exported by the caller; never embed it here.
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock
.venv/bin/python -m pip install --no-deps -e .
.venv/bin/python scripts/fetch_study_sources.py --include-full
.venv/bin/pytest -q
outdir=results/followup-20261030
mkdir -p "$outdir/recovery"
# Historical loci are a frozen audited input, not recomputed from absent caches.
cp docs/runs/followup-20261030/preregistration/legacy_exclusions.tsv "$outdir/recovery/"
cp docs/runs/followup-20261030/preregistration/recovery_audit.json "$outdir/recovery/"
# Obtain the exact ontology release using its committed checksum, not a silently
# updated master file. Cached source retrieval/audit instructions are in README.
.venv/bin/gxi-tenk10k followup audit --out "$outdir"
.venv/bin/gxi-tenk10k followup freeze --out "$outdir"
.venv/bin/gxi-tenk10k followup score-dev --out "$outdir"
.venv/bin/gxi-tenk10k followup diagnostics --out "$outdir"
.venv/bin/gxi-tenk10k followup select-context --out "$outdir"
.venv/bin/gxi-tenk10k followup score-test --out "$outdir"
.venv/bin/gxi-tenk10k followup geometry --out "$outdir"
.venv/bin/gxi-tenk10k followup evaluate --out "$outdir"
.venv/bin/gxi-tenk10k followup report --out "$outdir"
