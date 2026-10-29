# Codex execution runbook

## Stage 0: environment and discovery
1. Inspect repo and choose uv or pip with pinned lock; create `pyproject.toml`, CLI, package and tests. Install Python 3.11+, git, unzip, curl, build essentials, samtools/tabix and enough disk. Configure `GI_API_KEY`, optional `GI_BASE_URL`; never echo secrets.
2. Discover REST API from https://docs.genomicintelligence.ai/rest-api and https://api.genomicintelligence.ai/v1/openapi.json. Run authenticated model listing and minimal valid positive control, verifying response unit, sequence length, TSS index, scored-window metadata, rate-limit headers and likely costs. If network or API unavailable, implement mock client plus tests and explicitly report blockage.
3. Check disk `df -h`, free memory, CPU, HTTPS access, `ulimit -n`; enable log rotation, process resumability, timeouts. Record environment to results.
4. Query Zenodo metadata API, save record JSON and per-file manifest/checksums, display total compressed sizes. Identify source archive and permitted licence.

## Stage 1: data audit
Download only small ZIPs first, validate checksums, unzip list (without full extract), read representative table rows. Detect delimiters, filenames, compression, columns, chromosome/build, variant ID convention, beta scale, effect allele, gene IDs, cell types, sample sizes, PIP, SE, q-values. Map cell names explicitly and identify exact published source of reported stats. Write `data_audit.md`. Test candidate extraction from the small files. Only download the 16.6-GB full cis archive if necessary and disk permits; can stream members to DuckDB or Parquet by chunks.

## Stage 2: implementation
Build modular package:
`data/zenodo.py` manifest/checksum/resumable fetch (with safe archive listing);
`data/parse.py` autodetect schemas with explicit mappings/validation;
`genome/reference.py` GRCh38 FASTA + strand-aware TSS windows;
`genome/variants.py` SNV editor plus allele normalization;
`gi/client.py` SDK/REST adapter with retries, quotas, logging, response schema validation;
`pipeline/select.py` candidate/control selection;
`pipeline/score.py` batch/ref reuse, checkpoint, delta, coverage;
`analysis/metrics.py` uncertainty, CI, baselines, stratification;
`report.py` scientific result generator;
`cli.py` entrypoints. Code should use robust typing, concise docs, reliable exceptions, tests and deterministic seed.

## Stage 3: tests and smoke
Create fixtures for synthetic dataset, plus-strand and minus-strand SNVs and REF mismatch. Unit tests for one-based POS to zero-based interval, TSS-center, reverse complements, preserved length, cache key, identical context ref/alt, handling API 429/5xx and retry-after, malformed records, out-of-window exclusions, exact sign orientation and per-locus bootstrap. Add offline integration test with mock GI adapter. Then use 2 real variants, verify publication info and effect alleles, generate saved ref/alt predictions. Stop and repair any critical mismatch.

## Stage 4: budgeted trial
Select 20 distinct candidate SNVs across >=4 cells, where available, preferably high PIP and model coverage. Score REF/ALT, retain responses, validate distribution, evaluate coverage and sign orientation, compare to observed data. Check call count and quota: 20 pairs × 4 cell types yields at most ~160 requests before reference reuse; repeated genes lower this. If bad results, investigate input orientation and context; don't tune thresholds on primary test outcomes.

## Stage 5: full pilot
Scale to 100 candidate and 100 matched control variant–gene pairs, up to 4 cell contexts, zero-shot. The **worst-case** without reusable reference calls is 1,600 API calls; actual ref dedup reduces this. Implement `--dry-run` and `--max-requests`, queue bounded concurrency (e.g. 2–4, only after checking API limits), retry backoff and rate limit. Save checkpoint after every inference. Generate metrics, 95% clustered CI, data selection flow, plots (predicted vs observed effect, per cell type, reference vs alt distribution, eQTL coverage, permutation/null, fine-mapping if present), reproducibility manifest, and full narrative.

## Stage 6: interpretation and reporting
Report (1) sample counts at every filtering stage; (2) central result and 95% CI; (3) direction vs magnitude vs cell specificity; (4) error/coverage analysis; (5) controls; (6) limitations including observational eQTLs, LD, and transcript/assay differences; (7) next steps. Include statements only based on executed data. Write `results/<run_id>/report.md`. Keep every API error and unscored row in exclusions; leave no silent drops.

## Acceptance criteria
- `python -m pytest -q` passes without network using synthetic fixtures.
- `gxi-tenk10k --help` and `gxi-tenk10k audit/select/smoke/infer/evaluate/report/run --help` work.
- Selection deterministic by seed and config hash.
- Every table has documented columns, coordinate/allele orientation, source file, model/version and reference genome release.
- Every inference traceable by hash to REF and ALT, exact same cell context and scoring window, with resumable cache.
- Observations not misrepresented as causal evidence or true negatives.
- Summary metrics reproducible from saved predictions and observations without calling GI again.
- Report explains precisely what ran and what did not, including any quota, data licensing or access restrictions.
