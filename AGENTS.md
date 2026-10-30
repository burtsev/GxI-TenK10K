> **New study (30 October 2026):** If the task asks for the GI–TenK10K follow-up, follow [docs/FOLLOWUP_START.md](docs/FOLLOWUP_START.md), [docs/FOLLOWUP_PROTOCOL.md](docs/FOLLOWUP_PROTOCOL.md), [docs/FOLLOWUP_CODEX_RUNBOOK.md](docs/FOLLOWUP_CODEX_RUNBOOK.md), and [configs/followup.yaml](configs/followup.yaml). These override the original pilot's call limit and cohort size. GI_API_KEY is stated to have no total quota; keep rate limits, retries and checkpointing. Recover historical VM results before scoring an untouched holdout. Do not reset local uncommitted work.

# Instructions for Codex (read before making changes)

You are the implementing research engineer. Execute an auditable, reproducible **zero-shot** assessment of Genomic Intelligence (GI) DNA-to-expression predictions against TenK10K immune-cell cis-eQTL summary statistics. Read README.md, docs/BENCHMARK_PROTOCOL.md, docs/CODEX_RUNBOOK.md, configs/pilot.yaml and docs/VM_SETUP.md in that order. Implement working code and tests, not just a proposed plan.

## Operating requirements
- First inspect *actual* Zenodo archive inventory, CSV/TSV/Parquet schemas, effect allele/beta conventions, GENCODE release, cell naming, missing data, and fine-mapping fields. Save audit. Never assume filenames within ZIP or hardcode columns without detecting/validating them.
- Obtain data ONLY from Zenodo record 18221260 and documented official project sources; check DOI version, file SHA256 and provided checksums; include reproducible source manifest. Do not download all 170+ GB by default. Stream extraction and avoid zip bombs. No donor-level individual data are required.
- Use GI's *current* API documentation/OpenAPI and discover live expression models, context requirements, response formats and quotas. Prefer the supported GI SDK/client if documented; otherwise use httpx with documented REST endpoints. Do not invent request schemas. GI_API_KEY from environment only. Test a valid positive control before sampling.
- A FASTA reference must match GRCh38 and gene TSS, strand, sequence coordinates and reference alleles. Handle minus-strand genes with reverse-complemented input; normalize variants in genomic plus coordinates **before** strand transformation; explicitly check REF against FASTA, including indels (pilot SNVs only).
- Preserve exact model input geometry per discovered model; find TSS index/position in *oriented* sequence. An effective scored window smaller than input must be logged; exclude unsupported distal variants. Restrict primary to variants covered by GI's effective scoring interval.
- Do not infer donor-level expression from summary-level releases. Never treat non-significant associations as definite negatives. Avoid circular validation and leakage. Label secondary exploratory results as such.
- Cache inference using model ID, exact REF/ALT sequences hashes, cell-context string, API version, and hyperparameters; obey retry-after, concurrency, rate limits and API budget; no unbounded retries. Do not leak DNA sequences/API keys in logs unnecessarily.
- Use deterministic random seed; record source hashes, exact commands, Python package versions, model IDs, run timestamp and Git commit. Persist checkpoints so a cancelled study resumes without paying for prior predictions.
- Before 800+ calls, perform API rate-limit/quota and cost preflight. If blocked by API credentials/quota/network, do not fabricate findings; finish executable software/tests using mocks/fixtures and document blocker.
- Work sequentially through stage gates: preflight, audit, unit tests, 2-variant smoke, 20-variant trial, 100-positive/100-control pilot. Continue as far as runtime/quotas permit. Spend effort on correct inference before adding plots.
- Never commit raw archives, model secrets, thousands of inference JSON payloads, or huge result files. Commit a meaningful report and small aggregate outputs if permitted; otherwise provide output locations and counts.

## Deliverables
Python package `src/gxi_tenk10k/` (download, schemas, annotations, selection, sequences, gxi_client, inference, analysis, report); `tests/`; pinned dependency metadata; a CLI with subcommands `audit`, `select`, `smoke`, `infer`, `evaluate`, `report`, `run`; scripts as necessary; structured logs and a full result report.

## Scientific QA
Mandatory: testing SNV and minus-strand editing; REF mismatch failure; TSS coordinate / off-by-one; model length and scored-window coverage; cache determinism/resume; same allele orientation in observed and predicted effects; data split by locus/gene; zero-missing-data denominators and complete exclusions; bootstrap intervals clustered by gene/locus; calibration sensitivity. Use no target outcome to select inputs beyond prespecified eQTL significance/PIP criteria, and report ascertainment.

## Final communication
Report actual completed stages, sample sizes, pass/fail counts, metrics with uncertainty intervals, exact dataset/model versions, known limitations, commands to reproduce and residual blockers. Avoid describing planned work as done.
