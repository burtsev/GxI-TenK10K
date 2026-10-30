# Implementation, source contracts and reproduction

The final live scientific report is [docs/runs/live-20261029/report.md](runs/live-20261029/report.md). The earlier credentials-blocked report is retained as historical preflight evidence, not the final outcome.

## Setup and exact workflow

Python >=3.11 and samtools are required. Dependencies, including development tests and plotting, are pinned in requirements.lock.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install --no-deps -e .
.venv/bin/python -m pytest -q
.venv/bin/gxi-tenk10k audit --config configs/pilot.yaml --out results/pilot-20261029
.venv/bin/gxi-tenk10k select --config configs/pilot.yaml --out results/pilot-20261029
.venv/bin/python scripts/fetch_study_sources.py --include-full
.venv/bin/python scripts/audit_reference_smoke.py
# Supply GI_API_KEY in the command environment; never put it in Git.
.venv/bin/python scripts/smoke_sequence_diagnostics.py --flank 40960
.venv/bin/python scripts/smoke_sequence_diagnostics.py --flank 20480
.venv/bin/gxi-tenk10k run --verified-live --out results/live-20261029 --max-requests 760
.venv/bin/python scripts/compare_trial_cells.py
.venv/bin/gxi-tenk10k evaluate --out results/live-20261029
.venv/bin/gxi-tenk10k report --out results/live-20261029
```

The user's key was provided in an ignored local .env file. Execution exported that file into each command's process environment without printing it; the GI adapter reads GI_API_KEY only from environment. Credentials are absent from logs and tracked files.

`--include-full` is an explicit opt-in for the one required 16.6-GB common 100-kb ZIP. It does not fetch the other large releases. A verified complete source is reused; partial Zenodo downloads resume only with a matching HTTP Content-Range. Exact sizes, provided MD5 checksums and independent SHA256s are checked. Reference FASTA and GTF come from the official GENCODE v44 release and pass its MD5SUMS. FASTA expansion has a bounded size. ZIPs are streamed, with no filesystem extraction or traversal; targeted full-summary members are validated against expected names and size/ratio limits and read through CRC verification.

The first executed live workflow used `scripts/run_live_pilot.py` (now a thin wrapper for pipeline.live.run_live). A first pass had 100 candidates/76 matched controls. That pass is preserved under results/live-20261029/initial_unmatched_pass. Matching availability was subsequently moved before inference, using the same tolerances and source order, to construct a complete matched cohort without effect-based eligibility. The final report records this adjustment and the geometry amendment. The actual initial and subsequent commands are recorded in command_history.json. Re-running with existing sources and cache resumes without repeating successful predictions; evaluation and report never call GI.

## CLI and stage gates

All subcommands (`audit`, `select`, `smoke`, `infer`, `evaluate`, `report`, `run`) support `--help`. Default run executes conservative initial gates; `--verified-live` explicitly enables the audited full-summary workflow and its documented geometry/matching amendments. Exit 2 means blocked/incomplete scientific stages, 1 means command/data error, and 0 means the requested command completed; numerical success must be read from metrics, never the exit code.

The standard normalized-input path also accepts `--observations`, `--fasta`, `--reference-manifest`, and `--quota-preflight`. Observation TSV fields are gene, locus, cell, SNP, chrom, pos1, ref, alt, beta, effect_allele, effect_allele_source, PIP, source_file, genome_build, tss1, strand, maf, distance, gc, category. Reference manifest requires GRCh38, GENCODE v44 and exact FASTA SHA256. A standard quota preflight requires explicit cost approval and remaining-request evidence; the verified-live path instead stays below 800 requests and documents that monetary quota/pricing is unavailable. No quota balance is inferred from burst headers. A new normalized cohort needs a new output directory; the shared inference cache is reused.

The live workflow runs smoke before trial, trial before pilot expansion. Pilot inputs use one discovery context per pair, balanced across the four requested blood lineages; the first 20 candidates are additionally compared in all four available contexts. The cap is cumulative across cache-ledger attempts, including failures and retries. One initial HBB POST preceded that ledger and is counted separately. Serial inference stays below the discovered sustained rate; Retry-After, bounded retries and ambiguous network-failure stops are implemented. Only one process may write the same cache. No account balance/pricing endpoint is invented, and the report states actual attempt counts and remaining uncertainty about cost.

## Biological and data contracts

The two small archives use comma-delimited `.tsv` members. They contain gene-level p-values and SuSiE credible-set/PIP summaries, without beta or SE. The full ZIP is tab-delimited, with an initial row-index column and CHR, POS, MarkerID, Allele1, Allele2, AF_Allele2, BETA, SE, N and gene fields. Its exact observed schema is saved. The official SAIGE-QTL source documents Allele2 as ALT/effect allele and BETA as its effect; row identifiers and alleles are validated before paid scoring.

Coordinates: genomic pos1/TSS are one-based; FASTA slices and API indices/windows are zero-based half-open. Edits happen on genomic plus sequence, then minus-strand genes are reverse-complemented. Minus orientation shifts the fetched interval by one base to center the oriented TSS exactly; beta sign does not flip for gene strand. Only unique MANE_Select/Ensembl_canonical transcript TSSs from GENCODE v44 are accepted. Observed beta/SE must be finite and SE positive. REF must match FASTA. Missing PIP for controls is retained as missing, never zero.

Cell mapping: B naive -> B_naive; CD4 naive T -> CD4_Naive; CD14 classical monocyte -> CD14_Mono; CD16 nonclassical monocyte -> CD16_Mono. GI uses the human-readable names in the exact config context template. Biological cell description, assay and model are identical between alleles.

Prescreening uses per-cell BH q<=0.05 on gene ACAT, SuSiE PIP>=0.9, nonpalindromic biallelic SNVs, canonical annotation and nominal coverage, followed by deterministic seed/hash order and round-robin cells with distinct genes. Initial small-archive select results are provisional; the live workflow reads verified full observed betas and matches controls before inclusion. Controls share gene/cell and meet MAF difference<=0.05, TSS-distance difference<=5,000 bp and 201-bp reference GC difference<=0.1, with known PIP<0.5 or explicitly unreported PIP. They are uncertain comparators, not true negatives. No beta sign, delta, or model efficacy is used in matching.

## API contract and geometry amendment

The live OpenAPI/model catalog are archived and hashed. Supported GI integration-kit code was inspected; its predict method discards successful response headers, so this adapter uses documented httpx REST to retain rate evidence and persist budget before dispatch. It checks typed prediction.expression_log_tpm, unit, model, context, TSS, submitted length and scored-window provenance. The long model recommends 40,960 bp per flank but tokenization led to different effective endpoints for REF/ALT in the original 81,920-bp smoke. Strict equality rejected those comparisons. Before trial, a documented geometry amendment set 20,480-bp flanks (40,960 total), which the model read in full for both alleles. Those results are reported as amended-geometry measurements, with the original strict endpoint unresolved. A shorter input materially changed one smoke-effect sign; the report highlights this sensitivity.

Every prediction carries exact REF/ALT hashes, cache keys, context, requested model, source row/member, coordinates and reported scored interval. Keys also include API version/OpenAPI hash and base URL. The API does not expose an immutable weight revision; silent model deployment changes remain a reproducibility limitation. DNA stays out of committed tables.

## Statistics and exported evidence

Metrics distinguish all rows, missing estimates, near-zero abstentions and sign-eligible rows. Gene windows within +/-100 kb are merged into conservative locus clusters, shared across cells/controls. Bootstrap and null sign-flips preserve these clusters. Splits are disjoint by locus; no fitting uses test outcomes. Saved outputs include candidate/control metrics, per-cell results, trial results, PIP/distance/epsilon sensitivities, a distance-magnitude baseline, matched-control magnitude priority proxy, cell contrasts, plots and source/code/environment manifests. Missing/insignificant associations are not certain negatives. SAIGE count-model coefficients are not declared calibrated to log(TPM+1). Cross-cell contrasts are descriptive because donor covariance and power differ. External sequence-model baselines and causal fine-mapping validation were not performed.

Small sequence-free measured tables and aggregate outputs are committed under docs/runs/live-20261029. Large archives, reference FASTA, association pool, hundreds of response payloads and secrets remain ignored under data/ and results/. The historical initial report preserves the credentials-blocked stage; it does not describe the final live outcome.
