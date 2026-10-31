# Codex runbook: comprehensive, nonblocking TenK10K atlas

Read ATLAS_START.md, ATLAS_PROTOCOL.md, ATLAS_AVAILABILITY.md and ../configs/atlas.yaml. Implement the new study on the existing repository rather than rerunning the historical pilot or fixed-size followup.

**Scope of this commit:** instructions and proposed configuration. The commands below are implementation targets until Codex adds and tests them. Do not claim the new scheduler exists merely because these documents were committed.

## 0. Preserve existing work and inspect the implementation

Inspect git status, branch and diff before modifying anything. Fetch remote changes without a destructive reset, clean, forced checkout or blind overwrite of uncommitted VM artifacts. Preserve all historical results, frozen protocols, manifests and inference caches. Use a new results/atlas-<UTC timestamp>/ directory and a separate cache namespace when the API/model/precision/input contract differs.

Read the existing CLI, followup package, genome editing, GI client, analysis helpers and tests. The preceding followup already implements inference/cache/provenance utilities; extend or reuse validated components rather than unnecessarily rewriting them. Do not call the old followup freeze/select pipeline as the atlas entrypoint: its high-PIP, matched-control, sample-size, six-cell and all-or-nothing completion rules do not define this new study. Preserve those rules when reproducing historical runs.

Inspect live repository state instead of assuming the earlier missing-code problem remains. Use committed historical membership/registry metadata when raw VM caches are unavailable. Incomplete provenance disables untouched-holdout claims, not descriptive all-data scoring.

## 1. Implement and test skip-and-continue first

Implement an atlas-specific module runner with declared required and optional inputs, typed availability outcomes, durable status and validated output partitions. Parse and validate configs/atlas.yaml; explicit null optional paths must be accepted. Reject malformed configuration clearly, but a missing optional source is not a configuration error.

Resolve availability once per source/partition, using verified local artifacts before remote checks. Follow the bounded retries, licensing and status contract in ATLAS_AVAILABILITY.md. Never initiate access applications or wait for a response. Never print API keys, tokens or identifiable donor data. Record unmet requirements in skipped_analyses.md while proceeding.

Allow downstream modules to consume validated partial output; do not require every upstream module to have global status completed. Catch known missing-input cases at the smallest useful scope. Isolate unexpected exceptions, record failed_execution and continue unrelated modules. Quarantine failed integrity checks instead of silently labelling them missing. Final reporting/status writes must execute even when scientific modules skip or fail.

Add all 16 offline availability scenarios in ATLAS_AVAILABILITY.md as acceptance tests. Run these before large live inference. They must prove that missing expression does not stop common or rare summary analyses, missing PIP does not stop allele scoring, and a missing chromosome does not discard the remaining chromosomes. Test resume with newly supplied optional input and no replay of successful requests.

## 2. Audit data, sources and API independently

Reuse checksummed data/cache and the prior source manifest. Resolve the current Zenodo record, inventories and schemas; inventory common/rare, single-variant/gene-level and fine-mapping families, including overlapping cis windows. Check archive/member integrity, normalizable variant IDs, effect alleles, beta/SE definitions, cell IDs and transcript annotations. Record available and unavailable partitions; process available ones immediately.

Check approved aggregate or donor expression paths and authorization/normalization metadata. EGA access is not assumed; without an approved available source mark the expression branch skipped and continue. Do not start raw controlled-data ingestion merely to satisfy a protocol target. Optional annotations, covariance, LD, rare test membership and external scores have their own availability records.

Perform GI authentication/schema/model checks without exposing credentials. Owner reports no total GI quota: no historical request-count gate. Honor live rate/concurrency limits and bounded retry rules. If live inference is unavailable, continue auditing and evaluating verified matching cached predictions. Mark unscored coverage accurately; never invent model responses.

Use an explicit historical/development panel for raw precision, repeatability and geometry diagnostics. If raw model outputs/embeddings cannot be obtained, log the limitation and continue API-visible scoring. If base sequence geometry is invalid, prevent affected live biological scoring, but still run valid cached/data-only analyses and reporting.

## 3. Build the full eligible universe and freeze new choices

Construct a partitioned normalized ledger of every available variant-gene-cell-source observation. Include common and rare observations regardless of significance or PIP; do not require a matched control or expression target. Inventory genes outside the eQTL set for expression where observed data permit. Deduplicate requests and overlapping identical source records while retaining distinct source model fits.

For each row resolve exact coordinates, strand, alleles, TSS and model coverage. Use symbolic/out-of-window/ambiguous/unsupported statuses rather than dropping rows. Handle SNVs first; add indels only after geometry tests pass. A missing field that only supports a subgroup must not delete an otherwise valid score.

Freeze sources, canonical contexts, model/API/input/precision version, metric definitions and prespecified strata. Historical rows remain flagged in descriptive analysis. Calibration, aggregation thresholds and any model changes require development and held-out locus separation. If the information needed to prove independence is absent, label the analysis descriptive and proceed rather than declaring it confirmatory or aborting.

Create the unique request manifest before inference. Count source rows separately from unique REF/ALT sequence-context-model requests. Write intended counts and coverage, with no arbitrary cap or significance-based sample limit. Partition execution for manageable disk and memory. Do not estimate the whole workload by assuming every source row requires two new requests.

## 4. Score all available eligible inputs, checkpointed

Reference expressions are shared by exact gene/TSS/sequence/cell/model/context. ALT predictions are unique by actual edited sequence and context. Validate the genomic REF before editing, effect allele mapping, gene-sense orientation, TSS and identical genomic effective REF/ALT coverage. Do not reuse a cached rounded response as a new full-precision output.

Use bounded streaming and concurrency, durable per-request results, exact-key locking and exponential backoff/Retry-After. Record retries and ambiguous dispatches separately from successful outputs. No automatic 800/1600 total-call ceilings; rate limits remain mandatory. Do not claim a monetary price not supplied by the service.

A permanently failing request/cell/partition is recorded and skipped after bounded retries; valid other work continues. Persistent global auth/outage opens a circuit breaker for new live calls, not a loop across millions of requests. Evaluate caches and other modules and finalize a partial/unavailable report. A future explicit resume may retry. Disk/resource limits produce a partial ledger rather than silent data truncation or deletion of prior work.

## 5. Evaluate independent scientific modules

Implement each benchmark in ATLAS_PROTOCOL.md as an independently runnable function with explicit input schemas and outputs:

- Expression: donor-balanced reference when donor data exist; aggregate-only rank accuracy otherwise; per-cell across genes, per-gene across cells, gene-centered specificity, optional calibration/reliability. Missing observed expression skips accuracy, not reference predictions or variant metrics.
- Direction: all-data, discovery-cell, reliable-sign and uncertain-cell strata; correct/incorrect/abstaining/missing counts; proper zero-preserving nulls and coverage. Do not count a PIP-high discovery variant as proven active in every comparison cell.
- Magnitude: signed and absolute-effect Spearman, including zeros/ties; per-cell, per-gene-cell, pooled and equal-locus summaries. Missing SE only removes precision-dependent analyses.
- Cell specificity: all available cells per pair; profile and contrast correlations, sign/coverage, pair/locus weighting and cell-pair matrices. Missing covariance skips formal tests only.
- Fine-mapping: explicit known-PIP label discrimination with prevalence/AUPRC and AUROC; complete or labelled partial-set ranking; unknown PIP separate. Lack of a matched control never blocks core candidate evaluation.
- Rare variants: supported single-variant results, then compatible gene-set enrichment if inputs exist. Gene-level retesting is separately gated on sufficient approved genotype/score-covariance data.
- Subgroups: all prespecified available dimensions; skip unavailable expression/annotation-based dimensions individually. Small/constant groups are descriptive or not_estimable, not crashes.
- External comparators: use approved executable compatible data only, each independently; missing licence/track mapping/access skips that comparator. Do not purchase or initiate third-party paid workloads without authorization.

Use locus-cluster resampling and specified multiplicity. Keep utility gates distinct from p-values; no historical p=1 overwrite rule in the atlas. Preserve historical metrics unmodified. Export complete denominators, errors and skip reasons, not just the best-performing subsets.

## 6. Always report; commit only sanitized outputs

Generate report.md, availability.json, module_status.tsv/json, skipped_analyses.md, coverage.tsv, source/prediction manifests, metrics.json, environment.json and resume_manifest.json. Keep every benchmark heading even if its output is SKIPPED. State exactly which scientific questions were answered, partially answered, not estimable or skipped and why. An all-skipped run reports no scientific findings, not a successful benchmark.

A final summary must distinguish full intended release, available sources, model-scorable universe, successfully scored rows, pending/failed partitions and included metric denominators. Write plots only from actual outputs. Missing optional data must not prevent summary generation.

Run old tests for regression plus new atlas unit/integration tests; report actual test counts, not inherited followup counts. Commit code, config, tests and small authorized aggregates to a clearly identified branch/commit. Do not commit raw donor data, DNA, API credentials, large archives or entire response caches. Controlled aggregate outputs require authorization review before public publication.

## Proposed CLI to implement, not existing commands

```bash
# Expected interface after implementation. Validate help and offline tests first.
gxi-tenk10k atlas audit --config configs/atlas.yaml
gxi-tenk10k atlas plan --config configs/atlas.yaml
gxi-tenk10k atlas run --config configs/atlas.yaml
gxi-tenk10k atlas resume --config configs/atlas.yaml --run-dir results/atlas-<run-id>
gxi-tenk10k atlas evaluate --config configs/atlas.yaml --run-dir results/atlas-<run-id>
gxi-tenk10k atlas report --config configs/atlas.yaml --run-dir results/atlas-<run-id>
```

`run` orchestrates available modules, not an all-or-nothing shell chain. `evaluate` and `report` should work from compatible saved outputs without making fresh API calls. `resume` rechecks source availability once and handles only unfinished/newly enabled work; preserve previous snapshots.

## Completion criteria

The software implementation is acceptable when availability tests and regression tests pass, statuses/denominators are correct, inputs and caches are traceable, independent modules continue after unavailable inputs, and reporting completes. Scientific completion is a separate coverage statement: do not call a partial atlas a full-release evaluation. Leave executable resume commands and remaining-work manifests when anything is incomplete; never fabricate completion.
