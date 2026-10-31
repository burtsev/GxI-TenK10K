# Start the comprehensive GI x TenK10K study

This is the current task entrypoint for the **whole-data expression and variant-effect atlas**, distinct from the completed pilot and 1,000-pair followup. Read:

1. [Scientific protocol](ATLAS_PROTOCOL.md).
2. [Mandatory availability and skip contract](ATLAS_AVAILABILITY.md).
3. [Codex execution runbook](ATLAS_CODEX_RUNBOOK.md).
4. [Configuration](../configs/atlas.yaml).

These instructions/configuration are committed specifications. Codex must implement the atlas-specific CLI, scheduler and tests before executing them; no new atlas results are claimed here.

## Task for Codex

> Implement and execute the comprehensive GI x TenK10K study in docs/ATLAS_PROTOCOL.md using configs/atlas.yaml and docs/ATLAS_CODEX_RUNBOOK.md. Preserve existing work, historical results and validated caches. Evaluate gene expression when approved observations are available; score the complete available model-scorable common/rare variant universe across every available cell type; analyze direction, magnitude, cell specificity, fine-mapping, rare-variant prioritization and cell/gene/evidence subsets. Use my GI_API_KEY directly; it has no total quota. Keep server rate limits, bounded retries and durable resume. Crucially, unavailable data or capabilities must skip only the smallest dependent analysis: do not pause, ask me for optional inputs, or stop independent work. Follow docs/ATLAS_AVAILABILITY.md, implement its offline acceptance tests, and always generate a report with completed, partial, skipped, not-estimable and failed sections. Never invent missing observations, PIP, genotypes, raw precision or metrics. Commit implementation, tests, configuration and small authorized summaries; do not commit credentials, raw human data or DNA.

## On the study VM

Use the existing working directory so its data and caches can be reused. Before any integration of remote changes:

```bash
cd /path/to/GxI-TenK10K
git status
git fetch origin
# Review/reconcile local work with origin/main. Do not reset or clean the checkout.
# GI_API_KEY should already be provided through your secure environment.
```

Do not place the key in a prompt, tracked file or command example containing its literal value. No ChatGPT skill/MCP installation is required for this API-key workflow. Start Codex in this directory and provide the task above.

## Expected behavior when data are absent

No approved expression/count matrix: skip expression accuracy and dependent strata; run variant benchmarks. No rare files: skip rare modules only. No PIP/credible sets: retain core variant scoring and skip unsupported fine-mapping metrics. No donor genotypes: skip personalized expression and association retesting. No full-precision/embedding access or external comparator licence: record those limitations and proceed. One missing cell/chromosome does not discard remaining partitions.

Do not treat a corrupted reference or unresolved allele as valid data: quarantine that record/partition and continue unaffected work, reporting the failure. Missing inputs, integrity errors, undefined metrics and scientific null results are different states.

Every run must finish with an explicit availability/status report, even when no scientific metric is possible. Verified completed work is reused on resume; newly supplied optional data enable only the dependent unfinished modules.
