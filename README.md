# GxI × TenK10K — zero-shot cell-type-specific eQTL benchmark

This repository is the execution brief for **Codex running on a Linux VM**. The goal is to rigorously test whether Genomic Intelligence (GI) predicts naturally occurring variant effects on gene expression in TenK10K immune cell types, using the [TenK10K Phase 1 association statistics](https://zenodo.org/records/18221260) (DOI: [10.5281/zenodo.18221260](https://doi.org/10.5281/zenodo.18221260)).

**Codex should start with [AGENTS.md](AGENTS.md), then follow [docs/CODEX_RUNBOOK.md](docs/CODEX_RUNBOOK.md) and [docs/BENCHMARK_PROTOCOL.md](docs/BENCHMARK_PROTOCOL.md).** This repository initially contains specifications and configuration; Codex must implement, test, run, and document the pipeline. Do not misrepresent instructions or smoke tests as benchmark results.

## What to test

1. **Primary:** sign concordance between GI-predicted `ALT − REF` log(TPM+1) and observed allele-oriented TenK10K eQTL effects; stratify by cell type.
2. **Secondary:** association of predicted effects with observed beta magnitude/ranks and cross-cell-type effect contrasts.
3. **Exploratory:** discrimination of credible-set variants versus well-matched comparison variants; PIP-aware ranking, with LD caveats.
4. **Coverage:** report how many published associations are physically inside the model's *effective scored window*, not merely inside its submitted 81,920-bp sequence.
5. **Baselines:** distance-to-TSS and matched-frequency controls; two GI expression models if budget permits. External sequence-model comparisons only if executable, and label omissions honestly.

## Source of truth

- Zenodo record: https://zenodo.org/records/18221260 (version v2, published 2026-01-30); contains **summary associations, not donor-level genotype/individual expression data**.
- Paper: https://www.nature.com/articles/s41586-026-10721-2
- TenK10K analysis code: https://github.com/powellgenomicslab/tenk10k_phase1 and https://github.com/populationgenomics/saige-tenk10k
- GI live API: https://docs.genomicintelligence.ai/rest-api ; machine-readable contract: https://api.genomicintelligence.ai/v1/openapi.json ; models: https://api.genomicintelligence.ai/v1/tasks/expression/models
- Reference: GRCh38; verify exact TenK10K gene annotation and effect-allele conventions against actual files/code before evaluating.

## VM and access

Recommended: Ubuntu 24.04, **16 vCPU / 64 GiB RAM, 500 GB SSD, no GPU**. For only a small pilot, 8 vCPU / 32 GiB is sufficient; for the entire 170.6-GB compressed record, plan at least 2 TB disk and substantially more CPU/RAM. GI inference uses its remote API, not the VM's GPU.

Install Codex and Python 3.11+; provide the API key **as an environment variable**, not in Git:

```bash
git clone https://github.com/burtsev/GxI-TenK10K.git
cd GxI-TenK10K
export GI_API_KEY='...'
export GI_BASE_URL='https://api.genomicintelligence.ai'
# Start Codex and give it the instruction below:
```

> Read AGENTS.md and docs/CODEX_RUNBOOK.md. Implement the complete reproducible TenK10K eQTL benchmark per docs/BENCHMARK_PROTOCOL.md and configs/pilot.yaml. Execute the preflight, data audit, smoke tests, budgeted 20-variant trial, and (after checking actual quota/data/schema) the 100-positive/100-control pilot. Persist results, plots, exclusions, and a limitations-aware scientific report. Do not stop at a plan or fabricated outputs.

See [docs/VM_SETUP.md](docs/VM_SETUP.md) for infrastructure setup, environment variables, storage layout and expected commands. Actual CLI commands in that file are **implementation targets**, not pre-existing working scripts.

## Outputs

Codex should produce `results/<run_id>/` with a raw-data audit, manifest/checksums, deterministic candidate tables, a request/response cache, exclusion reasons, predictions, evaluation metrics with uncertainty intervals, figures, and `report.md`. Large archives and secrets must never be committed to Git. See the acceptance checklist in the runbook.

## Scientific scope

This is a **zero-shot observational benchmark**, not an intervention, a clinical model validation, or a claim of causal variant identification. An insignificant eQTL is not automatically a true negative, SuSiE PIP is not experimental causality, and published betas may not share units with GI's log(TPM+1).
