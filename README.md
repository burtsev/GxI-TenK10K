# GxI × TenK10K — expression and regulatory-variant benchmarks

## Current study: comprehensive atlas (31 October 2026)

The new study evaluates **gene-expression accuracy**, **allelic effect direction**, **effect magnitude**, **cell-type specificity**, **fine-mapping**, and **rare-variant prioritization**, with performance by cell type, gene and other prespecified subsets. It targets the complete available, model-scorable dataset rather than another selected-variant pilot.

**Codex entrypoint: [docs/ATLAS_START.md](docs/ATLAS_START.md).** Read [AGENTS.md](AGENTS.md), the [scientific protocol](docs/ATLAS_PROTOCOL.md), [availability/skip contract](docs/ATLAS_AVAILABILITY.md), [execution runbook](docs/ATLAS_CODEX_RUNBOOK.md), and [atlas configuration](configs/atlas.yaml).

**Unavailable data must not stop independent work.** Without approved expression observations, skip expression accuracy and dependent strata but continue variant benchmarks. Missing rare-variant inputs, PIP/credible sets, donor genotypes, annotations, raw server precision or external-model access skip only the corresponding analyses. Missing cells/chromosomes produce explicitly partial coverage. Always write a report listing completed, partial, skipped, not-estimable and failed modules with reasons. Do not wait for optional access or fabricate substitutes. Invalid alleles/references are quarantined, not silently treated as zeros.

Use the owner's unrestricted-total-quota `GI_API_KEY` via direct GI API/official client, with rate-limit-aware requests, bounded retries and resumable caches. Historical sample limits, six-cell caps and total-request gates do not apply to the atlas. Preserve previous code, frozen outputs and VM data.

**Status:** the atlas protocol, configuration and Codex implementation instructions are committed. The atlas-specific scheduler/CLI and acceptance tests still need implementation; this update does not claim a newly executed atlas or new results. Sixteen missing-data/error/resume acceptance scenarios are specified in the availability contract. Earlier studies and their measured results remain below for reproduction.

Suggested Codex task:

> Read docs/ATLAS_START.md and implement and execute docs/ATLAS_CODEX_RUNBOOK.md using configs/atlas.yaml. Run every available module of the whole-data study; skip only analyses whose required inputs are unavailable, without waiting or stopping independent work. Preserve previous outputs, add skip/resume tests, and always generate an explicit coverage/status report. Use my GI_API_KEY with no total quota, retaining server rate limits and caching. Commit code, tests and small authorized summaries; never fabricate missing data or results.

## Historical pilot: original scope and setup

This repository is the execution brief for **Codex running on a Linux VM**. The original goal was to rigorously test whether Genomic Intelligence (GI) predicts naturally occurring variant effects on gene expression in TenK10K immune cell types, using the [TenK10K Phase 1 association statistics](https://zenodo.org/records/18221260) (DOI: [10.5281/zenodo.18221260](https://doi.org/10.5281/zenodo.18221260)).

For historical pilot reproduction, follow [docs/CODEX_RUNBOOK.md](docs/CODEX_RUNBOOK.md) and [docs/BENCHMARK_PROTOCOL.md](docs/BENCHMARK_PROTOCOL.md). These are not the entrypoint or completion gates for the current atlas. Do not misrepresent instructions or smoke tests as benchmark results.

## Original pilot questions

1. **Primary:** sign concordance between GI-predicted `ALT − REF` log(TPM+1) and observed allele-oriented TenK10K eQTL effects; stratify by cell type.
2. **Secondary:** association of predicted effects with observed beta magnitude/ranks and cross-cell-type effect contrasts.
3. **Exploratory:** discrimination of credible-set variants versus well-matched comparison variants; PIP-aware ranking, with LD caveats.
4. **Coverage:** report how many published associations are physically inside the model's *effective scored window*, not merely inside its submitted 81,920-bp sequence.
5. **Baselines:** distance-to-TSS and matched-frequency controls; two GI expression models if budget permits. External sequence-model comparisons only if executable, and label omissions honestly.

## Source of truth

- Zenodo record: https://zenodo.org/records/18221260 ; inspect live Zenodo metadata for record version, date and checksum. It contains **summary associations, not donor-level genotype/individual expression data**.
- Locate and verify the original paper from the Zenodo metadata before citing it.
- TenK10K analysis code: https://github.com/powellgenomicslab/tenk10k_phase1 and https://github.com/populationgenomics/saige-tenk10k
- GI live API: https://docs.genomicintelligence.ai/rest-api ; machine-readable contract: https://api.genomicintelligence.ai/v1/openapi.json ; models: https://api.genomicintelligence.ai/v1/tasks/expression/models
- Reference: GRCh38; verify exact TenK10K gene annotation and effect-allele conventions against actual files/code before evaluating.

## VM and access: historical pilot sizing

Recommended for the earlier pilot: Ubuntu 24.04, **16 vCPU / 64 GiB RAM, 500 GB SSD, no GPU**. For only a small pilot, 8 vCPU / 32 GiB is sufficient; for the entire 170.6-GB compressed record, plan at least 2 TB disk and substantially more CPU/RAM. GI inference uses its remote API, not the VM's GPU. For the atlas, audit actual compressed/expanded partitions, cache size and disk before execution rather than assuming pilot sizing covers the complete release.

Install Codex and Python 3.11+; provide the API key **as an environment variable**, not in Git. Reuse the existing VM checkout when it holds study artifacts; do not overwrite it with a blind clone/reset.

```bash
# For a new checkout only:
git clone https://github.com/burtsev/GxI-TenK10K.git
cd GxI-TenK10K
# Supply GI_API_KEY securely through your environment/secret manager.
export GI_BASE_URL='https://api.genomicintelligence.ai'
```

See [docs/VM_SETUP.md](docs/VM_SETUP.md) for historical infrastructure setup and [docs/ATLAS_START.md](docs/ATLAS_START.md) for the current task.

## Outputs

Study outputs live under `results/<run_id>/` with source manifests/checksums, predictions, exclusions, evaluation metrics and reports. The new atlas additionally requires per-source availability, module statuses, full coverage denominators and an explicit skipped-analyses report. Large archives, raw donor data and secrets must never be committed to Git.

## Scientific scope

This is an observational benchmark, not an intervention, a clinical model validation, or a claim of causal variant identification. An insignificant eQTL is not automatically a true negative, SuSiE PIP is not experimental causality, and published betas may not share units with GI's log(TPM+1). Distinguish frozen no-fit inference from target-informed calibration/selection and verified experiment holdouts from unknown model-training overlap.

## Implemented benchmark and executed results

See [setup and reproduction](docs/IMPLEMENTATION.md) and the
[final live scientific report](docs/runs/live-20261029/report.md).
The live run completed a 20-pair trial and 100-candidate/100-control matched
pilot across four discovery cell types, with supplementary cross-cell queries.
The report records input-geometry and matching adjustments, abstentions,
uncertainty, source checksums and limitations. Original-protocol and amended
results are distinguished. The earlier credentials-blocked audit remains as
historical evidence.

## Follow-up study (30 October 2026): protocol archive

The initial 100-candidate exploratory pilot, its weak cell-specific contrast results, and subsequent prompt-conditioning experiments motivated a larger independent evaluation. Its specifications are [Follow-up Codex instructions](docs/FOLLOWUP_START.md), [preregistered protocol](docs/FOLLOWUP_PROTOCOL.md), [execution runbook](docs/FOLLOWUP_CODEX_RUNBOOK.md), and [follow-up configuration](configs/followup.yaml). Prior findings are summarized with caveats in [prior findings](docs/FOLLOWUP_PRIOR_FINDINGS.md).

The follow-up used the GI API through `GI_API_KEY` with no total quota cap for this key, while respecting network/server rate limits. Old pilot targets and `api.max_requests: 1600` belong only to the historical pilot. Its execution results are below; neither this protocol nor the pilot rules supersede the current atlas instructions.

## Completed independent follow-up (30 October 2026)

The frozen test includes 1,000 variant–gene pairs across 877 genes and 777 physical locus clusters, with 1,000 matched uncertain controls and 28 immune populations. None of the gene-free alternatives qualified on 250 development loci, so baseline alone was tested. All-pair direction was 40.8% for variants (95% locus CI 38.6–43.0%) and 42.9% for cell contrasts (41.4–44.4%); weak positive rank associations were detected, but neither endpoint crossed the published 50% utility gate. All 68 offline tests and 27 final QA checks passed.

Read the [five-page manuscript](docs/runs/followup-20261030/GI_TenK10K_followup_manuscript.pdf), [full summary](docs/runs/followup-20261030/summary.md), and [reporting/reproduction notes](docs/FOLLOWUP_REPORTING_NOTES.md). Frozen inputs and the pre-test decision are committed alongside the sanitized aggregates. Raw data and inference caches remain on the VM.
