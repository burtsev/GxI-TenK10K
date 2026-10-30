# GI × TenK10K zero-shot benchmark: executed status report

## Central finding

No live GI variant-effect measurements were obtained. H1 (direction), H2 (cell specificity), and H3 (PIP prioritization) remain untested. Sign concordance, magnitude correlation, confidence intervals and permutation p-values are unavailable, not zero. No inference about GI performance or causality is supported.

## Execution

Timestamp: 2026-10-29T22:31:11.919911+00:00. Source revision: 56951d77bb521dc8ec88bf5b2434ce814e271398.
Configuration SHA256: 5650ec2f1c49104a0ac073f069cc20a482645c9ff7e97822434147ed2f261087.

Stage status:

```json
{
  "audit": {
    "status": "complete"
  },
  "selection": {
    "status": "provisional_only"
  },
  "smoke": {
    "status": "blocked",
    "reasons": [
      "GI_API_KEY absent",
      "observations not supplied",
      "fasta not supplied",
      "reference_manifest not supplied",
      "quota_preflight not supplied"
    ],
    "n_scored": 0,
    "requests": 0
  },
  "trial": {
    "status": "blocked",
    "reasons": [
      "live smoke gate not passed",
      "GI_API_KEY absent",
      "observations not supplied",
      "fasta not supplied",
      "reference_manifest not supplied",
      "quota_preflight not supplied"
    ],
    "n_scored": 0,
    "requests": 0
  },
  "pilot": {
    "status": "blocked",
    "reasons": [
      "live smoke gate not passed",
      "GI_API_KEY absent",
      "observations not supplied",
      "fasta not supplied",
      "reference_manifest not supplied",
      "quota_preflight not supplied"
    ],
    "n_scored": 0,
    "requests": 0
  },
  "infer": {
    "status": "blocked",
    "reasons": [
      "live smoke gate not passed",
      "GI_API_KEY absent",
      "observations not supplied",
      "fasta not supplied",
      "reference_manifest not supplied",
      "quota_preflight not supplied"
    ],
    "n_scored": 0,
    "requests": 0
  },
  "offline_tests": {
    "status": "complete",
    "passed": 46,
    "failed": 0,
    "evidence": "test_results.txt"
  }
}
```

Selection flow:

```json
{
  "source_susie_rows": 4193682,
  "requested_cell_rows": 1012119,
  "gene_fdr_05_rows": 1012057,
  "pip_threshold": 0.9,
  "screened_pairs": 100,
  "screened_genes": 100,
  "screened_by_cell": {
    "B_naive": 25,
    "CD14_Mono": 25,
    "CD16_Mono": 25,
    "CD4_Naive": 25
  },
  "eligible_candidates": 0,
  "eligible_controls": 0,
  "reference_checked": 0,
  "excluded_screened_pairs": 100,
  "config_hash": "5650ec2f1c49104a0ac073f069cc20a482645c9ff7e97822434147ed2f261087",
  "note": "Screening only; all screened rows fail downstream eligibility. Source filtering counts are aggregate, not eligible efficacy denominators."
}
```

Metrics (only measured data are eligible):

```json
{
  "n_total": 0,
  "n_missing": 0,
  "n_abstain": 0,
  "n_scored": 0,
  "n_loci": 0,
  "n_genes": 0,
  "n_cells": 0,
  "sign_concordance": null,
  "ci95": null,
  "spearman": null,
  "permutation_p": null,
  "epsilon": 1e-06
}
```

Implementation scope is incomplete: automatic normalization of the unaudited full association archive, cross-cell differential analysis, fine-mapping ranking and scientific plot generation remain unimplemented. This report documents an executed preflight/data-audit/offline-validation run, not a completed efficacy benchmark. See docs/IMPLEMENTATION.md for these limits and the verified-input interface.

The 2-variant smoke, 20-variant trial and 100-positive/100-control pilot were gated before paid inference. GI_API_KEY was absent at preflight. Authenticated positive control, quotas and costs could not be validated. Zero paid requests were dispatched. Empty candidate/control and prediction files explicitly represent zero scientifically eligible observations; provisional PIP screening is saved separately and is not a completed cohort. Controls cannot be matched on frequency, GC or locus without the necessary data. No synthetic test outputs enter these metrics.

## Dataset and model contract

Dataset: [10.5281/zenodo.18221260](https://zenodo.org/records/18221260), licence cc-by-4.0. The immutable record inventory, MD5 validation and independently calculated SHA256 values are in manifest.json. Full member audit is in data_audit.md and archive_profiles.json. No raw donor data were used.

Live public catalog default: g0-expression-8192. Primary model g0-expression-8192 is confirmed available in the catalog, but no deployed weights/revision was obtained from an authenticated prediction. Catalog and OpenAPI bytes and SHA256 hashes are saved. Current contract requires bearer authentication, experimental context in options.description, explicit zero-based TSS for input beyond 9,198 bp, and transcriptionally oriented sequences. The long model recommends 40,960 bp per flank (81,920 total); actual scored_window must be read from each response. Nominal input coverage is insufficient. The legacy short model's reported slice is not guaranteed token coverage and is excluded from the primary adapter.

## Independent real reference checks

The separate smoke_reference_audit.json records official Ensembl GRCh38 matches for 1:36478697:C:A (CSF3R), 16:56971389:C:T (CETP), and 16:57020104:G:T (CETP coverage probe). All three genomic REF checks passed. CSF3R is minus-strand in the fetched gene record. Using current Ensembl gene boundaries only, distances are 4,625, 9,467 and 58,182 bp respectively. The distal probe exceeds the nominal 40,960-bp flank; this is not a measurement of GI's effective coverage. These checks do not verify published effect alleles, GENCODE v44 canonical transcript TSSs, or live inference. Snapshot hashes are retained; no DNA is committed.

## Methods implemented and validation scope

Streaming safe archive inspection, source checksum verification, deterministic PIP screening by seed, explicit cell mapping, canonical transcript TSS parsing with ambiguous TSS rejection, indexed FASTA retrieval, genomic-plus SNV edits followed by reverse complementation, explicit TSS indexing, strict scored-window checks, ALT minus REF effects without minus-strand sign flips, content-addressed caching and checkpointing, bounded retry handling and cumulative attempt budget, locus-preserving splits and bootstrap/sign-flip uncertainty are implemented. CLI gates prevent scoring unverified cohorts. Tests exercise synthetic alleles, geometry, errors, caching, retries, control matching and offline integration; their scores are not scientific findings. See test_results.txt for the executed test outcome.

## Interpretation and limitations

Small summary archives lack beta/SE/effect allele and allele frequency; published p-values and PIP cannot substitute for observed effects. Source code indicates GENCODE v44; exact reference assembly sequence and per-gene TSS remain unvalidated on the real cohort. PIP measures statistical fine-mapping uncertainty under LD, not intervention-confirmed causality. Insignificance is not a true negative. SAIGE-QTL effect units need not match log(TPM+1); rank/direction comparison and calibration sensitivity are necessary. Gene/locus clustering is required for uncertainty. Assay and cell-context mismatch, unequal cell power, multiple TSSs, coverage exclusions and LD can all bias results. G0 catalog training sources include GTEx, ENCODE and cellxgene; TenK10K donor/training overlap cannot be excluded. No formal training holdout is established.

No efficacy plots are generated because there are no eligible measured predictions. Distance/frequency baselines, fine-mapping ranks, cross-cell contrasts and external model comparisons remain unexecuted. The implementation is tested offline; it has not been validated by successful live inference.

## Reproduction and residual work

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install --no-deps -e .
.venv/bin/python -m pytest -q
.venv/bin/gxi-tenk10k run --config configs/pilot.yaml --out results/pilot-20261029 --max-requests 1600
```

Resume requires GI_API_KEY in the process environment, validated quota/cost evidence, verified GRCh38 reference and GENCODE v44 TSS annotation, and full association estimates with explicitly confirmed effect alleles. These are independent gates: supplying a key alone does not qualify provisional rows. The large archive remains an explicit opt-in. Stage attempts and exclusions are persisted. No thresholds were tuned against model outcomes.
