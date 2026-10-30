# GI × TenK10K: executed scientific benchmark

## Central result and scope

Under the explicitly amended 40,960-bp input geometry, 100 candidate pairs were scored across 4 discovery cell types. After 23 near-zero/tied-effect abstentions and 0 missing observations, sign concordance was **67.5%**, with locus-bootstrap 95% CI **56.4–77.9%** from 77 eligible effects and 74 clusters. Spearman rho across all 100 candidate predictions, including zero deltas, was 0.370; the locus sign-flip permutation p-value was 0.000999000999000999. These are observational, conditional-on-nonzero, amended-geometry results; they are not causal validation or an unchanged preregistered primary endpoint.

These amended candidate results suggest directional signal among non-abstaining effects. They do not establish robust cell specificity, which was not supported by the supplementary contrast analysis, or validate the unchanged original primary protocol.

The original 81,920-bp smoke REF/ALT queries covered both variants but reported endpoint differences of 4–7 bp. Strict same-window eligibility therefore failed. Before the trial, the input was shortened to 40,960 bp (20,480-bp flanks), which the same long model scored in full for both alleles. The original 100-candidate pass before resolving matching eligibility had concordance 73.1% (95% CI 62.2–83.3%; 78 non-abstaining); that sensitivity cohort is retained rather than hidden. Original-input and amended-input smoke predictions are saved separately. The CETP delta changed sign between these geometries; this is a material sensitivity and a reason against claiming robust success of the original benchmark. No efficacy threshold, locus ranking or cell choice was tuned using trial or pilot predictions.

## Completed stages and actual denominators

```json
{
  "smoke": {
    "status": "complete",
    "pairs": 2,
    "calls_cumulative": 443
  },
  "trial": {
    "status": "complete",
    "candidate_pairs": 20,
    "cells": 4,
    "calls_cumulative": 443
  },
  "pilot": {
    "status": "complete",
    "candidate_pairs": 100,
    "control_pairs": 100,
    "cells": 4,
    "calls_cumulative": 443,
    "new_calls_this_invocation": 0
  }
}
```

The trial used the first 20 eligible candidate pairs across four discovery cells. Its non-abstaining sign concordance was 0.6 (95% CI [0.3333333333333333, 0.8666666666666667]; n=15). The pilot reused trial results through the content-addressed cache. A single discovery context was used for each pilot pair; the supplementary 20-pair cross-cell analysis is described below. Candidate eligibility required gene-level ACAT BH q<=0.05, SuSiE PIP>=0.9, nonpalindromic biallelic SNVs, one canonical GENCODE TSS, verified REF, and exact matching REF/ALT scored intervals. Candidates were drawn in deterministic hash order with round-robin cells, with coverage-based exclusions rather than output-based selection.

| Discovery cell | All pairs | Non-abstaining pairs | Sign concordance | Locus 95% CI |
|---|---:|---:|---:|---|
| B_naive | 18 | 12 | 41.7% | 16.7–66.7% |
| CD14_Mono | 27 | 19 | 63.2% | 42.1–84.2% |
| CD16_Mono | 28 | 23 | 82.6% | 65.2–95.7% |
| CD4_Naive | 27 | 23 | 69.6% | 47.8–87.0% |

The first candidate-only pass produced 100 candidates but only 76 covariate-matched controls. That pass is retained under initial_unmatched_pass/. Match availability was then checked before inference using the same source order and matching tolerances, with no beta or GI-delta selection. Candidates without eligible controls were replaced from the prescreened queue to complete the joint matched cohort. This is an outcome-independent selection adjustment resolved after the first pass and is an additional protocol deviation, not an unchanged prespecified cohort.

All exclusions, including attempted input/coverage failures and bounded control-matching failures, are in exclusions.tsv. Prescreening retained 300 distinct genes before coverage/matching; expanded_prescreen.parquet and full_summary_audit.json preserve source denominators. Counts of source rows scanned in the four cells: 63,986,931. This is a targeted summary-data extraction, not a donor-level dataset.

## Controls, magnitude, baselines and uncertainty

The 100 controls are MAF-, TSS-distance-, local-GC- and gene/cell-matched uncertain comparators, not definite negatives. Matching tolerances were MAF difference <=0.05, absolute TSS-distance difference <=5,000 bp and 201-bp GC fraction difference <=0.1. Known PIP had to be <0.5; variants unreported in the 95% credible sets retained missing PIP, never an invented zero. No insignificance cutoff made a true-negative claim. PIP reporting: {'unreported': 98, 'known': 2, 'note': 'Unreported PIP is not zero.'}.

Control sign concordance: 0.5166666666666667 (95% CI [0.3933523200889136, 0.6440677966101694]; n=60, abstentions=40).

Matched magnitude prioritization:

```json
{
  "epsilon": 1e-06,
  "pairs": 100,
  "ties": 24,
  "non_tied": 76,
  "fraction_candidate_larger": 0.7894736842105263,
  "ci95": [
    0.6962025316455697,
    0.8783783783783784
  ],
  "interpretation": "Descriptive high-PIP candidate versus low/unreported-PIP uncertain comparator; not causal AUC"
}
```

This is a proxy comparison with uncertain controls, not a validated causal AUC or within-credible-set causal ranking. Distance baseline: {'spearman_negative_distance_vs_abs_observed_beta': -0.09248152559755274, 'note': 'Distance cannot predict allele sign; chance sign baseline is 0.5.'}. Chance directional agreement is 0.5. Additional external sequence-model baselines were not executed. GI effects and SAIGE-QTL beta have different measurement scales: observed effects come from a covariate-adjusted count model, while predictions are log(TPM+1). Robust regression and rank correlation are descriptive; no calibrated magnitude equivalence is asserted.

The primary analysis uses epsilon=1e-6. Saved metrics include sensitivity at 0, 1e-6, 1e-5 and 1e-4, distance strata, PIP>=0.99 sensitivity and disjoint locus-based split labels. API scores are rounded in the returned responses; zero deltas and abstention denominators must accompany sign agreement. Bootstrap resamples loci with all associated cells/variants retained (2,000 replicates); the null flips signs jointly by locus (1,000 replicates). Conservative loci merge overlapping +/-100-kb TSS windows, which is a genomic proxy rather than a measured LD block. There was no fitting or calibration on the holdout split. Multiple stratified analyses are descriptive and uncorrected.

## Cell specificity

Comparison stage:

```json
{
  "pairs": 20,
  "observations": 63,
  "cells": 4,
  "contrasts": 81,
  "exclusions": 17,
  "new_requests": 86,
  "cumulative_requests": 443
}
```

Cell-contrast metrics:

```json
{
  "n_total": 81,
  "n_missing": 0,
  "n_abstain": 10,
  "n_scored": 71,
  "n_loci": 17,
  "n_genes": 17,
  "n_cells": 6,
  "sign_concordance": 0.4507042253521127,
  "ci95": [
    0.2784634317862166,
    0.5897727272727271
  ],
  "spearman": -0.0402276640573293,
  "permutation_p": 0.7582417582417582,
  "epsilon": 1e-06,
  "robust_slope": -0.19553255813953613,
  "spearman_n": 81,
  "spearman_nonabstaining": -0.03511419677698817,
  "chance_baseline": 0.5,
  "calibration_sensitivity": {
    "0": {
      "retained": 71,
      "sign_concordance": 0.4507042253521127
    },
    "1e-06": {
      "retained": 71,
      "sign_concordance": 0.4507042253521127
    },
    "1e-05": {
      "retained": 71,
      "sign_concordance": 0.4507042253521127
    },
    "0.0001": {
      "retained": 70,
      "sign_concordance": 0.44285714285714284
    }
  },
  "n_biological_cells": 4
}
```

These compare allele-oriented differences in beta with differences in GI delta across prespecified blood lineages for the first 20 candidate pairs. Repeated contrasts are clustered by locus. Cell power, phenotype scale and within-donor covariance can differ; standard errors alone do not establish significant cross-cell heterogeneity. A four-cell cohort by itself is not evidence of cell specificity.

## Dataset, reference, inference and reproducibility

Dataset creator: Anna S. E. Cuomo (Zenodo metadata). Dataset: [Zenodo 18221260](https://zenodo.org/records/18221260), DOI 10.5281/zenodo.18221260, first version, CC-BY-4.0. The two initial ZIPs and the 16,559,137,990-byte full common 100-kb summary ZIP were downloaded and verified against their provided MD5 checksums, with SHA256 independently recorded. Full archive SHA256: 9be398a12c2b4e1b8118b03dca1123f7a1b3036f79e0322c1f392beb6505e5e0. Only four full-summary members were scanned; unrelated large releases were not downloaded. Small archives use comma-separated .tsv files; the full archive uses tab-separated tables with an extra row-index column.

Effect orientation was established from [official SAIGE-QTL documentation](https://weizhou0.github.io/SAIGE-QTL-doc/docs/single_step2.html): Allele1 is REF, Allele2 is ALT/effect allele, BETA is the effect of Allele2. Actual IDs and columns were checked against the genomic FASTA. Genomic edits occur before reverse complementation; minus-strand genes do not induce beta sign flips. Sample sizes N, beta, SE and missing values are retained in the local observation pool. The original record links a TenK10K preprint; no unverified bibliographic or total-donor claims are needed for these row-level results.

Reference: official GRCh38 primary assembly FASTA and GENCODE v44 basic GTF, both validated with the release MD5SUMS and independently SHA256-hashed (reference_manifest.json). Canonical transcript tags must identify a unique TSS. The annotation supports minus-strand CSF3R; independent real base checks passed at all three protocol smoke/coverage positions. The distal CETP probe is outside the chosen scored input and remains a coverage probe.

Model: g0-expression-8192, current live GI catalog and OpenAPI archived with hashes. Deployed immutable weight version is not exposed. Bearer authentication, context, oriented TSS index, response unit, sequence length and scored-window metadata were validated on real responses. The exact timing and values of both adjustments are in execution_amendments.json and configs/live_geometry_amendment.yaml. Each pilot response had identical 40,960-bp request geometry per allele and strict scored-window equality. Exact sequence hashes and cache keys are recorded; no DNA or key is committed. Exact biological context strings are present in cached local responses and the sequence-free measured table's model/coordinate fields.

The initial HBB K562 fixture returned 0.5954 log(TPM+1), and the longer reference query returned 1.508. Thus the quickstart's qualitative high-expression expectation was not independently calibrated. These calls verified the contract and exposed window dependence; they are not quantitative ground truth. CSF3R reference expression in the amended monocyte smoke was 4.9446, and its allele effect was zero. The smoke is never pooled with candidate efficacy results.

Total paid POST attempts including diagnostics: 444. No account-balance/pricing endpoint was published; observed burst capacity was 200 (documented sustained ceiling 1,200 rpm). This is rate capacity, not remaining monetary quota. The user authorized the pilot; calls were bounded below 800, with a 760-attempt cumulative cache-ledger cap plus one initial preflight call. All quota/authentication errors would stop expansion, and no retries were unbounded. No quota-balance or cost amount is invented.

Timestamp: 2026-10-29T23:07:46.372998+00:00. Starting Git revision: eb2b07ed5fba0661fb0260ae4fc0776ca137f930. Final implementation fingerprint: 66d1cdead6d3a521482cecb192822d44b8887c51108198d555cda6014659af50. Full package versions, per-file code hashes, source hashes and API response provenance are retained. Large sources and caches remain ignored by Git; measured_predictions.tsv is the small sequence-free reproducibility artifact.

## Limitations and interpretation

H1 is assessed descriptively under the geometry amendment; the original strict-geometry endpoint was not completed. H2 is exploratory on 20 pairs with repeated comparisons. H3 is limited by uncertain and frequently unreported control PIP; no causal ranking is established. LD, observational ascertainment, gene/cell power, multiple TSSs and assay mismatch remain limitations. GI training includes GTEx, ENCODE and cellxgene; donor/training overlap with TenK10K cannot be excluded, so this is zero-shot use without a formally proven training holdout. Gene-window clusters are not measured LD blocks. The result does not imply donor-level prediction or intervention-confirmed effects. A defensible follow-up would freeze an immutable model revision and input geometry, verify biological calibration, prespecify independent loci and collect a formal held-out replication.

## Reproduction

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install --no-deps -e .
.venv/bin/python -m pytest -q
# Initial discovery/audit is preserved in results/pilot-20261029.
# Export GI_API_KEY securely before live commands.
.venv/bin/gxi-tenk10k run --verified-live --out results/live-20261029 --max-requests 760
.venv/bin/python scripts/compare_trial_cells.py
.venv/bin/gxi-tenk10k evaluate --out results/live-20261029
.venv/bin/gxi-tenk10k report --out results/live-20261029
```

Use scripts and manifests for the source/reference preparation commands recorded in docs/IMPLEMENTATION.md. Cache reuse avoids repeating successful predictions. Evaluation, figures and the report require no GI calls. See test_results.txt for the executed offline test outcome. The standalone figure is figures/benchmark.png.
