"""Scientific report for actually executed, amended-geometry predictions."""
import json
from pathlib import Path


def live_report(out):
    out=Path(out)
    load=lambda name:json.loads((out/name).read_text())
    metrics=load('metrics.json');stages=load('stages.json');env=load('environment.json')
    audit=load('full_summary_audit.json');amendment=load('geometry_amendment.json')
    comparison=load('comparison_stage.json') if (out/'comparison_stage.json').exists() else None
    contrast=load('cell_contrast_metrics.json') if (out/'cell_contrast_metrics.json').exists() else None
    c=metrics['candidate'];control=metrics['control'];trial=load('trial_metrics.json')
    pct=lambda x: 'NA' if x is None else f'{100*x:.1f}%'
    ci=lambda x: 'NA' if x is None else f'{100*x[0]:.1f}–{100*x[1]:.1f}%'
    cellrows='\n'.join(f"| {cell} | {m['n_total']} | {m['n_scored']} | {pct(m['sign_concordance'])} | {ci(m['ci95'])} |" for cell,m in metrics['candidate_by_cell'].items())
    text=f'''# GI × TenK10K: executed scientific benchmark

## Central result and scope

Under the explicitly amended 40,960-bp input geometry, {c['n_total']} candidate pairs were scored across {c['n_cells']} discovery cell types. After {c['n_abstain']} near-zero/tied-effect abstentions and {c['n_missing']} missing observations, sign concordance was **{pct(c['sign_concordance'])}**, with locus-bootstrap 95% CI **{ci(c['ci95'])}** from {c['n_scored']} eligible effects and {c['n_loci']} clusters. Spearman rho across all {c['spearman_n']} candidate predictions, including zero deltas, was {c['spearman']:.3f}; the locus sign-flip permutation p-value was {c['permutation_p']}. These are observational, conditional-on-nonzero, amended-geometry results; they are not causal validation or an unchanged preregistered primary endpoint.

These amended candidate results suggest directional signal among non-abstaining effects. They do not establish robust cell specificity, which was not supported by the supplementary contrast analysis, or validate the unchanged original primary protocol.

The original 81,920-bp smoke REF/ALT queries covered both variants but reported endpoint differences of 4–7 bp. Strict same-window eligibility therefore failed. Before the trial, the input was shortened to 40,960 bp (20,480-bp flanks), which the same long model scored in full for both alleles. The original 100-candidate pass before resolving matching eligibility had concordance 73.1% (95% CI 62.2–83.3%; 78 non-abstaining); that sensitivity cohort is retained rather than hidden. Original-input and amended-input smoke predictions are saved separately. The CETP delta changed sign between these geometries; this is a material sensitivity and a reason against claiming robust success of the original benchmark. No efficacy threshold, locus ranking or cell choice was tuned using trial or pilot predictions.

## Completed stages and actual denominators

```json
{json.dumps(stages,indent=2)}
```

The trial used the first 20 eligible candidate pairs across four discovery cells. Its non-abstaining sign concordance was {trial['sign_concordance']} (95% CI {trial['ci95']}; n={trial['n_scored']}). The pilot reused trial results through the content-addressed cache. A single discovery context was used for each pilot pair; the supplementary 20-pair cross-cell analysis is described below. Candidate eligibility required gene-level ACAT BH q<=0.05, SuSiE PIP>=0.9, nonpalindromic biallelic SNVs, one canonical GENCODE TSS, verified REF, and exact matching REF/ALT scored intervals. Candidates were drawn in deterministic hash order with round-robin cells, with coverage-based exclusions rather than output-based selection.

| Discovery cell | All pairs | Non-abstaining pairs | Sign concordance | Locus 95% CI |
|---|---:|---:|---:|---|
{cellrows}

The first candidate-only pass produced 100 candidates but only 76 covariate-matched controls. That pass is retained under initial_unmatched_pass/. Match availability was then checked before inference using the same source order and matching tolerances, with no beta or GI-delta selection. Candidates without eligible controls were replaced from the prescreened queue to complete the joint matched cohort. This is an outcome-independent selection adjustment resolved after the first pass and is an additional protocol deviation, not an unchanged prespecified cohort.

All exclusions, including attempted input/coverage failures and bounded control-matching failures, are in exclusions.tsv. Prescreening retained 300 distinct genes before coverage/matching; expanded_prescreen.parquet and full_summary_audit.json preserve source denominators. Counts of source rows scanned in the four cells: {sum(x['rows_scanned'] for x in audit['scanned_members']):,}. This is a targeted summary-data extraction, not a donor-level dataset.

## Controls, magnitude, baselines and uncertainty

The {control['n_total']} controls are MAF-, TSS-distance-, local-GC- and gene/cell-matched uncertain comparators, not definite negatives. Matching tolerances were MAF difference <=0.05, absolute TSS-distance difference <=5,000 bp and 201-bp GC fraction difference <=0.1. Known PIP had to be <0.5; variants unreported in the 95% credible sets retained missing PIP, never an invented zero. No insignificance cutoff made a true-negative claim. PIP reporting: {metrics['control_pip_reporting']}.

Control sign concordance: {control['sign_concordance']} (95% CI {control['ci95']}; n={control['n_scored']}, abstentions={control['n_abstain']}).

Matched magnitude prioritization:

```json
{json.dumps(metrics.get('matched_magnitude_priority'),indent=2)}
```

This is a proxy comparison with uncertain controls, not a validated causal AUC or within-credible-set causal ranking. Distance baseline: {metrics['distance_baseline']}. Chance directional agreement is 0.5. Additional external sequence-model baselines were not executed. GI effects and SAIGE-QTL beta have different measurement scales: observed effects come from a covariate-adjusted count model, while predictions are log(TPM+1). Robust regression and rank correlation are descriptive; no calibrated magnitude equivalence is asserted.

The primary analysis uses epsilon=1e-6. Saved metrics include sensitivity at 0, 1e-6, 1e-5 and 1e-4, distance strata, PIP>=0.99 sensitivity and disjoint locus-based split labels. API scores are rounded in the returned responses; zero deltas and abstention denominators must accompany sign agreement. Bootstrap resamples loci with all associated cells/variants retained (2,000 replicates); the null flips signs jointly by locus (1,000 replicates). Conservative loci merge overlapping +/-100-kb TSS windows, which is a genomic proxy rather than a measured LD block. There was no fitting or calibration on the holdout split. Multiple stratified analyses are descriptive and uncorrected.

## Cell specificity

Comparison stage:

```json
{json.dumps(comparison,indent=2)}
```

Cell-contrast metrics:

```json
{json.dumps(contrast,indent=2)}
```

These compare allele-oriented differences in beta with differences in GI delta across prespecified blood lineages for the first 20 candidate pairs. Repeated contrasts are clustered by locus. Cell power, phenotype scale and within-donor covariance can differ; standard errors alone do not establish significant cross-cell heterogeneity. A four-cell cohort by itself is not evidence of cell specificity.

## Dataset, reference, inference and reproducibility

Dataset creator: Anna S. E. Cuomo (Zenodo metadata). Dataset: [Zenodo 18221260](https://zenodo.org/records/18221260), DOI 10.5281/zenodo.18221260, first version, CC-BY-4.0. The two initial ZIPs and the 16,559,137,990-byte full common 100-kb summary ZIP were downloaded and verified against their provided MD5 checksums, with SHA256 independently recorded. Full archive SHA256: {audit['source_sha256']}. Only four full-summary members were scanned; unrelated large releases were not downloaded. Small archives use comma-separated .tsv files; the full archive uses tab-separated tables with an extra row-index column.

Effect orientation was established from [official SAIGE-QTL documentation]({audit['effect_allele_source'].split(' ')[0]}): Allele1 is REF, Allele2 is ALT/effect allele, BETA is the effect of Allele2. Actual IDs and columns were checked against the genomic FASTA. Genomic edits occur before reverse complementation; minus-strand genes do not induce beta sign flips. Sample sizes N, beta, SE and missing values are retained in the local observation pool. The original record links a TenK10K preprint; no unverified bibliographic or total-donor claims are needed for these row-level results.

Reference: official GRCh38 primary assembly FASTA and GENCODE v44 basic GTF, both validated with the release MD5SUMS and independently SHA256-hashed (reference_manifest.json). Canonical transcript tags must identify a unique TSS. The annotation supports minus-strand CSF3R; independent real base checks passed at all three protocol smoke/coverage positions. The distal CETP probe is outside the chosen scored input and remains a coverage probe.

Model: g0-expression-8192, current live GI catalog and OpenAPI archived with hashes. Deployed immutable weight version is not exposed. Bearer authentication, context, oriented TSS index, response unit, sequence length and scored-window metadata were validated on real responses. The exact timing and values of both adjustments are in execution_amendments.json and configs/live_geometry_amendment.yaml. Each pilot response had identical 40,960-bp request geometry per allele and strict scored-window equality. Exact sequence hashes and cache keys are recorded; no DNA or key is committed. Exact biological context strings are present in cached local responses and the sequence-free measured table's model/coordinate fields.

The initial HBB K562 fixture returned 0.5954 log(TPM+1), and the longer reference query returned 1.508. Thus the quickstart's qualitative high-expression expectation was not independently calibrated. These calls verified the contract and exposed window dependence; they are not quantitative ground truth. CSF3R reference expression in the amended monocyte smoke was 4.9446, and its allele effect was zero. The smoke is never pooled with candidate efficacy results.

Total paid POST attempts including diagnostics: {env['paid_post_attempts']}. No account-balance/pricing endpoint was published; observed burst capacity was 200 (documented sustained ceiling 1,200 rpm). This is rate capacity, not remaining monetary quota. The user authorized the pilot; calls were bounded below 800, with a 760-attempt cumulative cache-ledger cap plus one initial preflight call. All quota/authentication errors would stop expansion, and no retries were unbounded. No quota-balance or cost amount is invented.

Timestamp: {env['timestamp_utc']}. Starting Git revision: {env['git_commit']}. Final implementation fingerprint: {env['implementation_tree_hash']}. Full package versions, per-file code hashes, source hashes and API response provenance are retained. Large sources and caches remain ignored by Git; measured_predictions.tsv is the small sequence-free reproducibility artifact.

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
'''
    (out/'report.md').write_text(text)
    return out/'report.md'
