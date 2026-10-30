"""Generate a scientific status report; absent measurements stay null."""
import json
from pathlib import Path


def report(out):
    out = Path(out)
    if (out/'live_run_manifest.json').exists():
        from .live_report import live_report
        return live_report(out)
    stages = json.loads((out/'stages.json').read_text()) if (out/'stages.json').exists() else {}
    metrics = json.loads((out/'metrics.json').read_text()) if (out/'metrics.json').exists() else {'n_scored': 0}
    manifest = json.loads((out/'manifest.json').read_text())
    selection = json.loads((out/'selection.json').read_text()) if (out/'selection.json').exists() else {}
    models = json.loads((out/'gi_models.json').read_text())
    env = json.loads((out/'environment.json').read_text())
    text = f'''# GI × TenK10K zero-shot benchmark: executed status report

## Central finding

No live GI variant-effect measurements were obtained. H1 (direction), H2 (cell specificity), and H3 (PIP prioritization) remain untested. Sign concordance, magnitude correlation, confidence intervals and permutation p-values are unavailable, not zero. No inference about GI performance or causality is supported.

## Execution

Timestamp: {env['timestamp_utc']}. Source revision: {env['git_commit']}.
Configuration SHA256: {env['config_sha256']}.
Final implementation tree fingerprint: {env.get('implementation_tree_hash', 'unavailable')}.
The code was uncommitted at execution; environment.json records per-file hashes to identify the actual implementation independently of the starting revision.

Stage status:

```json
{json.dumps(stages, indent=2)}
```

Selection flow:

```json
{json.dumps(selection, indent=2)}
```

Metrics (only measured data are eligible):

```json
{json.dumps(metrics, indent=2)}
```

Implementation scope is incomplete: automatic normalization of the unaudited full association archive, cross-cell differential analysis, fine-mapping ranking and scientific plot generation remain unimplemented. This report documents an executed preflight/data-audit/offline-validation run, not a completed efficacy benchmark. See docs/IMPLEMENTATION.md for these limits and the verified-input interface.

The 2-variant smoke, 20-variant trial and 100-positive/100-control pilot were gated before paid inference. GI_API_KEY was absent at preflight. Authenticated positive control, quotas and costs could not be validated. Zero paid requests were dispatched. Empty candidate/control and prediction files explicitly represent zero scientifically eligible observations; provisional PIP screening is saved separately and is not a completed cohort. Controls cannot be matched on frequency, GC or locus without the necessary data. No synthetic test outputs enter these metrics.

## Dataset and model contract

Dataset: [{manifest['doi']}](https://zenodo.org/records/18221260), licence {manifest['license']['id']}. The immutable record inventory, MD5 validation and independently calculated SHA256 values are in manifest.json. Full member audit is in data_audit.md and archive_profiles.json. No raw donor data were used.

Live public catalog default: {models['default_model']}. Primary model {models['default_model']} is confirmed available in the catalog, but no deployed weights/revision was obtained from an authenticated prediction. Catalog and OpenAPI bytes and SHA256 hashes are saved. Current contract requires bearer authentication, experimental context in options.description, explicit zero-based TSS for input beyond 9,198 bp, and transcriptionally oriented sequences. The long model recommends 40,960 bp per flank (81,920 total); actual scored_window must be read from each response. Nominal input coverage is insufficient. The legacy short model's reported slice is not guaranteed token coverage and is excluded from the primary adapter.

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
'''
    (out/'report.md').write_text(text)
    return out/'report.md'
