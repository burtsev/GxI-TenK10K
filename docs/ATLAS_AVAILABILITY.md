# Atlas availability and nonblocking execution contract

Version: atlas-v1. This is a mandatory implementation requirement for Codex, not a claim that the existing followup CLI already implements it.

## Governing requirement

**If an input is unavailable, skip only the smallest dependent analysis and continue every independent runnable part of the study. Always write the final report, including skipped sections.** Do not pause the whole run to ask the user for optional files, credentials, permissions, server instrumentation, external model access or clarification that only affects one branch. Do not initiate permission requests or wait for collaborators. No fabricated substitute data, unauthorized access or silent relabelling of another cohort as TenK10K.

Optional means conditional on availability, not optional to execute when available. Audit each configured input and run every supported analysis. A missing cell, chromosome, column, gene annotation or archive should not skip a whole benchmark when valid remaining partitions suffice.

## Availability checks

Before expensive inference, build a persistent source/capability manifest. Use verified local files first, then the documented public sources or explicitly authorized configured inputs. Null paths mean not supplied; they are not a reason to wait. Existence alone does not establish authorization or schema validity. Record check time, access/licence basis, schema, checksum, version, attempted source and reason for unavailable status. Redact secrets and controlled paths/identifiers from public exports.

Bound transient network checks to at most three attempts with timeout/backoff and Retry-After. Distinguish transient failure from 401/403 permission denial, nonexistent files, absent fields, unsupported features and integrity failures. Do not repeatedly request a denied/absent source for each variant. Cache unavailability for the run, and recheck once on a deliberate resume/recheck stage. If the endpoint is temporarily unavailable, use already verified compatible caches; otherwise record the affected live branch as unavailable and continue.

## Dependency matrix

| Analysis | Minimum usable inputs | Missing-data behavior |
|---|---|---|
| Inventory, universe and coverage audit | Available manifests/schema and any available source partitions | Always run; describe inaccessible sources; no inference required |
| Reference-expression predictions | Valid sequence/TSS/cell contexts and live GI or compatible reference cache | Score available genes; predictions without measured expression are not expression accuracy |
| Expression accuracy | Approved aggregate expression OR approved donor counts plus usable metadata/normalization, and matched GI predictions | Skip expression accuracy only; continue variant benchmarks |
| Donor reliability / donor bootstrap | Donor-resolved observed expression and donor labels | Skip these statistics; aggregate expression accuracy can still run |
| Donor-personalized expression | Approved donor genotypes/haplotypes plus expression and validated personalized inputs | Skip this extension; no effect on population-average or variant modules |
| Common variant scoring | Available common variant definitions, validated genomic inputs, contexts and live GI or compatible caches | Use all supported partitions; no PIP, matched controls or expression means required |
| Rare variant scoring | Available rare variant definitions with sufficient allele information, valid genomic inputs, contexts and GI/cache | Missing rare source skips rare scoring only; common analyses continue |
| Allelic direction / signed or absolute magnitude | Matched observed ALT-oriented beta and valid predicted delta | Run available pairs; missing SE skips precision-dependent strata, not these base metrics |
| Cell-effect profiles / contrasts | The same variant-gene in >=2 compatible measured and predicted cells (>=5 for profile correlation) | Missing cells remove only their comparisons; do not require all 28 |
| Formal cross-cell heterogeneity or equivalence | Compatible uncertainty/covariance and prespecified test/margin | Skip formal tests; retain descriptive contrasts, profiles and precision sensitivity when possible |
| Statistical-label AUROC/AUPRC | Explicit high and low PIP for the applicable gene-cell association, plus scores | Exclude unknown PIP from binary truth; if either class absent, report not_estimable; do not block direction/magnitude |
| Credible-set / within-signal ranking | Candidate membership and relevant scores, plus valid probabilities when a probability-weighted metric is used | Skip absent membership/probability metrics only; label partial-set coverage explicitly |
| Matched comparison ranking | Available matched comparison variants and scores | Unmatched candidates remain in the core atlas; no rematching based on GI outputs |
| Rare gene-cell enrichment | Rare gene tests, compatible tested-set membership and scored rare variants | Skip set-specific enrichment when membership/tests absent; retain rare single-variant work |
| GI-weighted rare association retesting | Approved donor-level inputs OR sufficient compatible score/covariance summaries | Skip retesting only; gene p-values alone are not sufficient |
| Expression-based or annotated subgroups | Relevant observed expression or annotation fields, plus base metrics | Skip the absent dimension only; other cell/gene/frequency/evidence strata continue |
| Joint expression vs variant performance | Expression-accuracy outputs and variant-performance outputs | Skip joint analysis if expression is unavailable |
| External model comparison | Legitimate executable/licensed scores and compatible targets | Skip each unavailable comparator independently; no GI block |
| Full-precision/internal representation diagnostics | Authorized raw precision/token/embedding capability | Continue API-visible scoring with precision limitation; do not infer hidden representations |
| Holdout-specific claims | Verifiable historical and new split provenance | Retain descriptive atlas; label holdout independence unverified rather than stopping analysis |
| Final report and status manifest | Any completed, partial, skipped or failed module records | Always run, even when no scientific metric is estimable |

## Module status model

Every enabled module/partition ends in one of:
- `completed`: all eligible available inputs handled, outputs validated.
- `partial`: usable outputs exist, but coverage, partitions or computations are incomplete; residual work and failed/skipped inputs are explicit.
- `skipped_unavailable`: no usable required input/capability for this unit after bounded checks.
- `skipped_dependency`: a required producer has no usable outputs; name the producer and reason.
- `not_estimable`: inputs exist but a statistic is mathematically undefined or lacks sufficient observations/classes; not a zero result.
- `failed_qc`: data/model integrity checks failed; quarantine affected outputs. Unrelated modules continue.
- `failed_execution`: unexpected code/operational error; preserve trace and checkpoints. Unrelated modules continue, but the overall result must expose the error.

Transient `pending`/`running` states are permitted while executing; a final report may retain `pending` only if the run was explicitly interrupted and labelled partial. An error is not silently converted into an unavailable-data skip.

A partial producer may supply validated partitions to downstream analyses. The dependency engine must inspect output availability, not simply require every upstream status to equal completed. Missing optional covariates or annotations cannot make an entire producer mandatory. Report which analysis assumptions change when a covariate is absent.

Run-level statuses: `completed`, `completed_with_skips`, `partial`, `completed_with_errors`. A run with planned unavailable-data skips can exit normally after reporting; unexpected execution/QC errors must be surfaced in structured status and a nonzero final exit code after independent work completes. When nothing scientific is runnable, still emit a report stating that no scientific result was produced; do not claim successful validation. Process interruption/resource exhaustion writes a partial report and resumable pending ledger when possible.

Do not globally suppress exceptions with `except Exception: pass`. Catch known unavailable-input cases as typed outcomes, record unexpected failures at module boundaries, and finalize the report in a `finally` path. Fatal disk corruption or unsafe execution may require stopping writes, but never publish unvalidated outputs. This is different from an ordinary missing optional dataset.

## Required status/provenance outputs

`availability.json`: source/capability ID, available/partial/unavailable/failed_qc, reason code, checked time, source/version/hash, authorization basis, available partitions, unavailable partitions and verified cache substitutes.

`module_status.tsv` and machine-readable JSON: run_id, module_id, partition_id, status, required_input_ids, missing_inputs, reason_code, detail, eligible_n, processed_n, scored_n, excluded_n, failed_n, started_utc, finished_utc, output_paths, output_hashes, next_resume_action. Undefined counts are null, never invented zeros.

`coverage.tsv`: source row counts and unique gene/variant/cell counts separately; scorable, outside window, unsupported, ambiguous, failed, successfully scored. Preserve source-analysis IDs and denominators; no silent denominator changes after a skipped partition.

`skipped_analyses.md`: what could not be answered, why, evidence of the availability check, and what input would enable it. This is a record, not an instruction to wait.

`report.md`: all benchmark headings retained. An unavailable expression source yields an explicit expression-accuracy SKIPPED section, not a blank omission. A rare-enrichment skip does not remove rare single-variant metrics. Report intended full-release coverage and actually attained coverage separately.

## Resume behavior

Persist successful predictions, validated intermediate tables, status and exact cache/protocol hashes after every completed partition. On resume, verify provenance, reuse successes and re-evaluate availability once; run newly enabled or unfinished units only. Never replay all completed work simply because one optional input arrived. New input/outcome data may require a new analysis version and disclosure; do not silently reopen a frozen confirmatory decision.

## Required offline acceptance scenarios

Codex must implement fixture-based tests without paid calls:

1. No expression file / no EGA approval: expression accuracy, donor metrics and expression strata are skipped; common direction, magnitude, cell specificity, available fine-mapping and rare modules complete.
2. Rare archive missing: only rare branches skip; common results and final report are written.
3. SuSiE/PIP absent: variant scoring, direction, magnitude and cell contrasts run; missing PIP is never imputed zero; fine-mapping labels skip.
4. Credible-set membership absent but explicit PIPs present: classification runs where labels exist; set-specific ranking skips.
5. Aggregate expression present without donors: available expression rank metrics run; donor bootstrap/reliability skip.
6. SE or covariance absent: basic beta comparisons run; uncertainty-dependent subsets/formal tests skip.
7. One cell/chromosome/source partition absent: remaining partitions run, counts and partial coverage are preserved.
8. GI auth/network unavailable with some valid cached predictions: evaluate matching cached rows and report partial inference; without predictions still audit sources and write a no-results report.
9. Raw server precision unavailable: API-visible inference continues with quantization warning; no invented raw values.
10. An external comparator is unavailable/licence-unapproved: only that comparator skips.
11. Historical cache mapping incomplete: descriptive atlas proceeds and untouched-holdout claims are disabled.
12. REF mismatch or checksum/allele failure: quarantine affected records/partition as failed_qc, never score them as zero; unaffected work continues and the failure is visible.
13. Constant vectors, too few cells or a single AUROC class: metric is not_estimable with a reason, not zero or an exception that kills the run.
14. An unexpected module exception: record failed_execution; continue independent modules; final report includes the exception and overall errors status.
15. Resume after an optional source appears: only newly enabled work executes, successful request count is unchanged elsewhere, prior caches remain intact.
16. No optional inputs at all: the orchestrator terminates without waiting for user action and writes a complete status report with no fabricated findings.

These scenarios are acceptance requirements for the new atlas implementation. Existing 68-test followup results do not establish that these new scenarios pass.
