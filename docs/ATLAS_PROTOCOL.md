# Comprehensive GI x TenK10K atlas: scientific protocol

Protocol version: atlas-v1, 31 October 2026. Status: implementation specification, not executed results.

## 1. Purpose and governing rule

Evaluate (1) gene-expression prediction, (2) variant effects across the complete available and model-scorable TenK10K release, and (3) performance across cells, genes and other prespecified subsets. This is not another fixed-size high-PIP pilot.

**Missing or inaccessible data must skip only the dependent analysis. Continue every independent runnable module and always produce a report.** Do not wait for EGA permission, a collaborator response, an external-model licence, donor genotypes, full-precision server instrumentation or complete historical caches. See [ATLAS_AVAILABILITY.md](ATLAS_AVAILABILITY.md) for the executable behavior Codex must implement. Never substitute invented data or unverified labels.

Start with [ATLAS_START.md](ATLAS_START.md), implement [ATLAS_CODEX_RUNBOOK.md](ATLAS_CODEX_RUNBOOK.md), and use [../configs/atlas.yaml](../configs/atlas.yaml). Existing pilot/followup protocols and executed outputs remain immutable historical evidence. This protocol supersedes their sample caps, universal high-PIP eligibility, six-cell restriction, global completion gates and quota checks for the new atlas only.

## 2. Benchmark questions and estimands

| Benchmark | Principal metrics | Question |
|---|---|---|
| Gene expression | Spearman across genes within each cell; across cells within each gene; gene-centered cell-profile agreement | Which genes are expressed, at what relative levels, and in which immune cells? |
| Allelic effect direction | Sign concordance with ALT-oriented observed beta, with coverage, zero/abstention and uncertainty counts | Does GI predict whether a variant increases or decreases expression? |
| Effect magnitude | Signed Spearman(delta, beta), separately Spearman(abs(delta), abs(beta)) | Does GI capture direction and rank strong versus weak regulatory effects? |
| Cell-type specificity | Correlation of predicted/observed beta contrasts; within-variant-gene effect profiles across cells | Can GI identify which immune cell types are affected? |
| Fine-mapping | AUROC/AUPRC on explicit statistical labels; within-locus/credible-set rankings | Can GI prioritize candidate causal variants over linked alternatives? |
| Rare variant prioritization | Gene-cell enrichment adjusted for variant-set properties; separately association-test improvement when inputs permit | Can GI improve functional selection of rare regulatory variants? |

All inferential claims must name the analysis population, independent sampling unit, model/input version, available observations and uncertainty. A data-source row is not an independent biological replicate.

## 3. Sources, permissions and complete coverage

Reuse checksummed local sources and source manifests before downloading again. Discover the current inventory of Zenodo record 18221260; record exact version, file/member names, hashes, schemas and analysis identifiers. Relevant source families include common gene-level results, full common cis statistics, SuSiE summaries/credible-set information, rare single-variant results and rare gene-level tests. Inventory both +/-100-kb and +/-1-Mb releases; stream/deduplicate overlapping records and preserve distinct association specifications. Download partitions that add information, not duplicate archives indiscriminately.

Use verified GRCh38 sequence and GENCODE v44 to retain comparability. Record gene IDs, transcript/TSS choice, genomic strand, effect allele, variant normalization and source-analysis identity. Reject ambiguous records with explicit reasons. Do not mix coefficients from incompatible fits.

Expression observations require a separate source: an authorized TenK10K aggregate expression reference or authorized processed counts/metadata. EGA accession EGAD50000002379 is a source lead identified during planning, **not evidence that this environment has access or permission**. Check current terms and authorization for this use. Without an available approved source, skip expression accuracy, donor analyses and observed-expression-based strata; continue summary-statistic benchmarks immediately. Do not request access automatically, wait for approval, or use a different dataset as unlabeled TenK10K ground truth. Never infer expression means or individual donors from eQTL betas.

The atlas must include all available cell types (expected 28; reconcile actual metadata), eligible genes, and common/rare variant-gene-cell observations, independent of significance, PIP, availability of a control or membership in old studies. Include genes without eQTLs in the expression arm when the expression reference permits. Historical rows stay in the descriptive atlas with provenance flags.

Record two denominators separately: the entire inventory and the model-scorable subset. Each record has a terminal coverage status: scored, unavailable observation, unsupported variant, ambiguous mapping, out_of_window, failed_qc or failed_inference. **Out-of-window is never a zero effect.** Missing PIP is never PIP=0; no observed beta is never beta=0. A missing archive/cell/chromosome does not remove the other partitions.

Canonical identity includes assembly, chromosome, position, REF, ALT, gene, cell and source-analysis ID. Deduplicate inference by actual input hashes, not merely by rsID. No significance-based downsampling of the all-data atlas. Execution may be partitioned, but a convenient subset must not be described as all data.

## 4. Common model, precision and allele controls

Use the existing direct GI REST/official-client implementation with GI_API_KEY from environment. The owner reports no total quota for this key. No automatic 790/800/1600-call cap applies. Authentication, documented rate limits, bounded retries, durable checkpoints and operational safeguards still apply; unlimited quota does not imply unlimited throughput or free service.

The reproducibility arm uses g0-expression-8192, 40,960 bases, oriented 0-based TSS index 20,480 and the prior exact per-cell baseline descriptions, subject to schema/model preflight. Verify model ID, actual effective spans and full-precision fields rather than assuming API behavior. Keep REF/ALT model, cell description, TSS anchoring and genomic coverage identical. Freeze descriptions before evaluation, with no gene identifiers or new target-informed prompt search in the main arm.

Apply edits to verified genomic-plus REF, then reverse complement for minus-strand genes. Compare observed effects only after mapping their effect allele to ALT. Preserve source beta and flip provenance. SNVs are the initial validated variant class; add normalized indels only through a separately tested TSS-anchored builder that preserves the required model geometry. Unsupported indels or other classes are reported and skipped without stopping SNVs.

Compute delta = log(TPM_ALT + 1) - log(TPM_REF + 1) from the API's documented expression output, not from rounded display text. Store precision and units explicitly. Compute REF once per exact gene/TSS/sequence/cell/model/context version and reuse it. Persistent cache keys include sequence hashes, context, model/revision, API contract, TSS, geometry and scoring options. Reuse old caches only after equivalence checks; a changed precision/model implementation gets a new namespace.

Before broad inference, test historical development loci for repeatability, mutation survival through preprocessing, effective-window equality and rounding. Request raw regression output if available through authorized server instrumentation. If not available, label raw_precision_unavailable and continue with API-visible values. Do not delay the atlas waiting for internals or claim rounded ties prove encoder insensitivity. Distinguish exact returned ties, known rounding-induced ties, small nonzero effects, missing outputs and invalid geometry. An API-returned tie cannot be labelled an exact raw-model tie without raw output.

Use epsilon=1e-6 for historical comparability and 1e-4/1e-3 sensitivity; preserve raw numeric deltas. Record empirical output resolution. Keep 81,920-bp or repaired-geometry experiments separate, and exclude unequal genomic REF/ALT spans from biological effects. If a model changes or preflight fails, isolate affected scoring; continue data audit, verified-cache evaluation and unrelated modules. Do not pool unknown revisions. Prefer immutable weights/version identifiers; if unavailable, use sentinel repeats and document the residual limitation.

## 5. Expression accuracy module

Run only with authorized observed expression and compatible GI reference predictions. Accept either vetted aggregate gene-cell expression or counts with donor/cell metadata and an explicit normalization description. A counts file without usable labels is not sufficient.

For counts, construct donor-by-cell-type pseudobulk UMI sums, normalize consistently for library size, then form a donor-balanced expression reference. Freeze QC using data quality rather than GI outcomes. Retain donor counts, cell counts, detection fraction and variability. Aggregate-only data can support available rank comparisons; skip donor bootstraps, split-donor reliability or unavailable detection metrics individually.

Do not relabel normalized 10x UMI expression as measured TPM. Primary rank comparisons do not require identical units. Fit any affine/monotone calibration on development genes/loci only, and evaluate MAE/RMSE/R2 on held-out compatible quantities, with a development-derived constant baseline. No post-hoc calibration on the test targets.

Report:
- Across-gene Spearman for each cell type, and an equal-cell summary.
- Across-cell Spearman for each gene with >=5 measured nonconstant cell values, plus the count of ineligible/constant profiles. Require >=30 measured nonconstant genes for a per-cell primary correlation; report smaller strata descriptively or as not_estimable.
- Gene-centered cell-profile correlation and cell-pair expression-difference agreement, testing specificity beyond a gene's overall abundance.
- Expression/detection/variability strata where observed fields exist. Include low-expression genes in the coverage accounting.
- Generic-PBMC-context and shuffled-cell-description diagnostic baselines on a frozen panel; use donor splits to estimate target reproducibility where donor data permit.

Population-average reference prediction is not donor-level prediction. Personalized donor expression is an optional separate module requiring approved genotypes/haplotypes and expression, validated sequence construction and a within-gene/cell deviation metric. Missing donor inputs never block the aggregate or variant modules.

## 6. Allelic direction and magnitude modules

Analyze all scored observations and separate prespecified evidence strata: discovery-cell effects, supported comparison-cell effects, uncertain comparison-cell effects, and high-PIP statistical candidates. High PIP in one discovery cell is not proof of a nonzero effect in every other cell.

Define reliable observed signs using a preregistered source-supported interval excluding zero; by default use beta +/-1.96 SE only when that normal approximation is valid for the source fit. Sparse rare-variant estimates may not support it; report separately. These are statistical-evidence strata, not experimentally causal labels. Missing SE removes precision stratification, not otherwise valid beta-based metrics. Freeze release-specific flags and significance conventions before analysis.

For each stratum report eligible rows, observed exact zeros, valid scores, correct/incorrect signs, model zeros, missing/failed scores, conditional sign concordance, direction coverage and all-eligible correct fraction. Observed exact-zero betas have no direction and need their own count. Failures are not biological zeros. Explain denominators explicitly.

Use locus-cluster bootstrap and zero-preserving randomizations. Statistical evidence and a utility threshold are separate outputs: **never set a p-value to 1 merely because all-pair success is <=50%.** For a symmetric random-sign null at directional coverage q, expected all-pair success is q/2, not automatically 1/2. Use empirical nulls preserving dependence/coverage and investigate observed sign imbalance; a simple q/2 calculation is an illustration, not a substitute for the null. Report majority-sign baselines learned on development data when relevant. Do not retroactively rewrite historical conclusions.

Magnitude: report both signed Spearman(delta,beta) and Spearman(abs(delta),abs(beta)), with ties/zeros retained. Analyze pooled, per-cell, per-gene-cell (>=20 valid variants and nonconstant values), and equal-gene/locus summaries. Mark undefined correlations not_estimable. Large variant counts in one LD-correlated gene do not establish independent precision. Report beta-precision strata and calibrate scales only on development; do not interpret SAIGE beta and GI log(TPM+1) as inherently interchangeable.

## 7. Cell-specific variant effects

Use every available measured cell for the same variant-gene pair, not a six-cell cap. Compute within-pair cross-cell profile Spearman (>=5 nonconstant measured cells), then summarize equally across pairs/loci. Also calculate all compatible pairwise contrasts C_GI=delta_a-delta_b and C_obs=beta_a-beta_b, with correlation, direction, zero rate and coverage. Summarize within variant-gene before equal-pair/locus pooling; many pairwise contrasts from one locus are not independent observations.

Provide cell-pair heatmaps and related-lineage versus distant-lineage comparisons. Rank cells by abs(delta) versus abs(beta) and distinguish supported affected cells from inconclusive estimates. No significance in a cell is not evidence of no effect. Formal unaffected labels require a prespecified equivalence margin and adequate uncertainty information.

Use covariance/sharedness/heterogeneity estimates only when compatible released inputs support them. With unknown cross-cell covariance, retain descriptive contrasts and the prior conservative |beta_a-beta_b|>2(SE_a+SE_b) sensitivity, explicitly not a calibrated interaction test. Missing covariance skips formal heterogeneity inference, not correlation of available contrasts. A missing cell only removes comparisons involving that cell.

## 8. Fine-mapping and ranking

Do not require controls to score candidates in the core atlas. Fine-mapping is independently gated on available statistical labels and candidate lists.

A. Discrimination: default high-PIP >=0.90 versus explicitly reported low-PIP <=0.01; <=0.10 is a sensitivity arm. Match within locus/gene/cell where possible and adjust/match MAF, TSS distance and local GC. Report eligibility, positive prevalence, AUROC, AUPRC/average precision and locus-bootstrap intervals. Unknown PIP stays outside binary truth. Name this statistical-label discrimination, not experimental causal classification. Never borrow a discovery-cell PIP as the label for another cell.

B. Complete loci/credible sets: score all released competing variants that the model covers, not only lead/high-PIP members. Report within-locus ranks, top-1/top-5 retrieval of prespecified high-confidence candidates and, if component probabilities permit, posterior support captured by top-ranked variants. Marginal PIPs need not sum to one across multiple signals; normalize only when the statistical object justifies it. Treat PIP as evidence of inclusion, not an effect-magnitude label.

C. Coverage: report each credible set's scorable fraction and, when defined, omitted posterior support. Whole-set claims require complete relevant candidates; partial sets receive a separate analysis. Missing component membership skips whole-set analysis, while known-label classification can continue.

D. Controls: low reported PIP, missing PIP, and experimentally assessed negatives are separate categories. Report magnitude wins/ties/losses and tie-aware ranking (standard AUROC tie treatment); do not equate strict wins with AUROC. Existing uncertain-control comparisons remain exploratory. Genuine experimental labels are an optional independent validation module, not assumed available.

## 9. Rare-variant prioritization

Run rare single-variant scoring when a normalized rare-variant list and valid sequence inputs exist; use release-defined frequency/allele-count strata, including singletons where supported. Missing gene-level rare tests does not prevent scoring or available single-variant direction/magnitude analysis.

Gene-cell enrichment additionally requires compatible rare-gene test results and exact or adequate variant-to-tested-set membership. Freeze aggregation on development: mean top-k absolute delta (default k=3, use/report available k) and fraction exceeding a threshold learned on development; maximum is a sensitivity, not the only score. Account for tested variant count, minor-allele counts, covered regulatory span, callability and donor/cell sample size when supplied. Report unavailable covariates and omit only the corresponding adjusted model. Observed-expression adjustment is optional; missing expression must not block the entire rare analysis.

Compare high GI gene scores with rare-variant-associated gene-cell labels using matched enrichment/regression and available baseline annotations, on held-out loci where tuning occurred. Keep null-result genes and untested genes distinct. Match variants to the exact mask, frequency cutoff and gene/test definition; do not attach a gene-level p-value to an incompatible arbitrary collection of variants. If membership is absent, skip set-specific enrichment rather than guess it.

Actual improvement in association power is a separate optional retesting module requiring approved donor genotypes/counts or sufficient score-and-covariance summaries. Compare GI-informed masks/weights with original tests under a fixed error-control procedure and independent validation. Gene-level p-values alone cannot support a rerun. Missing inputs skip retesting, not enrichment or common-variant work.

## 10. Subgroups, baselines and inference discipline

Predefine bins from annotations or data-quality fields, never from favorable GI results: every cell and immune family; each estimable gene/gene-cell group; observed expression/detection/variability where available; biotype/TSS confidence; pathway/promoter annotations; distance to TSS; variant class; frequency/allele count; PIP and credible-set size; discovery versus comparison cell; beta precision; sample/cell counts; numerical-zero class; coverage and geometry.

Unavailable annotations skip their own strata. A gene symbol/ontology lookup failure must not remove an otherwise valid stable-ID score. Report sample sizes, missingness, coverage, conditional/all-eligible metrics, equal-gene/locus and observation-weighted summaries. Use >=30 independent locus clusters for a formal subgroup inference by default; show smaller strata descriptively. Per-gene rankings remain descriptive unless adequate independent signal units support inference. Compare subgroup effects directly; do not infer differences from one significant and one nonsignificant result. Control exploratory subgroup multiplicity (BH FDR by predefined family), use locus-cluster intervals, and label post-hoc discoveries.

Investigate the joint question: do accurately predicted expression profiles coincide with accurate variant effects? This joint analysis needs both inputs; missing expression skips the joint analysis only.

External model comparisons require actual matching executable/licensed scores and compatible cell/track/output mapping. Compare on a shared intersection while separately reporting each model's native coverage. Do not force every architecture to GI's short context; specify native versus matched-context estimands. Missing PromoterAI/AlphaGenome/Borzoi access or incompatible cell tracks skips that comparator only. Published correlations from different cohorts are not head-to-head results.

The descriptive atlas includes historical loci. New calibration/aggregation/model choices use locus-disjoint development and held-out evaluation, frozen before outcomes; unknown historical provenance removes untouched-holdout claims, not atlas rows. Use prior development for engineering where practical. No new unrestricted prompt search. For raw no-fit summaries report all-data results, clearly distinguished from held-out claims. GI training-data independence is not established by an API holdout.

## 11. Execution and reporting

Partition by source family/cell/chromosome; use bounded streaming, explicit disk checks and durable caches. Count unique REF and ALT input-context requests before execution; the owner has not imposed a total call cap. No arbitrary sample truncation. A run interrupted by runtime/storage/API limits is partial with an exact remaining-work ledger, never full coverage. Resume verified completed partitions without replay.

Pin resolved sources, model/API/input version, context registry, normalized universe, split rules and metric definitions before evaluating new choices. Default bootstrap 3,000 locus replicates; randomization 2,000, with documented Monte Carlo resolution. Expensive analyses may be partitioned but not silently reduced. Donor bootstrap is additionally supported for expression when authorized input exists. Keep computational failures separate from planned unavailable-data skips.

Always produce module_status.tsv, availability.json, coverage/exclusions, metrics for available modules, and report.md with one section per benchmark (including explicit SKIPPED sections). Output design and status contracts are in ATLAS_AVAILABILITY.md. Save under results/atlas-<UTC timestamp>/ and export sanitized summaries to docs/runs/atlas-<date>/. No raw human data, DNA or credentials in this public repository. A successful execution may have completed_with_skips status; it is not evidence that every goal was answered.

## Sources and provenance

- Existing measured evidence and limitations: [independent follow-up summary](runs/followup-20261030/summary.md), [implementation](FOLLOWUP_IMPLEMENTATION.md), and [metrics](runs/followup-20261030/metrics.json). These remain unchanged.
- Public association inventory: https://zenodo.org/records/18221260 ; verify live source metadata when executing.
- Expression access lead only: https://ega-archive.org/datasets/EGAD50000002379 ; verify availability and permission before use, otherwise skip.
- GI contract: https://docs.genomicintelligence.ai/rest-api ; use the audited official schema/client, not guessed endpoints.
- This protocol translates the project owner's comprehensive-study plan and nonblocking missing-data requirement. New thresholds and design choices are proposed specifications, not new measured findings.
