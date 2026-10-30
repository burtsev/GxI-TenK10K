# GI × TenK10K follow-up — frozen scientific protocol

**Status:** proposed and not yet executed. This document is the authoritative specification for the new study. Codex must preregister its resolved model, dataset, candidate lists, endpoints, exclusion thresholds, prompt matrix and code revision **before** scoring any confirmatory locus. All deviations require a dated amendment in `amendments.md`, with a transparent separation of affected exploratory results.

## 0. Motivation and evidence, not new results

Read `docs/FOLLOWUP_PRIOR_FINDINGS.md` for the detailed historical results and caveats. The prior 100 candidate / 100 matched-control pilot found **52/77 = 67.5%** correct effect direction among non-abstaining candidates (23/100 near-zero), Spearman rho **0.370** across all 100, candidate/control non-tied magnitude ranking **60/76**, and **32/71 = 45.1%** correct cross-cell effect-contrast directions. It used a post-start protocol amendment to **40,960-bp** windows and a control-matching amendment, so it was exploratory. Text optimization improved in-sample development metrics but yielded no convincing advantage on held-out 20-locus cohorts (latest baseline vs selected: 10/20 vs 10/20). This follow-up must NOT treat any of these as independent confirmation. Historical reports are referenced in Slack #general (October 29–30, 2026); repo currently may not contain their underlying artifacts.

## 1. Questions and prespecified outcomes

**Primary scientific question:** can GI predict whether the allelic effect of the *same variant on the same gene* differs between two immune cell types?

- **Co-primary A (cell-effect contrasts):** across independently selected held-out locus clusters, calculate observed `C_obs = beta_Allele2_ALT(cell_a) - beta_Allele2_ALT(cell_b)` and predicted `C_GI = delta_GI(cell_a) - delta_GI(cell_b)`, where `delta_GI = expression_log_tpm1(ALT) - expression_log_tpm1(REF)`. Primary statistics: *all-pair correct-contrast direction proportion* (zero/tied model effects count as incorrect, but report abstentions separately) and Spearman rho including model zeros. Estimate 95% locus-cluster bootstrap CIs, within-locus permutation null, and total eligible / scored / informative counts. Specifically assess whether the locus-bootstrap lower CI for rho exceeds zero; treat this as stronger evidence than small-sample direction alone.
- **Co-primary B (variant effect direction):** all-pair correct-direction proportion among held-out high-PIP variants, with model zeros counted as incorrect; additionally report conditional direction on non-abstainers. Compare to oriented observed BETA, clustered 95% CIs and a locus sign-flip null. No claiming 67.5% without stating 77/100 and abstentions.
- Both are co-primary, so provide multiplicity adjustment (e.g. Holm on the two preregistered primary p-values), effects, CIs and decision rule. Exploratory analyses must not change primary selection.

**Secondary:** relative change vs the fixed baseline under *held-out* context changes; candidate-versus-matched-control |delta| ranking and PIP-aware ranking (controls are uncertain, not true negatives); variation by cell type, MAF, PIP, TSS distance, expression level; zero-effect prevalence; observed-to-predicted effect-scale calibration; ablation of cell and gene labels; sequence-window geometry; and shared-model baseline. **Optional** independent caQTL corroboration and PromoterAI matched comparison, treated as separate analyses.

## 2. Sources and exact compatibility

- TenK10K Phase 1 Zenodo **record 18221260**, https://zenodo.org/records/18221260 (v2; discover live metadata/checksums): `common_variant_gene_level_results_with_annotated_fails.zip`, `susie_summary.zip`, `common_all_cis_pvalues_100kb.zip`. The original report scanned 63,986,931 common full-summary rows for four populations; new work should cover **as many of the 28 populations as data/comparability permit**. Use `results_iscovariateoffset_true.zip` only with explicit protocol provenance; don't silently merge different model fits.
- Annotate GRCh38 primary assembly with GENCODE v44. Verify official reference checksums and document transcript/TSS choice. Biallelic SNVs only for confirmatory analysis; indels/rare variants separately.
- SAIGE-QTL `Allele2` is historically described as the effect allele: recheck published file metadata/source code, normalize to **ALT effect** by explicit allele mapping, record any flip. One-based genomic variant position vs zero-based intervals, plus/minus strands, REF allele identity, TSS strand orientation must be tested and logged. Report missing BETA/SE/MAF and convergence errors.
- Original study used `g0-expression-8192`, with `model` pinned explicitly on each call, 40,960-bp genomic sequence in **gene-sense orientation** with **0-based `tss_index=20480`** and same input geometry for REF/ALT. Treat this as frozen primary geometry **after verifying live API model response**; fail preflight or amend *before* holding-out test if actual model no longer supports it. The 81,920-bp input formerly gave 4–7-bp REF/ALT differences in effective scored boundaries and even an example sign reversal; include 81,920 as a separately reported geometry sensitivity only. For confirmatory input, variant must be inside **both** returned effective scored intervals with identical genomic coverage. Store returned intervals.
- Response unit: predicted `log(TPM+1)`; observed SAIGE-QTL BETA may be on a different statistical scale. Avoid reporting direct calibrated effect magnitude without an explicitly fitted, held-out evaluated calibration.
- Capture version of GI model, REST contract, exact context text, genomic build, source checksums and code SHA. GI training data may overlap other datasets; no claim of proven independent *training* holdout.

## 3. Holdout design and power targets

The key is **locus isolation from everything previously inspected**. All original 100 candidates, their controls, two smoke and 20-trial loci, 16-pair description-development loci, 20 first prompt-validation loci, 20 later prompt-validation loci, every repeated-paraphrase locus, and any other previously scored locus must be grouped into the historical exclusion list, using *any* model/context and both alleles. Recover exact locus IDs from the prior VM's local results, request/cache ledger and scripts. Previous final-report GitHub commit `b7ecbed` was mentioned in Slack but was not found on `main` during inspection. Do not assume prior artifacts are already pushed. First recover and commit code/aggregate results to a separate archival branch if available. **If the historical scored-locus set cannot be recovered, do not claim a truly untouched confirmatory holdout**: continue the pipeline and clearly mark the study exploratory or blocked on holdout independence.

Cluster all candidate gene TSSs by overlapping +/-100-kb windows (200-kb centers separation as a first approximation), and if feasible merge by actual cis-LD blocks from a legally accessible matched reference panel; never call physical distance proof of independence. Also place excluded historical loci in the same clusters. Split entire locus clusters into (A) fresh **development** pool, (B) final untouched **test** pool. No variant or gene can cross pools through overlapping windows or LD block. Freeze both lists and their hashes before *any* new GI predictions, with no outcome-driven reassignment after inspection.

**Targets, contingent on eligible records:**
- Fresh development: **200–300 independent loci**, stratified across cell families, used ONLY to choose text/context ablations and debug model sensitivity.
- Confirmatory held-out: **at least 1,000 distinct high-PIP variant–gene pairs from >=600 independent locus clusters** if feasible; expand to 2,000 if available. PIP >= 0.9 in at least one cell is preferred. If not enough, report available high-PIP sample and separately label a PIP 0.5–0.9 secondary tier; do not silently relax a preregistered primary cutoff.
- Paired uncertain controls: target **one control per candidate** where matched (same gene/cell, |MAF gap| <=0.05, |absolute TSS-distance gap| <=5 kb, 201-bp GC gap <=0.10). PIP missing vs low reported separately; matched eligibility is frozen before predictions.
- At least **12 immune populations**, ideally all 28; report exact metadata mapping and per-population power. At least 4 main families (B, CD4/CD8 T, NK, monocyte) if data permit. For each variant–gene pair use up to **six** cell types (discovery plus prespecified related/distant types), where per-cell observed BETA/SE exist, avoiding cherry-picking based on GI predictions. Ensure a substantive number of cell contrasts (aim >=500 informative, >=300 independent loci).
- Proportion of effects sufficiently far from zero and high information content should be reported as a selection stratum, not silently made the only results. If sample targets prove unavailable, report attained n and do not fabricate statistical power.

**Uncertainty-aware contrast subset:** primary includes *all preselected cross-cell comparisons with complete betas*, including weak contrasts and model-zero effects. Additionally prespecify a conservative "informative observed contrast" subset where `|beta_a-beta_b| > 2*(SE_a+SE_b)` (upper bound for standard error of a difference when covariance unknown), with BETA and SE on the same scale. This is a descriptive high-confidence stratum, **not a formal cross-cell significance test**. Also analyze expected precision and unequal cell counts. If full covariance/sharedness results become available, use a more defensible interaction or sharedness statistic and retain the initial subset as preregistered.

## 4. Fixed reference context and text-conditioning experiments

**Primary baseline:** preserve the original biological description template exactly, using verified cell identity:
`Human peripheral blood {cell}, unstimulated, 10x 3-prime single-cell RNA-seq`
The phrase "unstimulated" and chemistry must be validated against TenK10K metadata; if unverifiable, preregister corrected text and explain that the prior baseline cannot be reproduced verbatim. Use *one* exact string per cell; REF/ALT must always use identical text for a given comparison. Do not use a different cell-context template for each locus.

Evaluate **prespecified prompt arms** on the *development pool*, with identical DNA/variants/model/TSS:
1. Original baseline text.
2. Cell identity only (canonical human-readable name).
3. Canonical cell identity + verified Cell Ontology (`CL:...`) ID.
4. Exact Cell Ontology ID only.
5. Baseline plus structured assay/organism/biosample fields; confirmed values only.
6. Baseline plus gene symbol.
7. Baseline plus Ensembl stable gene ID (without gene symbol).
8. Baseline plus both gene symbol and Ensembl ID.
9. The previously selected exact repeated paraphrase (historical comparison, not presumed best).
10. Simple deliberately duplicated baseline (mechanistic repetition control).

Ablation contrasts should isolate one semantic field at a time; for example baseline vs baseline+symbol vs baseline+Ensembl, and canonical label vs canonical+ontology. Do not misrepresent a composite of several changed fields as a one-factor ablation. Have a fixed `description_registry.tsv` with exact strings and SHA256 before inference. If a population lacks validated Cell Ontology mapping, mark the ontology arm ineligible; never invent IDs.

**Leakage/conditioning diagnostics:** use nonsense gene identifiers, deliberately swapped gene symbol / Ensembl ID (clearly nonbiological negative controls), and shuffled cell labels on a development-only diagnostic subset to test if predicted deltas rely on text lookup rather than DNA. These are *model-behaviour* tests, not valid biological inputs for a benchmark. Try exact repeated identical requests to measure numerical determinism; no cherry-picking based on API noise.

**Selection of one frozen alternative**: on dev only, maximize **all-pair correct cross-cell contrast directions** as the first criterion (with minimum effect coverage and no selective abstention), followed by dev contrast rho and then all-pair variant sign. Use one-sided/repeated-arm corrections and report the entire candidate family. Only freeze an alternative if it improves the *preregistered development primary endpoint* by at least 3 percentage points without reducing eligible coverage by more than 3 percentage points; otherwise freeze baseline as the "winner". The 3 pp threshold is a **decision heuristic**, not evidence of statistical significance. Lock selection script and string registry before confirmatory scoring. Evaluate baseline and **one** frozen alternative on the same untouched test variants, using paired locus-bootstrap CI. Do not retune after seeing test outcomes.

Important: adding gene symbol/ID risks learned-gene lookup and should be reported as a **separate gene-aware conditioning track**, not replace the core DNA+cell-context zero-shot result. Use the same cohort in its matched comparison but label that result separately.

## 5. Computational workflow

Implement deterministic, streaming processing and a GI API client based on the official published REST contract and/or pinned integration-kit client. Use `GI_API_KEY` from environment, `GI_BASE_URL=https://api.genomicintelligence.ai`. Verify model catalog endpoint and schema at runtime. Expression inference uses documented `tss_index` and `options={"description": text}`; pin model ID explicitly. Save only sanitized request/response metadata, sequences by hashes. Reuse REF scores per gene+cell+description+model+geometry. Cache successes/failures with proper provenance; keep raw/reference DNA on the VM, not Git. Separate per-inference output from evaluation outcomes.

**This key has no quota cap** per study owner. Set `api.enforce_total_quota=false`, `api.max_total_requests=null`; never reintroduce previous ceilings of 790, 800 or 1,600. The model endpoint may still impose server rate limits, concurrency constraints, network timeouts and monetary usage: honour HTTP 429 `Retry-After`, bounded retries and exponential backoff; discover rate limits rather than claiming infinite throughput, track usage, provide optional user-configurable cost guard but don't block at arbitrary old attempt count. Do not manufacture an API price. Run with adaptive concurrency (start 2, increase conservatively if accepted; no uncontrolled flood), batch grouping for cache efficiency and robust checkpoint resume. Stop only on scientific/data integrity failures, inaccessible APIs, sustained errors, explicit user constraints or operational safety conditions; report partial results honestly.

Geometry QA: compare 40,960-bp REF/ALT effective scored spans and allele match on every call. Separately run preselected 81,920-bp sensitivity (e.g. 100–200 distinct loci), same TSS/cell prompts. Never pool geometry arms in confirmatory primary.

## 6. Analysis: prespecified estimands and controls

**Variant sign:** `sign(delta_GI)` vs `sign(beta_ALT)`. Primary denominator all valid candidate variant-cell observations; predicted |delta| <=1e-6 is incorrect/abstain, not discarded. Secondary conditional accuracy excludes abstentions and exposes eligible n. Report additional epsilons 1e-4, 1e-3 as sensitivity only. Count full missing/error denominator separately.

**Cross-cell effect:** predefine cell pairs for each variant before GI inference. `C_GI=delta_A-delta_B` and `C_obs=beta_A-beta_B`; primary all-pair correct sign and rho with zeros; report pairwise n, 95% locus CI, permutations shuffling cell descriptions *within each locus* as a null for context information. Gene/cell-specific reference expression checks are descriptive, not valid eQTL variant effects.

**Calibration:** observed BETA is not directly the GI output scale. Fit any slope/intercept or monotone calibration on development only, evaluate on frozen test; raw rank and direction remain primary. Report per-cell/precision-weighted sensitivity if SE available but not infer cross-cell statistical significance without covariance.

**Matched controls:** rank |delta| candidate > control, ties not hidden, plus PIP-stratified within-locus ranking. Controls with absent PIP are *unknown*, not PIP=0 or causal negatives. Control matching must be pre-inference and blinded to GI predictions.

**Baselines:** sign=chance 50% for non-abstaining; permuted cell labels and permuted allele effect signs (respect locus/cell dependence) to quantify context/allele information; TSS-distance and MAF/GC heuristics only for effect magnitude and prioritization (they cannot predict allelic sign without additional allele features). For a cell-agnostic baseline, its predicted cross-cell contrast is zero; report it as zero-coverage instead of assigning artificial 50% direction success. Compare to published PromoterAI 0.66 only narratively unless **re-executing PromoterAI on exactly the same allele/cell pairs**; the published median correlation is not same metric. PromoterAI precomputed scores or trained model may require a distinct **commercial licence**; check and obtain legitimate access before running. Absence is not a reason to block GI primary analysis.

**Statistics:** bootstrap >=2,000 iterations by independent locus clusters; predeclared 95% confidence intervals; clustered permutation; 2 co-primary corrected p-values; within-cell type exploratory performance with complete denominators, counts and CIs; compare model variants by paired locus-bootstrap not disconnected proportions; do not hide abstentions or failed inputs. Document selection bias and unknown GI training overlap. Do not interpret observational fine-mapping as experimentally causal.

## 7. Required deliverables and report format

In `results/followup-<timestamp>/` produce:
- `source_manifest.json` (download URLs, checksums, GENCODE/FASTA IDs); `legacy_exclusions.tsv`; `cohort_flow.tsv`; `cohort_development.tsv`; `cohort_test.tsv`; `cohort_hashes.json`; `frozen_protocol.json`; `amendments.md`.
- `description_registry.tsv`; `cell_mapping.tsv` with exact TenK10K labels and verified ontology; `allele_checks.tsv`; `coverage_exclusions.tsv`; `requests.jsonl` (redacted metadata); `predictions.parquet`; `observed_effects.parquet`; `metrics.json`; `figures/`; `report.md`; `reproduce.sh`; `environment.json`.
- A study comparison table: legacy exploratory report, frozen baseline independent test, frozen alternative independent test, and matched controls, **each with explicit counts, origin, uncertainty, and no unsupported cross-experiment relative claims**.
- Figures: cross-cell GI vs observed contrast scatter and per-cell heat map; all-pair versus conditional accuracy / abstention; calibration and null distributions; prompt-ablation effect sizes with CIs; effective-window/geometry diagnostics; per-family coverage and power.
- Create small sanitized, reviewable protocol and summary files in Git; never push private data, API keys, large archives or unnecessary raw DNA. Data and full results stay on VM; explicit manifest permits reproduction.
- Document exact CLI commands, source versions, git SHA, model version, API/schema version, request counts, rate-limit/timeout events, full result denominators and all deviations.

## 8. Stop and report conditions

Never call an exploratory historical cohort a holdout. If previous results cannot be recovered, correct legacy exclusions first, or execute exploratory work without confirmatory claims. If the GI API changes sequence geometry/schema or a source's allele orientation is unresolved, pause that analysis and fix data QA rather than computing meaningless metrics. Don't adjust protocol on test results. If exact sample targets are unattainable, report attained n and CI rather than inventing results. All optional external baselines must be explicitly marked run/not run and why.
