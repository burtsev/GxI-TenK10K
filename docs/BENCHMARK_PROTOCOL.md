# Prespecified scientific benchmark protocol

## Hypotheses
H1: Zero-shot GI predicted variant effects (ALT minus REF in log(TPM+1)) have better-than-chance sign concordance with TenK10K oriented eQTL betas in comparable immune cell types.
H2: GI can predict differential allelic effects across cell types.
H3 (exploratory): sequence-based prioritization scores enrich high-PIP variants relative to matched lower-PIP variants within a locus.

## Files and version
Primary record: https://zenodo.org/records/18221260. Query Zenodo REST API `https://zenodo.org/api/records/18221260` for authoritative filenames, sizes, links, checksums, licence and version. First download `common_variant_gene_level_results_with_annotated_fails.zip` (~10 MB) and `susie_summary.zip` (~60 MB), *if inventory confirms these exact names*. Download `common_all_cis_pvalues_100kb.zip` (~16.6 GB) only once schemas and disk are verified; stream extract to a single file at a time. Optional `rare_variant_gene_level_results.zip` after pilot; avoid ±1 Mb and 23.5 GB rare-single-variant archives by default. Document if observed inventory differs.

Published TenK10K summaries cover approximately 1,925 WGS-matched individuals and 28 immune cell types; verify paper/metadata rather than infer per-row counts. Dataset is NOT raw donor-level expression. Check assay protocol, genome build, phenotype-normalization, variant IDs, association software and beta interpretation from accompanying publications/project code.

## Candidate definition
- Pilot target: 100 **distinct** high-confidence candidate variant–gene pairs from different genes/loci where possible, selected reproducibly from SuSiE high-PIP credible sets combined with significant common-variant eQTL summary data. Prefer PIP >=0.9 if sufficient; relax to >=0.5 and label if needed. Primary analysis requires orientable non-palindromic biallelic SNVs (or explicitly strand-resolved ambiguous alleles), non-ambiguous gene/TSS mapping, variant within actual model scoring interval, and unambiguous effect allele.
- Spread across at least B-naive, CD4-naive T, CD14 classical monocyte and CD16 nonclassical monocyte if available; map to *exact* TenK10K labels from metadata. For each pair, score the discovery cell and up to three comparison populations with comparable measured betas. Pre-specify comparison cell selection by lineage, measurable gene expression, observation power (when available), and matched assay; do not select the comparison by GI output.
- Controls: 100 allele-frequency, TSS-distance, local-GC and locus-matched variants with clearly defined control category. An insignificant eQTL is an uncertain negative; use PIP contrasts for ranking and label proxy interpretation. All controls must have allele-oriented eQTL summary estimates when used for beta metric.
- Initial 2-variant smoke: CSF3R rs3917929 at GRCh38 chr1:36,478,697 (check allele IDs and beta in actual data) and CETP rs2052880 chr16:56,971,389 as a second candidate; verify alleles in file before editing. A second CETP variant at chr16:57,020,104 is potentially out of effective GI range: record as coverage test, *not* primary efficacy test.
- Do not use initial IL7R qualitative model outputs as scientific ground truth.

## Data processing and allele orientation
Read actual columns with explicit data dictionary. Standardize coordinates (GRCh38, one-based variant POS, zero-based BED intervals), contig naming, indel normalization and variant identifier. Variant annotation and TSSs should be based on consistent GENCODE/Ensembl releases. Handle multiple TSS via a deterministic canonical rule and record exclusions or sensitivity analysis. Minus-strand sequence oriented transcriptionally; if genomic REF/ALT are substituted in + strand before reverse complement, effect remains ALT minus REF and does **not** change sign simply due to gene strand. Check if published beta uses ALT, effect allele or reference allele and flip observed beta only when necessary; preserve provenance.

For each gene fetch exactly the model's documented input window around correct oriented TSS from GRCh38 (FASTA recommended via Ensembl REST or downloadable indexed GRCh38), record ref sequence SHA256. Implement `apply_snv` that checks the expected one-based genomic REF, edits exactly one nucleotide, preserves sequence length, validates alt base and reverse complement symmetry. Reject variants outside the window or the actual scored interval. If model scores only a central subwindow, calculate variant coverage in the same coordinate frame and verify with response metadata. If coverage metadata unavailable, mark as unknown and analyze separately.

## GI inference
Use `g0-expression-8192` as primary if live model catalog confirms; optional `g0-expression` comparison with its own independent window and context. Query docs/OpenAPI rather than assume these model IDs never change. Model conditions: replicate biological cell type, assay (10x scRNA-seq), human peripheral blood, unstimulated state, and relevant annotation; **same context string for REF and ALT**. Store exact requested text. Score reference once per gene/cell/model; alternative once per variant/cell/model. Cache by full content hash, record API response units and errors. Primary delta is `pred_alt_log_tpm1 - pred_ref_log_tpm1`. Never compare raw TPM deltas to logged β without explicit rationale.

Previously observed G0 API tests returned nominal 81,920-bp inputs with effective scored windows smaller than 81,920 bp; validate actual coverage from metadata, not submitted input length.

## Evaluation and statistics
1. Main metric: observed allele-oriented beta sign vs GI delta sign, ignoring a pre-specified near-zero delta threshold (report threshold and abstentions). Report binomial/cluster-bootstrap 95% CI, n variant-gene, n loci, n cells, and a permutation test by locus.
2. Secondary: Spearman rho and robust regression of GI delta vs beta, stratify by study-assay and cell type; note SAIGE-QTL beta vs log(TPM+1) unit mismatch.
3. Cell specificity: compare GI `delta_cell_A - delta_cell_B` with similarly oriented observed beta contrast; calculate Spearman, sign and uncertainty, with careful recognition of differing precision and cell-type power. Do not call differences significant without standard errors and model.
4. Fine-mapping: within independent credible sets, evaluate ranking of PIP-high vs PIP-low candidates and PIP-weighted metrics; avoid treating PIP as certain ground truth. Use matched random and distance-to-TSS baseline.
5. Baselines: chance (50% sign), distance/TSS (can be useful for magnitude/priority but has no direct allelic sign), local allele-frequency and sequence heuristics. Optionally established models if software accessible; say explicitly if skipped.
6. Bias/robustness: restrict to high-confidence single-variant credible sets; evaluate rare variants separately; stratify by expression abundance, beta magnitude, TSS distance, cell power, SNP class, and possible LD. De-duplicate by locus across splits; bootstrap clusters by locus or gene.
7. Leakage: G0 training includes human GTEx/ENCODE/cellxgene RNA data; a TenK10K holdout is not formally guaranteed. Document if overlapping donors or training data cannot be excluded.

## Go/no-go thresholds (descriptive targets, not performance claims)
- End-to-end live smoke produces both REF/ALT scores, correct strand and observed effect-allele match on >=2 variants.
- Pilot produces >=100 valid candidate pairs, >=100 control pairs, >=4 study cell types and valid denominators if possible; otherwise quantify shortfall.
- Primary sign concordance >0.50 is suggestive, not sufficient by itself; seek locus-bootstrap lower 95% CI above 0.50 and improvement over pre-registered baselines. No inferred success without results.
- All filtered/unscored cases must appear in an exclusions table with named reasons.

## Deliverables
`manifest.json`, `data_audit.md`, `candidates.tsv`, `controls.tsv`, `requests.jsonl` (redacted), `predictions.parquet`, `observations.parquet`, `exclusions.tsv`, `metrics.json`, `report.md`, `figures/*.png`, `environment.json`. Use schema versioning and include variants, model ids, cell IDs, coordinates, observed beta/SE/PIP, effect allele, ref/alt log scores/delta, model effective coverage, request success status, and provenance. Keep all large or sensitive artifacts out of Git.
