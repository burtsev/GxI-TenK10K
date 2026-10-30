# Independent assessment of DNA-to-expression variant effects and immune-cell specificity with Genomic Intelligence

Computational research report | 30 October 2026 | Manuscript for project review

## Abstract

We assessed GI predictions against TenK10K summary effects using 250 development loci and a frozen test of 1,000 variant–gene pairs across 877 genes and 777 physical locus clusters, with matched uncertain controls and 28 immune populations. Previously scored loci were excluded using the recovered request registry. No gene-free description alternative qualified in development, so baseline alone was tested. Variant direction was correct in 2,392/5,863 observations: 40.8% (95% locus CI 38.6%–43.0%) with zeros counted incorrect, and 57.9% among 4,131 nonzero predictions. Cross-cell direction was correct in 6,195/14,446 contrasts: 42.9% (CI 41.4%–44.4%), or 52.2% among 11,860 nonzero contrasts. Cross-cell Spearman rho was 0.061 (CI 0.030–0.091). Raw clustered randomizations and the published all-pair decision gate are distinguished below. The assessment does not establish independence from GI training data or provide a matched head-to-head comparison with PromoterAI.

## Introduction

A useful regulatory sequence model should predict both the sign of an allelic effect and how that effect differs between cell types. TenK10K provides summary-level estimates for this assessment. The Nature version of record [1] supersedes the linked medRxiv preprint [2]; they are versions of one study. Our previous pilot and text experiments were exploratory and motivated this independent cohort rather than serving as confirmatory evidence.

## Methods: data and isolation

We used only Zenodo record 18221260 [3] for association summaries, official GRCh38 primary sequence and GENCODE v44 annotations [4]. Three archive checksums and SHA256 values were verified; all 28 full-summary ZIP members were streamed with schema and CRC checks, without downloading donor-level data. Canonical annotation required an unambiguous TSS among tagged canonical/MANE transcripts. Candidate SNVs required PIP ≥0.90 in a significant discovery gene–cell association (per-cell BH q≤0.05), MAF≥0.01, complete finite effect/SE and verified genomic REF. We recovered 2,783 cached historical request keys plus one uncached positive control, spanning 167 genes. Historic and eligible prospective TSS ±100-kb windows were merged transitively, and entire historical-containing clusters excluded. Development and test share no cluster or gene. These physical clusters approximate dependence; matched LD blocks were unavailable.

## Methods: matched inputs and GI inference

Each candidate received one fixed uncertain control in the same gene, with MAF difference≤0.05 in each selected cell, absolute TSS-distance difference≤5 kb and 201-bp reference GC difference≤0.10. Reported control PIP had to be below 0.5; missing PIP stayed unknown. Discovery plus related/distant cells (up to six) and all cell contrasts were selected before inference. SAIGE-QTL Allele2=ALT determined observed BETA orientation. SNVs were edited on genomic-plus sequence before reverse complementation for minus-strand genes. GI g0-expression-8192 [5] used 40,960 bases and oriented TSS index 20,480. REF/ALT text and effective scoring spans had to match, with the variant covered. Predicted effect was ALT−REF in log(TPM+1). A distinct cache binds exact sequence hashes, model, description, API contract, geometry and options; bounded retries and adaptive pacing preserve resumability without a total-call limit. The 3-prime assay is confirmed in the Nature methods; unstimulated refers to the absence of a described experimental stimulation arm, not absence of endogenous activation. The catalog supplies an ID, but no immutable model-weight digest; silent server updates cannot be ruled out.

## Methods: context selection and statistics

All ten description arms were evaluated on development only. Gene-symbol and Ensembl arms were separated as gene-aware exploratory conditioning. A gene-free alternative required ≥3 percentage-point gain in all-pair cross-cell direction and ≤3 percentage-point coverage loss; ties used contrast rho then variant direction. Otherwise baseline was retained. The exact registry, source/cohort hashes, analysis plan and implementation hashes were locked before test inference. Zeros (|delta|≤10^-6) and unscored rows remain in all-pair direction denominators; rho retains numerical zeros and reports missing coverage. Co-primary and development-family intervals used 3,000 physical-locus bootstrap resamples; exploratory per-cell, family and covariate/pathway strata used 2,000. Clustered randomizations used 2,000 replicates. Co-primary A uses the larger of direction and rho cell-label permutation p-values; B uses a locus effect-sign null. Holm adjusts A/B; the prespecified direction decision additionally requires >50% all-pair accuracy; below this floor the co-primary p-value is conservatively set to1 before Holm. The raw zero-preserving randomization p-values remain in metrics.json. The informative contrast stratum |beta_a−beta_b|>2(SE_a+SE_b) is descriptive because cross-cell covariance is unknown. Affine scale calibration was fitted only on development.

## Results: independent test

Acceptance checks passed 27/27, with 2,500 reference variants verified and 0 failures; 68 passed in 7.68s. The frozen test comprised 1,000 candidates across 877 genes, the same number of uncertain controls, 777 clusters and 11,726 candidate/control variant–cell rows. The available matched reserve did not reach the 2,000-pair expansion target; no PIP relaxation was used. Baseline variant effects: 2,392/5,863 = 40.8% (95% locus CI 38.6%–43.0%); conditional 57.9% among 4,131; zero 1,732, missing 0; rho 0.183 (CI 0.137–0.230). Baseline cross-cell contrasts: 6,195/14,446 = 42.9% (95% locus CI 41.4%–44.4%); conditional 52.2% among 11,860; zero 2,586, missing 0; rho 0.061 (CI 0.030–0.091). The preregistered all-pair gate was not crossed, giving Holm-adjusted gated decision values A=1 and B=1. Raw randomization evidence is distinguished from that gate below. The cross-cell rho CI lower bound exceeded zero. The descriptive precision stratum yielded 1,458/3,333 = 43.7% (95% locus CI 41.0%–46.6%); conditional 53.6% among 2,722; zero 611, missing 0; rho 0.091 (CI 0.032–0.152).

## Results: frozen description comparison

No gene-free development arm satisfied the frozen gain/coverage rule; baseline was retained. Consequently only baseline was scored on test, with no second description chosen after viewing test outcomes.

## Results: controls and diagnostics

Uncertain controls under baseline: 1,879/5,863 = 32.0% (95% locus CI 30.3%–33.8%); conditional 52.2% among 3,601; zero 2,262, missing 0; rho 0.064 (CI 0.025–0.105). Paired ranking used 5,863 scored comparisons; control PIP was unreported in 5,601 comparisons. These are prioritization summaries, not causal classification accuracy. Six input-only Reactome98 immune regulatory strata, fixed before test inference, are reported with locus intervals in secondary_metrics.json; no pathway-specific significance is claimed. The full development arm matrix with locus intervals and nine-arm Holm comparisons is supplied in development_inference.json and Figure2A. Swapped gene identifiers, a nonsense gene identifier and shuffled cell labels were tested on a fixed 25-locus development subset; their complete denominators and paired intervals are supplied separately. Exact historical-locus repeat checks agreed in 2/2 cases. Larger geometry diagnostics: 150/150 prespecified development loci evaluated; 106 had identical covered REF/ALT spans. Invalid long-window deltas are shown only as geometry diagnostics and excluded from biological sensitivity estimates.

## Secondary matched ranking uncertainty

Candidate |delta| exceeded its matched uncertain control in 2,752/5,863 comparisons, with 1,683 ties and 1,428 losses. Strict all-pair success was 46.9% (95% locus CI 44.7%–49.3%); conditional on a non-tie, it was 65.8% among 4,180 comparisons (CI 63.5%–68.3%). This is a prioritization proxy, not causal classification. Continuous reported-PIP ordering was assessable only where both probabilities were released; its full counts and locus intervals are supplied in reporting_supplement.json. Missing PIP was never assigned zero. Both PIPs were reported for 178 comparisons across 111 clusters. Predicted magnitude ordering agreed in 86/178: 48.3% (95% locus CI 40.4%–56.1%), with 46 predicted ties. This compares released posterior probabilities, not experimental causal labels.

## Discussion: relation to the Nature study

Cuomo et al. [1] mapped genetic effects in 28 immune populations and reported substantial cell specificity. Their PromoterAI analysis reported median correlation 0.66 between promoter scores and lead-eQTL effects. PromoterAI uses a specialized promoter-variant task and a different selected variant population. Our pooled Spearman estimates, high-PIP ascertainment, zeros, cell contrasts and input geometry differ; subtracting these correlations or claiming model superiority would be invalid. No licensed matched PromoterAI predictions were supplied, and AlphaGenome/Borzoi were not executed on this cohort. The comparison is narrative, not head-to-head. The broader implication is that variant-effect agreement and prediction of cell-specific differences require separate evidence; a favorable conditional sign score alone cannot establish the latter.

## Discussion: signal, abstention and the decision gate

The distinction between statistical signal and the prespecified utility gate is material. Variant-effect rho was 0.183 (95% locus CI 0.137–0.230); the cross-cell estimate and its interval are reported above. Raw zero-preserving randomization p-values were 0.0004998 for variant direction and 0.0009995 for the joint cell direction/rho statistic. The gated co-primary decision values were 1.0 and 1.0, using the rule published before test inference. Small raw null p-values indicate orientation information even when the all-pair gate fails; a gate failure must not be interpreted as proof that predictions contain no information. Zero predictions accounted for 1,732/5,863 variant observations and 2,586/14,446 cell contrasts. The conditional percentages expose these abstentions rather than replacing the primary denominator. None of the other development descriptions exceeded baseline's point cross-cell accuracy. This is not a held-out comparison of the unselected alternatives, because baseline alone qualified for test.

## Discussion: limitations

High-PIP significant-gene selection enriches detectable associations and limits generalization to all variants, distal regulation, indels and rare alleles. Requiring both candidate and control estimates in several cells also favors expressed genes. Observational fine-mapping is not experimental causality. Cell labels are approximations, and exact ontology mappings were unavailable for some populations. Physical locus clusters do not establish LD independence; unequal cell sample sizes and unknown cross-cell covariance complicate precision comparisons. GI training overlap is unknown, and development prompt selection is target-informed even though frozen test genes are unseen by these API experiments. API output quantization can create ties. SAIGE Poisson coefficients and GI log(TPM+1) changes have different scales; empirical calibration does not repair all biological differences. Six immune regulatory pathway strata used fixed Reactome98 membership and exact GENCODE44 symbols; overlapping sets and low gene counts make these descriptive. Independent caQTL replication was not run because no matching authorized source was included. All absent analyses are explicitly marked in secondary_metrics.json. Observed gene-expression means were unavailable in the permitted association summaries, so observed-expression-level stratification was not run and no donor-level expression was inferred.

## Conclusion

The two prespecified co-primary criteria were not jointly met. At the 0.05 Holm threshold, the cross-cell criterion did not pass, and the variant-direction criterion did not pass. All-pair and conditional results must be distinguished because abstentions affect the primary denominator. The findings apply to this frozen summary-statistic cohort and API model identifier; they do not establish experimentally causal effects or independence from model training data.

## Reproducibility and availability

Code and sanitized aggregates: https://github.com/burtsev/GxI-TenK10K . Full local evidence: results/followup-20261030. Preflight revision: f57f31505e629232dd602a8fa4f76586200b32c5; pre-test selection code revision: a45d6884b55e24b01e849e7d7c1e1b8bcbcd65f4. Exact pre-test implementation hashes are in code_manifest.json; GitHub publication metadata is in published_freeze.json. Model g0-expression-8192, API2026.10.19.2 (afabae26), contract SHA256 c4c388f0c48bb9bdfd79b2f0c94fc08e04493bff4589ce5ca387fd386b9dc91a; GRCh38/GENCODE44; seed 20261030. Recorded API attempts: 59,688; statuses {'200': 59682}; attempts without a recorded response 6 (interruption/ambiguous dispatches, not assumed successes). Raw archives, DNA and response caches remain local. Reproduce using the pinned requirements and followup CLI audit, freeze, score-dev, diagnostics, select-context, score-test, geometry, evaluate, report. This computational manuscript was generated for project review and has not been submitted to bioRxiv or peer reviewed.

## References

[1] Cuomo ASE et al. Rare and common genetic variants shape single-cell immune gene expression. Nature (2026), published 28 October 2026. https://doi.org/10.1038/s41586-026-10721-2

[2] Source preprint v2: https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2

[3] TenK10K Phase1 summary release. https://doi.org/10.5281/zenodo.18221260

[4] GENCODE human release44. https://www.gencodegenes.org/human/release_44.html

[5] GI API documentation. https://docs.genomicintelligence.ai/rest-api.md

[6] PromoterAI software and license terms. https://github.com/Illumina/PromoterAI

## Entire development family

| Arm | Gene-aware | Cell correct/all | Cell rho | Variant correct/all |
|---|---|---|---|---|
| baseline | False | 1658/3616 | 0.124 | 676/1467 |
| canonical_cell_only | False | 1611/3616 | 0.088 | 661/1467 |
| cell_plus_verified_cl_id | False | 921/3616 | 0.190 | 517/1467 |
| verified_cl_id_only | False | 801/3616 | 0.037 | 400/1467 |
| structured_verified_assay | False | 1608/3616 | 0.101 | 681/1467 |
| baseline_plus_gene_symbol | True | 1627/3616 | 0.093 | 665/1467 |
| baseline_plus_ensembl_id | True | 1612/3616 | 0.105 | 667/1467 |
| baseline_plus_symbol_and_ensembl | True | 1585/3616 | 0.064 | 652/1467 |
| historical_repeated_paraphrase | False | 1593/3616 | 0.086 | 711/1467 |
| duplicated_baseline | False | 1600/3616 | 0.070 | 686/1467 |





## Legacy and independent study comparison

| Cohort | Status | All-pair variant direction | Conditional | Cross-cell all-pair |
|---|---|---|---|---|
| Legacy original100 | Exploratory, amended | 52/100 | 52/77 | 32/81 (17 genes, trial cell contrasts) |
| Fresh test baseline | Frozen independent API experiment | 2392/5863 | 2392/4131 | 6195/14446 |

Legacy uncertainty, recomputed with the same locus-bootstrap code: variant 52/100 = 52.0% (95% locus CI 42.6%–61.6%); conditional 67.5% among 77; zero 23, missing 0; rho 0.370 (CI 0.197–0.531); trial cell contrasts 32/81 = 39.5% (95% locus CI 24.1%–55.7%); conditional 45.1% among 71; zero 10, missing 0; rho -0.040 (CI -0.378–0.252). The legacy cell contrasts use the comparison-trial subset, not all 100 pilot genes; their ascertainment differs from the new test.




