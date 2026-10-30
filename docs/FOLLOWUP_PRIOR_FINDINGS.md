# Prior evidence and reproducibility caveats (read before follow-up)

Historical findings are taken from the three PDF reports attached to the 29–30 October 2026 thread in [GI Slack #general](https://genomicintelligence.slack.com/archives/C0BCCNHV6UU/p1793313628242929). This text is **source-derived historical evidence**, not a fresh independent reanalysis.

## Original exploratory pilot
- 100 high-PIP candidates, 100 matched uncertain controls, four cell populations; 63,986,931 summary rows scanned.
- Candidate directional agreement 52/77 = 67.5% with 23 model-zero/abstentions; all-candidate rho=0.370.
- Controls 31/60 = 51.7%, with 40/100 abstentions; candidate versus control larger magnitude in 60/76 non-tied comparisons.
- Cross-cell effect contrasts 32/71 = 45.1%, rho = -0.040; 81 comparisons from 63/80 complete cell observations, 17 absent observations recorded.
- 40,960-bp sequence geometry and matched-control eligibility introduced following initial investigation; hence not untouched confirmation. 81,920-bp REF/ALT scored-window endpoints varied and a CETP prediction changed sign across geometries.
- Historical report says 49 tests passed and 444 GI inference calls. It names Git commit `b7ecbed`; this commit was not accessible on GitHub `main` when checked on 30 October 2026. Recover artifacts from the VM before any true holdout selection.

## Context optimization 1
- 16 development pairs, four description candidates; selected compact description.
- Original 100-pair result: conditional sign agreement 67.5% -> 67.6%, Spearman 0.370 -> 0.418; all-pair correct 52/100 -> 50/100 (abstentions 23 -> 26).
- Fresh 20-locus comparison 11/20 baseline vs 12/20 selected; paired 95% CI for +5 percentage point change -10 to +25 pp; not compelling.

## Context optimization 2
- Nine alternative descriptions searched on the same 16 development pairs, selected a repeated sentence paraphrase: 9/16 correct; high overfit risk.
- Original 100 previously seen variants: baseline conditional agreement 67.5% -> 75.0%; rho 0.370 -> 0.442; all-pair correct 52/100 -> 60/100. Post-selection in-sample result, not new validation.
- Independent 20-locus validation: baseline 10/20 versus selected 10/20; paired gain 0 pp (CI -15 to +15 pp).
- Proposed but not executed there: gene symbol, Ensembl ID, Cell Ontology IDs, canonical labels and assay ablations; gene identity may enable memorized gene lookup.

## Critical follow-up interpretation
Main unresolved question is cross-cell allelic differences, not just average gene expression. Do not overclaim direction agreement conditional on nonzero effects. Distinguish model abstentions and missing data, matched uncertainty and genuine causal experiments, pooled Spearman vs per-gene PromoterAI reported median 0.66 (not same benchmark), non-overlapping physical windows vs true LD independence, and zero-shot inference vs training-set nonoverlap.

**Sources**: GI_TenK10K_short_report.pdf, GI_context_optimization_report.pdf, GI_repetition_paraphrasing_report.pdf in the linked Slack thread; doi:10.5281/zenodo.18221260. Track actual previous prediction files and their SHA hashes separately.
