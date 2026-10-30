# Resolved protocol and operational notes

- Physical transitive TSS clusters are locus proxies; matched LD data are unavailable.
- Frozen cohort: 250 development loci; 1,000 test pairs / 777 test clusters. Eligible matched test reserve did not reach 2,000; PIP threshold unchanged.
- Historical-only smoke executed after input-only cohort freeze, before any new cohort request. No holdout outcome used.
- Concurrency scheduling increased conservatively from 2 to at most 32 in-flight requests with 8 requests/s ceiling, below documented 10 requests/s edge limit; halves on 429. Total calls unrestricted. Interrupted development once to apply this scheduling change; cache/checkpoints retained.
- Exact co-primary p-value definitions and calibration frozen in analysis_plan.json before test; conservative 50% direction floor required in addition to clustered randomization.
- Ontology arms ineligible for unmatched exact labels. Gene-aware arms and 25-locus nonbiological diagnostics exploratory.
- Before test inference, added six fixed immune regulatory pathway strata from official Reactome98, documented as a supported annotation database in the official scDeepID source. Exact symbol membership and checksum in configs/followup_pathways.json. No test predictions inspected; descriptive overlapping strata only.
- Assay source check: Nature Capture and sequencing methods confirms 10x 3-prime GEX HT v3.1 / Chromium X. Unstimulated denotes absence of a described experimental stimulation arm in the sample-processing protocol; it does not assert biologically resting cells or healthy donors.
- 2026-10-30T16:03:46.324713+00:00, before test inference: offline QA strengthened shared-locus cell permutations to canonical cell order, canonicalized contrast arithmetic to10 decimal places (below primary epsilon), added descriptive paired-rho gain intervals, and verified marker/allele columns across the retained source pool. No development efficacy metrics or test predictions were inspected to choose these changes.
- Before any diagnostic requests, aligned wrong-gene suffix punctuation and shuffled-cell canonical labels with the correct-label arms, froze exact diagnostic strings, and required completion of all diagnostics before context selection.

- 2026-10-30T18:10:37.178953+00:00: after inference and frozen primary analysis, a reporting-only I/O wrapper suppressed intermediate checksum-manifest rewrites while retaining all hashes/checks and identical final contents. Frozen source files and primary metrics were unchanged. Supplemental reporting added the prespecified matched-ranking uncertainty/PIP-ordering summaries and clarified raw-null versus gated decisions; its post-analysis timing is explicit. See docs/FOLLOWUP_REPORTING_NOTES.md and reporting_io_amendment.json.
