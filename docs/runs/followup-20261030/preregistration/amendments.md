# Resolved protocol and operational notes

- Physical transitive TSS clusters are locus proxies; matched LD data are unavailable.
- Frozen cohort: 250 development loci; 1,000 test pairs / 777 test clusters. Eligible matched test reserve did not reach 2,000; PIP threshold unchanged.
- Historical-only smoke executed after input-only cohort freeze, before any new cohort request. No holdout outcome used.
- Concurrency scheduling increased conservatively from 2 to at most 32 in-flight requests with 8 requests/s ceiling, below documented 10 requests/s edge limit; halves on 429. Total calls unrestricted. Interrupted development once to apply this scheduling change; cache/checkpoints retained.
- Exact co-primary p-value definitions and calibration frozen in analysis_plan.json before test; conservative 50% direction floor required in addition to clustered randomization.
- Ontology arms ineligible for unmatched exact labels. Gene-aware arms and 25-locus nonbiological diagnostics exploratory.
- Before test inference, added six fixed immune regulatory pathway strata from official Reactome98, documented as a supported annotation database in the official scDeepID source. Exact symbol membership and checksum in configs/followup_pathways.json. No test predictions inspected; descriptive overlapping strata only.
- Assay source check: Nature Capture and sequencing methods confirms 10x 3-prime GEX HT v3.1 / Chromium X. Unstimulated denotes absence of a described experimental stimulation arm in the sample-processing protocol; it does not assert biologically resting cells or healthy donors.
