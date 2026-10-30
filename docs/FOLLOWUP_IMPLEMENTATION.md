# Independent follow-up implementation

The resolved protocol, input-only cohort membership and recovery audit are in
`docs/runs/followup-20261030/preregistration/`. The authoritative scientific
requirements remain `FOLLOWUP_PROTOCOL.md`. Legacy pilot code and outputs remain
unchanged and were archived before merging the new instructions.

Use the pinned Python environment, export `GI_API_KEY`, and run:

```bash
gxi-tenk10k followup audit
gxi-tenk10k followup freeze
gxi-tenk10k followup score-dev
gxi-tenk10k followup diagnostics
gxi-tenk10k followup select-context
gxi-tenk10k followup score-test
gxi-tenk10k followup geometry
gxi-tenk10k followup evaluate
gxi-tenk10k followup report
```

The default output is `results/followup-20261030`. `scripts/reproduce_followup.sh`
adds official-source retrieval and the full development-family/secondary analyses.
The raw summaries, FASTA and historical artifacts are cached on the study VM.
On a clean clone, the reproduction script seeds the audited historical exclusion
input from the committed preregistration instead of pretending to recover missing
raw caches. `audit` retrieves and verifies the exact pinned ontology file.

`freeze` deterministically builds a surplus reserve, streams all 28 full-summary
members, validates source records and GRCh38 REF, matches controls without model
outcomes, and freezes whole physical TSS clusters. Existing frozen cohorts are
hash-checked rather than overwritten. The attained test target is 1,000 pairs
from 777 clusters; the eligible matched reserve did not support expansion to2,000.
The test and development cohorts share no genes or transitive ±100-kb clusters;
all 167 historical genes and their connected prospective clusters are excluded.
Physical clustering is a dependence proxy, not LD-based proof of independence.

The new async client has no total-call ceiling. It starts at2 concurrent requests,
increases gradually to at most32, paces at no more than8 requests/second and honors
Retry-After with bounded retries. It halves concurrency on429. This leaves room
under the documented10 requests/second edge rate. A request cache key binds the
model, base URL, API version and complete contract digest, exact sequence hash,
TSS geometry, context string and options. Exact duplicate requests are locked and
shared. Successful responses and failure sentinels persist separately; ambiguous
network failures can incur duplicate billing and are recorded. Persistent failures
require operational review rather than indefinite replay. No DNA or API secret is
written into the request ledger.

SQLite WAL checkpoints preserve every finished observation; raw responses are
under `data/cache/followup-v1`. There is a separate smoke and long-window cache.
No successful historical cache entry is overwritten. `score-test` verifies frozen
cohort, prompt, context selection, API contract and source-code hashes. It scores
baseline and at most one frozen gene-free alternative. Incomplete development
families cannot select a context, and incomplete test arms cannot produce a full
confirmatory evaluation. A deliberately partial evaluation requires an explicit
scientific amendment; it is not the default success path.

Primary all-pair direction counts numerical zeros and unscored comparisons as
incorrect; conditional accuracy exposes its smaller denominator. Spearman retains
model zeros. All locus resamples include every row from the sampled cluster. Cell
randomization draws common random cell priorities for each locus, inducing a
permutation within each variant's frozen measured-cell set. Different measured
sets can therefore have different induced permutations, while sharing the locus
random priorities. Allele signs flip jointly by locus. Co-primaryA uses the larger
one-sided direction/rho permutation p, co-primaryB uses a one-sided sign-flip p;
Holm adjusts the two. A direction claim additionally requires all-pair accuracy
above50%. The lower bootstrap bound for cross-cell rho is explicitly reported.
Paired baseline/alternative gains use a separate two-endpoint Holm family. The
nine development comparisons receive a separate nine-arm Holm adjustment and do
not constitute independent confirmation.

Gene labels and deliberately wrong labels are exploratory model-behavior
experiments, not evidence of a DNA-only benchmark. Ontology arms require exact
verified CL label/ID matches; missing mappings remain in the full coverage
denominator. Six fixed Reactome98 pathway strata use exact GENCODE44 symbols and
are descriptive overlapping groups; their membership was fixed before test.
Long-window sensitivity uses150 prespecified development loci. Unequal REF/ALT
effective spans are recorded as invalid apparent deltas and never pooled with
primary effects. Calibration is fitted on development only. Controls have low or
unknown PIP and are not established causal negatives. External model and caQTL
analyses explicitly report absence of executable matching data or licensed access.
