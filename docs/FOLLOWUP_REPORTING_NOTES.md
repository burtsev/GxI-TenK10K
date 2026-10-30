# Follow-up reporting and reproducibility notes

The confirmatory code, cohort and context decision were published before the first
test request at commit `3727141f1f020e40d02ab807326ed610487dbf68`. The original
implementation hashes remain in the committed context freeze. None of those
source files, cohorts, contexts, predictions or primary analysis parameters were
changed after the test began.

All intended inference completed. The final source and inference audit passed
27 checks, including 2,500 REF validations; the offline suite passed 68 tests.
The full VM evidence is under `results/followup-20261030/`.

Two explicitly post-analysis reporting helpers finish the manuscript:

```bash
.venv/bin/python scripts/followup_report_io.py
.venv/bin/python scripts/followup_reporting_supplement.py
```

The first executes the frozen report/QA code with an I/O-only workaround: the
original checksum loop rewrote its complete growing JSON dictionary after every
response. Intermediate writes are suppressed; every response is still hashed,
every original acceptance check runs, and the identical final dictionary is
written once. All 59,676 inference/geometry response hashes are retained. The
primary `metrics.json` hash is checked before and after. The workaround and its
source hash are recorded in `reporting_io_amendment.json`.

The second adds locus-bootstrap uncertainty for the prespecified uncertain-control
ranking, and reports PIP ordering only for the 178 comparisons with both posterior
probabilities available. It reuses the frozen bootstrap implementation, seed and
3,000 replicates; it does not modify primary results or introduce a new primary
hypothesis. This supplementary reporting code was written after primary analysis
and is identified as such in `reporting_supplement.json`. Observed expression
means are absent from the allowed association summaries; observed-expression
stratification was not run and donor-level expression was not reconstructed.

Editorial changes distinguish raw zero-preserving randomizations from the
published 50% all-pair direction gate. The raw tests detect weak orientation
information, while the gated co-primary criteria fail. The manuscript does not
reinterpret those criteria as having passed. The final renderer improves heading
placement and the branching cohort diagram; all numbers come from saved outputs.

To rerun inference and analyses, use the original follow-up CLI stages documented
in `FOLLOWUP_IMPLEMENTATION.md`; replace its last `report` stage with the two
reporting commands above. Raw archives and response payloads are not committed.
The public source manifest and frozen cohort membership allow source recovery;
the VM cache enables exact reanalysis without repeating predictions. Fresh API
queries cannot guarantee identical weights because the service exposes a model
ID but no immutable weight digest. External PromoterAI/AlphaGenome/Borzoi runs and
independent caQTL corroboration remain explicitly unrun; the Nature comparison
is narrative, not a matched head-to-head experiment.
