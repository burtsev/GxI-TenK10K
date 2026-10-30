# Start Codex for the GI × TenK10K independent follow-up

Read `AGENTS.md` and `docs/FOLLOWUP_PROTOCOL.md`, then follow `docs/FOLLOWUP_CODEX_RUNBOOK.md` and `configs/followup.yaml`. **The primary objective is cross-cell-type allelic-effect prediction**, followed by independent direction replication and context-description ablations.

Before editing or starting fresh API scoring, **recover the preceding pilot code, VM-local artifacts and all historically scored loci**. The historical report mentions commit `b7ecbed`, but GitHub main inspected at protocol-writing time did not include the benchmark code. Keep old output intact and commit it to a separate recovery branch if accessible. Do not reuse old test loci as new confirmatory validation. If the ledger cannot be recovered, perform data preparation/exploration but flag independent confirmation as blocked.

Run with a user-supplied `GI_API_KEY` via the documented GI REST/official client, not the ChatGPT skill or MCP. **Key has no total API quota**: remove the old 790/800/1600 attempt ceilings; preserve rate limits, retries, caching, adaptive concurrency and operational logging. Preflight auth, model, 40,960-bp identical effective windows and source checksums. Freeze selection, prompt registry, thresholds and test cohort before scoring.

Implement and execute all feasible stages now. Target development 250 loci and held-out >=1000 high-PIP variant–gene pairs across >=600 independent locus clusters and >=12 immune populations where possible; use up to six cell descriptions/pair with complete measured betas. Keep a fully independent test, optimize only on dev, select one fixed gene-free context alternative (or baseline) and evaluate test only once. Add separate gene-ID conditioning diagnostics, controls, full zero-effect denominators, matched-model baselines if legally executable, bootstrap/permutation inference and an audit-ready scientific report. Do not tune the holdout based on outcomes. Save checkpoints so the run can be resumed. Never fabricate results.

On VM:
```bash
cd ~/GxI-TenK10K
git status
git fetch origin
# IMPORTANT: do not reset/overwrite a VM checkout containing uncommitted study outputs.
export GI_API_KEY='<set securely in your shell / secret manager>'
export GI_BASE_URL='https://api.genomicintelligence.ai'
# Start Codex in this working directory, then give it this file as a task.
```

Suggested Codex message:
`Read docs/FOLLOWUP_START.md and fully execute docs/FOLLOWUP_CODEX_RUNBOOK.md following docs/FOLLOWUP_PROTOCOL.md. Recover previous trial artifacts, freeze disjoint development and independent confirmatory cohorts, then perform direct-GI-API scoring with my unquotaed GI_API_KEY, rate-limit-aware adaptive concurrency and checkpointing. Run full validation and write rigorous results. Commit code, configuration, tests and sanitized summary to this repo. Never invent results or reuse any previously tested locus in the holdout.`
