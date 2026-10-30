"""Recompute all descriptive results and figures from saved observations only."""
import json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gxi_tenk10k.analysis.metrics import evaluate,split_loci
from gxi_tenk10k.analysis.figures import figures
from gxi_tenk10k.util import save,file_hash,digest
from gxi_tenk10k.data.parse import read_table


def summarize_live(out):
    out=Path(out);frame=pd.read_parquet(out/'predictions.parquet')
    frame['split']=split_loci(frame,20261029)
    frame.to_parquet(out/'predictions.parquet',index=False)
    observed=[c for c in frame.columns if c not in ['delta','ref_score','alt_score','ref_sequence_hash','alt_sequence_hash','ref_cache_key','alt_cache_key']]
    frame[observed].to_parquet(out/'observations.parquet',index=False)
    import shutil
    shutil.copyfile('data/cache/requests.jsonl',out/'requests.jsonl')
    # Small, sequence-free reproducibility table; preserve exact numeric scores.
    columns=['gene','locus','split','SNP','chrom','pos1','ref','alt','effect_allele','cell','category','beta_alt','SE','PIP','maf','distance','ref_score','alt_score','delta','model','tss_index','variant_index','scored_start','scored_end','ref_sequence_hash','alt_sequence_hash','ref_cache_key','alt_cache_key','context','source_file']
    frame[columns].to_csv(out/'measured_predictions.tsv',sep='\t',index=False)
    summary={c:evaluate(g) for c,g in frame.groupby('category')}
    candidate=frame[frame.category=='candidate']
    summary['candidate_by_cell']={c:evaluate(g) for c,g in candidate.groupby('cell')}
    summary['candidate_by_split']={c:evaluate(g) for c,g in candidate.groupby('split')}
    summary['candidate_pip_99']=evaluate(candidate[candidate.PIP>=.99])
    summary['candidate_by_distance']={label:evaluate(candidate[mask]) for label,mask in [('within_5kb',candidate.distance<=5000),('beyond_5kb',candidate.distance>5000)]}
    controls=frame[frame.category=='control']
    pairs=controls.merge(candidate[['gene','SNP','delta','distance']],on='gene',suffixes=('_control','_candidate'),validate='one_to_one')
    differences=abs(pairs.delta_candidate)-abs(pairs.delta_control)
    nonzero=abs(differences)>1e-6
    if nonzero.any():
        rng=np.random.default_rng(20261029)
        valid=pairs[nonzero].copy();valid['win']=(differences[nonzero]>0).to_numpy().astype(float)
        groups=[g for _,g in valid.groupby('locus')]
        boots=[pd.concat([groups[i] for i in rng.integers(len(groups),size=len(groups))]).win.mean() for _ in range(2000)]
        summary['matched_magnitude_priority']={'epsilon':1e-6,'pairs':len(pairs),'ties':int((~nonzero).sum()),'non_tied':len(valid),'fraction_candidate_larger':float(valid.win.mean()),'ci95':[float(x) for x in np.quantile(boots,[.025,.975])],'interpretation':'Descriptive high-PIP candidate versus low/unreported-PIP uncertain comparator; not causal AUC'}
    summary['distance_baseline']={'spearman_negative_distance_vs_abs_observed_beta':float(spearmanr(-candidate.distance,abs(candidate.beta_alt)).statistic), 'note':'Distance cannot predict allele sign; chance sign baseline is 0.5.'}
    missing_pip=int(controls.PIP.isna().sum())
    summary['control_pip_reporting']={'unreported':missing_pip,'known':len(controls)-missing_pip,'note':'Unreported PIP is not zero.'}
    save(out/'metrics.json',summary)
    figures(frame,summary,out/'figures')
    import importlib.metadata,sys,subprocess,datetime
    env={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'packages':{p.metadata['Name']:p.version for p in importlib.metadata.distributions()},'config_sha256':file_hash('configs/pilot.yaml'),'code_sha256':{str(p):file_hash(p) for root in ['src','scripts','tests','configs'] for p in Path(root).rglob('*') if p.is_file() and p.suffix in ('.py','.yaml')},'paid_post_attempts':sum(1 for _ in open('data/cache/requests.jsonl'))+1,'quota_note':'One initial positive-control POST outside the persistent ledger; every later attempt is in the ledger.'}
    env['code_sha256'].update({name:file_hash(name) for name in ['pyproject.toml','requirements.lock','AGENTS.md','docs/BENCHMARK_PROTOCOL.md','docs/CODEX_RUNBOOK.md']})
    env['implementation_tree_hash']=digest(env['code_sha256']);save(out/'environment.json',env)
    print(json.dumps({k:v for k,v in summary.items() if k in ['candidate','control','matched_magnitude_priority']},indent=2))

    return summary
