"""Frozen-description repeat of the fixed matched study, then all available cells."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from gxi_tenk10k.gi.client import Client,Blocked
from gxi_tenk10k.pipeline.score import score_pair
from gxi_tenk10k.pipeline.select import CELL_MAP
from gxi_tenk10k.analysis.metrics import evaluate
from gxi_tenk10k.analysis.contrasts import cell_contrasts
from gxi_tenk10k.util import save,file_hash

OUT=Path('results/full-best-description-20261030');BASE=Path('results/live-20261029')


def paired(a,b):
    f=a[['gene','SNP','cell','locus','beta_alt','delta']].merge(b[['gene','SNP','cell','delta']],on=['gene','SNP','cell'],suffixes=('_baseline','_best'),validate='one_to_one')
    old=((abs(f.delta_baseline)>1e-6)&(np.sign(f.delta_baseline)==np.sign(f.beta_alt))).astype(float)
    new=((abs(f.delta_best)>1e-6)&(np.sign(f.delta_best)==np.sign(f.beta_alt))).astype(float)
    f['change']=new-old
    g=f.groupby('locus').change.agg(['sum','count']);rng=np.random.default_rng(20261030);boot=[]
    for _ in range(2000):
        x=rng.integers(len(g),size=len(g));boot.append(float(g['sum'].to_numpy()[x].sum()/g['count'].to_numpy()[x].sum()))
    return {'n_pairs':len(f),'n_loci':len(g),'baseline_correct_all_pairs':int(old.sum()),'best_correct_all_pairs':int(new.sum()),'all_pair_gain':float((new-old).mean()),'gain_ci95':np.quantile(boot,[.025,.975]).tolist()}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    source=pd.read_parquet(BASE/'predictions.parquet');original=source[source.category.isin(['candidate','control'])].copy()
    ids=original[['gene','SNP','category','locus','cell']].rename(columns={'cell':'discovery_cell'})
    pool=pd.read_parquet(BASE/'association_pool.parquet')
    available=ids.merge(pool,on=['gene','SNP'],validate='one_to_many')
    available['beta_alt']=available.beta;available['is_discovery']=available.cell==available.discovery_cell
    if available.duplicated(['gene','SNP','cell']).any():raise ValueError('Duplicate source observations')
    bad=~np.isfinite(available.beta)|~np.isfinite(available.SE)|(available.SE<=0)
    if bad.any():raise ValueError('Invalid source estimates')
    labels={v:k for k,v in CELL_MAP.items()}
    baseline='Human peripheral blood {cell}, unstimulated, 10x 3-prime single-cell RNA-seq'
    best=json.load(open('docs/runs/context-repetition-20261030/frozen_selection.json'))['template']
    config={'best_template':best,'baseline_template':baseline,'model':'g0-expression-8192','flank_bp':20480,'input_bp':40960,'max_cumulative_requests':4000,'expected_complete_observations_per_description':len(available),'potential_observations_per_description':800,'missing_observed_betas_per_description':800-len(available),'quota_preflight':'User/API developer confirmed quota lifted for key and explicitly requested full scoring. Serial bounded requests; rate headers retained. Monetary price not published, no cost estimate claimed.','gate':'Discovery-candidate correct/all-pairs and Spearman must both improve over original baseline, before all-cell expansion.','interpretation':'Outcome-tuned prompt repeat, not a new zero-shot study or independent validation.'}
    save(OUT/'design.json',config);available.to_parquet(OUT/'observations.parquet',index=False)
    expected=ids.assign(join=1).merge(pd.DataFrame({'cell':list(CELL_MAP.values()),'join':1}),on='join').drop(columns='join')
    absent=expected.merge(available[['gene','SNP','cell']],on=['gene','SNP','cell'],how='left',indicator=True)
    exclusions=absent[absent._merge=='left_only'].drop(columns='_merge');exclusions['reason']='missing_observed_beta';exclusions.to_csv(OUT/'exclusions.tsv',sep='\t',index=False)
    client=Client('data/cache',json.load(open(BASE/'gi_openapi.json')),config['model'],max_requests=4000,retries=0)
    initial=client.spent;records=[]
    def score(rows,description):
        template=best if description=='best' else baseline
        for row in rows.to_dict('records'):
            result=score_pair(client,'data/raw/gencode44/GRCh38.primary_assembly.genome.fa',row,20480,template.format(cell=labels[row['cell']]))
            result['description']=description;records.append(result)
            pd.DataFrame(records).to_parquet(OUT/'predictions.parquet',index=False)
            if len(records)%50==0:print('Scored records',len(records),'ledger',client.spent,flush=True)
    discovery=available[available.is_discovery]
    score(discovery,'baseline');score(discovery,'best')
    scored=pd.DataFrame(records);old=scored[(scored.description=='baseline')&(scored.category=='candidate')];new=scored[(scored.description=='best')&(scored.category=='candidate')]
    gate={'paired':paired(old,new),'baseline':evaluate(old),'best':evaluate(new)}
    gate['passed']=gate['paired']['best_correct_all_pairs']>gate['paired']['baseline_correct_all_pairs'] and gate['best']['spearman']>gate['baseline']['spearman']
    save(OUT/'discovery_gate.json',gate)
    print('Discovery gate',gate['passed'],flush=True)
    if gate['passed']:
        rest=available[~available.is_discovery]
        score(rest,'baseline');score(rest,'best')
    f=pd.DataFrame(records);metrics={}
    for description,g in f.groupby('description'):
        metrics[description]={category:evaluate(x) for category,x in g.groupby('category')}
        metrics[description]['candidate_by_cell']={cell:evaluate(x) for cell,x in g[g.category=='candidate'].groupby('cell')}
        contrasts=cell_contrasts(g[g.category=='candidate'])
        contrasts.to_csv(OUT/(description+'_cell_contrasts.tsv'),sep='\t',index=False)
        metrics[description]['cell_contrasts']=evaluate(contrasts)
    metrics['paired']={category:paired(f[(f.description=='baseline')&(f.category==category)],f[(f.description=='best')&(f.category==category)]) for category in ['candidate','control']}
    metrics['paired']['cell_contrasts']=paired(pd.read_csv(OUT/'baseline_cell_contrasts.tsv',sep='\t'),pd.read_csv(OUT/'best_cell_contrasts.tsv',sep='\t'))
    save(OUT/'metrics.json',metrics)
    f.to_csv(OUT/'measured_predictions.tsv',sep='\t',index=False,na_rep='NA')
    save(OUT/'execution.json',{'completed_records':len(f),'description_counts':f.description.value_counts().to_dict(),'new_requests':client.calls,'ledger_requests':client.spent,'total_requests_including_initial_preledger':client.spent+1,'resume_start_ledger':initial,'gate_passed':gate['passed'],'completed_all_available_cells':len(f)==2*len(available),'source_hash':file_hash(BASE/'association_pool.parquet')})
    print(json.dumps(json.load(open(OUT/'execution.json')),indent=2),flush=True)

if __name__=='__main__':main()
