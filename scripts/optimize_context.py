"""Exploratory context tuning; freeze selection before fresh-locus validation."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from gxi_tenk10k.util import save,digest
from gxi_tenk10k.gi.client import Client
from gxi_tenk10k.pipeline.score import score_pair
from gxi_tenk10k.pipeline.select import CELL_MAP
from gxi_tenk10k.analysis.metrics import evaluate

OUT=Path('results/context-optimization-20261030')
BASE=Path('results/live-20261029')
TEMPLATES={
'baseline':'Human peripheral blood {cell}, unstimulated, 10x 3-prime single-cell RNA-seq',
'structured':'assay term name is single-cell RNA sequencing. library preparation is 10x Genomics Chromium 3-prime RNA-seq. biosample summary is Homo sapiens peripheral blood mononuclear cells, {cell}, unstimulated.',
'compact':'Homo sapiens; peripheral blood mononuclear cells (PBMCs); {cell}; unstimulated; 10x Genomics Chromium 3-prime single-cell RNA sequencing.',
'biosample_first':'biosample summary is unstimulated human {cell} from peripheral blood mononuclear cells. assay term name is single-cell RNA sequencing using 10x Genomics 3-prime poly(A) capture.'}
NAMES={v:k for k,v in CELL_MAP.items()}

def panel(f,n,prior=()):
    chosen=[];used=set();prior=list(prior)
    for cell,g in f.groupby('cell',sort=True):
        g=g.copy();g['order']=[digest([20261030,r.gene,r.SNP,r.cell]) for r in g.itertuples()]
        count=0
        for r in g.sort_values('order').to_dict('records'):
            if r['gene'] in used:continue
            if any(r['chrom']==x['chrom'] and abs(r['tss1']-x['tss1'])<=200000 for x in prior+chosen):continue
            chosen.append(r);used.add(r['gene']);count+=1
            if count==n:break
        if count!=n:raise ValueError('Insufficient independent loci for '+cell)
    return pd.DataFrame(chosen).drop(columns='order')

def utility(f):
    return int(((abs(f.delta)>1e-6)&(np.sign(f.delta)==np.sign(f.beta_alt))).sum())

def paired(old,new):
    f=old[['gene','SNP','cell','delta','beta_alt']].merge(new[['gene','SNP','cell','delta']],on=['gene','SNP','cell'],suffixes=('_old','_new'),validate='one_to_one')
    a=((abs(f.delta_old)>1e-6)&(np.sign(f.delta_old)==np.sign(f.beta_alt))).astype(float)
    b=((abs(f.delta_new)>1e-6)&(np.sign(f.delta_new)==np.sign(f.beta_alt))).astype(float)
    d=(b-a).to_numpy();rng=np.random.default_rng(20261030)
    boot=[float(d[rng.integers(len(d),size=len(d))].mean()) for _ in range(2000)]
    return {'n':len(f),'baseline_correct_all_pairs':int(a.sum()),'selected_correct_all_pairs':int(b.sum()),'difference':float(d.mean()),'ci95':np.quantile(boot,[.025,.975]).tolist(),'note':'Fresh validation loci are nonoverlapping +/-100kb clusters; original-cohort interval is gene-level descriptive.'}

def main(extra_prior=None, incumbent="baseline", cap=790, planned=344, preflight=None):
    OUT.mkdir(parents=True,exist_ok=True)
    allbase=pd.read_parquet(BASE/'predictions.parquet');original=allbase[allbase.category=='candidate'].copy()
    training=panel(original[original.split=='development'],4)
    historical=pd.read_parquet(BASE/'initial_unmatched_pass/predictions.parquet')
    prior=pd.concat([allbase,historical])[['gene','chrom','tss1']].drop_duplicates().to_dict('records')
    if extra_prior is not None:prior += extra_prior.to_dict('records')
    prescreen=pd.read_parquet(BASE/'expanded_prescreen.parquet')
    prescreen=prescreen[abs(prescreen.pos1-prescreen.tss1)<20380]
    fresh_ids=panel(prescreen,5,prior)
    pool=pd.read_parquet(BASE/'association_pool.parquet')
    fresh=fresh_ids[['gene','SNP','cell']].merge(pool,on=['gene','SNP','cell'],validate='one_to_one')
    fresh['beta_alt']=fresh.beta;fresh['locus']=fresh.gene
    if len(fresh)!=20 or fresh.beta.isna().any():raise ValueError('Invalid fresh cohort')
    training[['gene','SNP','cell','locus']].to_csv(OUT/'development_panel.tsv',sep='\t',index=False)
    fresh[['gene','SNP','cell','tss1']].to_csv(OUT/'fresh_panel.tsv',sep='\t',index=False)
    config={'templates':TEMPLATES,'seed':20261030,'development_n':16,'fresh_n':20,'selection':'Maximize correct directions / all 16 pairs; abstentions=incorrect; retain the declared incumbent on primary ties; Spearman breaks ties only among strictly better alternatives.','incumbent':incumbent,'cap':cap,'quota_preflight':preflight,'initial_preledger_calls':1,'worst_case_new_calls':planned,'interpretation':'Outcome-tuned exploratory experiment; original zero-shot results unchanged.'}
    save(OUT/'design.json',config)
    c=Client('data/cache',json.load(open(BASE/'gi_openapi.json')),'g0-expression-8192',max_requests=cap,retries=0)
    start=c.spent
    if start+planned>cap and not (OUT/'frozen_selection.json').exists():raise ValueError('Budget preflight fails')
    results=[]
    def score(frame,prompt,stage):
        output=[]
        for row in frame.to_dict('records'):
            r=score_pair(c,'data/raw/gencode44/GRCh38.primary_assembly.genome.fa',row,20480,TEMPLATES[prompt].format(cell=NAMES[row['cell']]))
            r.update(prompt=prompt,stage=stage);output.append(r);results.append(r)
            pd.DataFrame(results).to_parquet(OUT/'predictions.parquet',index=False)
        return pd.DataFrame(output)
    metrics={};frames={}
    for prompt in TEMPLATES:
        f=score(training,prompt,'development');frames[prompt]=f
        metrics[prompt]=dict(evaluate(f),correct_all_pairs=utility(f))
        save(OUT/'development_metrics.json',metrics);print(prompt,metrics[prompt]['correct_all_pairs'],metrics[prompt]['spearman'],flush=True)
    alternatives=[p for p in TEMPLATES if metrics[p]['correct_all_pairs']>metrics[incumbent]['correct_all_pairs']]
    winner=max(alternatives,key=lambda p:(metrics[p]['correct_all_pairs'],metrics[p]['spearman'] or -2)) if alternatives else incumbent
    frozen={'winner':winner,'template':TEMPLATES[winner],'design_hash':digest(config),'development_metrics':metrics[winner],'calls_at_freeze':c.spent}
    if (OUT/'frozen_selection.json').exists() and json.load(open(OUT/'frozen_selection.json'))['winner']!=winner:raise ValueError('Frozen winner changed')
    save(OUT/'frozen_selection.json',frozen)
    old=score(fresh,'baseline','fresh');new=score(fresh,winner,'fresh')
    validation={'baseline':evaluate(old),'selected':evaluate(new),'paired':paired(old,new)}
    save(OUT/'fresh_metrics.json',validation)
    full=score(original,winner,'original_cohort_post_tuning')
    save(OUT/'summary.json',{'winner':winner,'baseline_original':evaluate(original),'selected_original':evaluate(full),'original_paired':paired(original,full),'fresh':validation,'new_calls':c.calls,'ledger_calls':c.spent,'total_calls_including_preledger':c.spent+1,'fresh_validated_once_after_freeze':True})
    pd.DataFrame(results).drop(columns=['CHR'],errors='ignore').to_csv(OUT/'measured_predictions.tsv',sep='\t',index=False,na_rep='NA')
    print('Complete',winner,c.calls,flush=True)

if __name__=='__main__':main()
