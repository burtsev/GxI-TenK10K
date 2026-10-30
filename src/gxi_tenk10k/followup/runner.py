"""Resumable fixed-cohort scoring. SQLite stores one durable row per job."""
import asyncio,json,sqlite3,time,subprocess
from pathlib import Path
from functools import lru_cache
import numpy as np
import pandas as pd
import pysam
from ..util import save,file_hash,digest
from ..genome.reference import window
from ..genome.variants import oriented_pair,covered
from ..gi.client import Blocked
from .client import AsyncGI
from .freeze import freeze
from .prompts import registry,ARMS
FASTA='data/raw/gencode44/GRCh38.primary_assembly.genome.fa'
def database(out):
    db=sqlite3.connect(Path(out)/'jobs.sqlite');db.execute('PRAGMA journal_mode=WAL')
    db.execute('CREATE TABLE IF NOT EXISTS jobs (id TEXT PRIMARY KEY, split TEXT, arm TEXT, observation_id TEXT, result TEXT)');db.commit();return db

def export(out,split):
    db=database(out);rows=[json.loads(x[0]) for x in db.execute('SELECT result FROM jobs WHERE split=?',(split,))];db.close()
    f=pd.DataFrame(rows);f.to_parquet(Path(out)/f'predictions_{split}.parquet',index=False);return f

def freeze_analysis(out):
    out=Path(out)
    spec={'epsilon':1e-6,'secondary_epsilons':[1e-4,1e-3],'bootstrap_replicates':3000,'permutation_replicates':2000,'seed':20261030,'primary_population':'candidate observations; fixed preselected cell pairs','baseline_primary':True,'co_primary_A':'intersection-union of positive all-pair cell-contrast direction above 0.5 and positive Spearman; p_A=max(one-sided direction permutation p, one-sided rho permutation p)','co_primary_B':'all-pair variant direction above 0.5; one-sided locus sign-flip permutation p','correction':'Holm across A and B; separate Holm family for paired alternate improvements','zero_or_missing':'absolute delta<=epsilon wrong in all-pair direction; numerical zeros retained in rho; missing remain denominator and excluded rho with explicit coverage','selection':'gene-free development arms only; >=3pp all-pair cell contrast direction gain, <=3pp scored-pair coverage loss; tie rho then variant direction; else baseline','calibration':'development-only affine least squares; test evaluation without refitting','independence':'transitive overlapping 200kb TSS windows, LD blocks unavailable','precision_secondary':'abs(beta_a-beta_b)>2*(SE_a+SE_b), descriptive; covariance unavailable'}
    path=out/'analysis_plan.json'
    if path.exists() and json.load(path.open())!=spec:raise ValueError('Analysis plan changed')
    if not path.exists():save(path,spec)
    return spec

async def score(out,split,arms=None,diagnostics=False):
    out=Path(out);freeze(out);freeze_analysis(out);reg=registry(out)
    preflight=json.load(open(out/'api_preflight.json'))
    if not preflight['passed']:raise Blocked('API smoke gate failed')
    if split=='test':
        lock=json.load(open(out/'context_selection.json'))
        for name,h in lock['file_hashes'].items():
            if file_hash(out/name)!=h:raise ValueError('Selection binding changed: '+name)
        for name,h in json.load(open(out/'code_manifest.json')).items():
            if file_hash(name)!=h:raise ValueError('Frozen code changed: '+name)
        arms=list(dict.fromkeys(['baseline',lock['selected_arm']]))
    arms=arms or ARMS
    f=pd.read_parquet(out/f'cohort_{split}.parquet')
    if diagnostics:
        f=f[f.locus.isin(sorted(f.locus.unique(),key=lambda x:digest([20261030,x]))[:25])].copy();arms=['swapped_gene_label','nonsense_gene_label','shuffled_cell_label']
    reg=reg[reg.split==split];contexts={(r.gene,r.cell,r.arm):r.description for r in reg.itertuples()}
    fa=pysam.FastaFile(FASTA)
    @lru_cache(maxsize=4096)
    def seq(chrom,tss,strand):
        start,end=window(tss,strand,20480);return start,fa.fetch(chrom,start,end).upper()
    client=AsyncGI('data/cache/followup-v1',json.load(open(out/'gi_openapi.json')))
    db=database(out);completed={r[0] for r in db.execute('SELECT id FROM jobs')};tasks=[]
    genes=sorted(f.gene.unique());cells=sorted(f.cell.unique())
    for arm in arms:
        for row in f.to_dict('records'):
            context=contexts.get((row['gene'],row['cell'],arm))
            if diagnostics:
                context=contexts[row['gene'],row['cell'],'baseline']
                if arm=='swapped_gene_label':context+=' . Gene: '+genes[(genes.index(row['gene'])+1)%len(genes)]
                elif arm=='nonsense_gene_label':context+=' . Gene: GXI_NONBIOLOGICAL_LABEL_20261030'
                else:context=contexts.get((row['gene'],cells[(cells.index(row['cell'])+1)%len(cells)],'baseline')) or 'Human peripheral blood '+cells[(cells.index(row['cell'])+1)%len(cells)]
            jobid=digest([split,arm,row['observation_id'],context,40960])
            if jobid not in completed:tasks.append((jobid,arm,row,context))
    started=time.time();done=0
    async def one(task):
        jobid,arm,row,context=task
        result=dict(row,arm=arm,status='complete',delta=None,context_hash=digest(context) if pd.notna(context) else None)
        if not isinstance(context,str) or not context:
            result.update(status='ineligible_context',reason='No verified ontology mapping')
        else:
            start,s=seq(row['chrom'],int(row['tss1']),row['strand'])
            ref,alt,tss,index=oriented_pair(s,start,int(row['pos1']),row['ref'],row['alt'],int(row['tss1']),row['strand'])
            if len(ref)!=40960 or tss!=20480:raise ValueError('Primary geometry mismatch')
            a,b=await asyncio.gather(client.predict(ref,context,tss),client.predict(alt,context,tss))
            if a['window']!=b['window'] or not covered(index,a['window'],len(ref)):
                raise ValueError('REF/ALT scored geometry mismatch or uncovered variant')
            result.update(delta=b['score']-a['score'],ref_score=a['score'],alt_score=b['score'],ref_cache_key=a['cache_key'],alt_cache_key=b['cache_key'],scored_window=a['window'],variant_index=index)
        clean=json.loads(pd.Series(result).to_json())
        db.execute('INSERT OR REPLACE INTO jobs VALUES (?,?,?,?,?)',(jobid,split,arm,row['observation_id'],json.dumps(clean)));db.commit()
    try:
        while tasks:
            batch=tasks[:max(1,client.concurrency//2)];tasks=tasks[len(batch):]
            results=await asyncio.gather(*(one(t) for t in batch),return_exceptions=True)
            errors=[e for e in results if isinstance(e,BaseException)]
            if errors:raise errors[0]
            done+=len(batch);client.adapt()
            if done%100<len(batch) or not tasks:
                progress={'split':split,'diagnostics':diagnostics,'completed_now':done,'remaining':len(tasks),'requests_this_process':client.calls,'cache_hits':client.hits,'rate_events':client.rate_events,'concurrency':client.concurrency,'elapsed_seconds':round(time.time()-started,1),'timestamp':time.time()}
                save(out/f'progress_{split}.json',progress);print(json.dumps(progress),flush=True)
    finally:
        db.close();fa.close();await client.close();export(out,split)
    return export(out,split)
