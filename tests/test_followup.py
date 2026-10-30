import asyncio,json
import numpy as np
import pandas as pd
import pytest,httpx
from gxi_tenk10k.followup.data import cluster_genes
from gxi_tenk10k.followup.prompts import descriptions,ARMS,GENE_ARMS
from gxi_tenk10k.followup.client import AsyncGI
from gxi_tenk10k.gi.client import Blocked
from gxi_tenk10k.genome.variants import oriented_pair

def test_transitive_history_clusters():
    a=pd.DataFrame([('a','chr1',100000),('b','chr1',290000),('c','chr1',480000),('d','chr1',800000)],columns=['gene','chrom','tss1'])
    groups,old=cluster_genes(a,a.iloc[:1]);assert groups['a']==groups['c'];assert groups['c'] in old;assert groups['d'] not in old

def test_arm_registry():
    m={'B':dict(canonical='naive B cell',cl_id='CL:0000788',ontology_eligible=True)}
    r=dict(cell='B',gene='ENSG1',symbol='ABC');p=descriptions(r,m)
    assert list(p)==ARMS and len(p)==10
    for arm in set(ARMS)-GENE_ARMS:assert 'ABC' not in p[arm] and 'ENSG1' not in p[arm]
    assert p['duplicated_baseline']==p['baseline']+' '+p['baseline']
    m['B']['ontology_eligible']=False;assert descriptions(r,m)['verified_cl_id_only'] is None

def test_minus_palindromic_orientation():
    ref,alt,tss,index=oriented_pair('ACGTA',100,102,'C','G',103,'-')
    assert ref=='TACGT' and alt=='TACCT' and tss==2 and index==3
    with pytest.raises(ValueError,match='REF mismatch'):oriented_pair('ACGTA',100,102,'A','T',103,'-')

def test_async_cache_retry_and_no_total_cap(tmp_path,monkeypatch):
    monkeypatch.setenv('GI_API_KEY','test-only');seq='A'*9198;context='human';contract={'info':{'version':'test'}}
    async def run():
        client=AsyncGI(tmp_path,contract,retries=1);client.interval=0;calls=[]
        def response(request):
            calls.append(1)
            if len(calls)==1:return httpx.Response(429,headers={'retry-after':'0'})
            return httpx.Response(200,json={'data':{'task':'expression','model':client.model,'input':{'description':context,'tss_index':4599},'prediction':{'unit':'log(TPM+1)','expression_log_tpm':1}},'meta':{'sequence_length':len(seq),'task_specific_counts':{'tss_index':4599,'scored_window':[0,len(seq)]}}})
        await client.http.aclose();client.http=httpx.AsyncClient(transport=httpx.MockTransport(response),base_url=client.base)
        a,b=await asyncio.gather(client.predict(seq,context,4599),client.predict(seq,context,4599))
        assert a==b and len(calls)==2 and client.hits==1 and client.rate_events==1
        assert not hasattr(client,'max_requests');await client.close()
    asyncio.run(run())

def test_persisted_failure_bounds_retry(tmp_path,monkeypatch):
    monkeypatch.setenv('GI_API_KEY','test-only')
    async def run():
        c=AsyncGI(tmp_path,{'info':{'version':'test'}},retries=0);c.interval=0
        await c.http.aclose();c.http=httpx.AsyncClient(transport=httpx.MockTransport(lambda r:httpx.Response(403)),base_url=c.base)
        for _ in range(2):
            with pytest.raises(Blocked):await c.predict('A'*9198,'human',4599)
        assert c.calls==1;await c.close()
    asyncio.run(run())
from gxi_tenk10k.followup.analysis import correct,summary,holm,contrasts,primary_permutation

def test_all_pair_and_missing_denominators():
    m=summary([1,-1,1,1],[.1,0,np.nan,-.1],['a','a','b','b'],B=100)
    assert m['n_all']==4 and m['n_correct']==1 and m['n_zero']==1 and m['n_missing']==1
    assert m['all_pair_accuracy']==.25 and m['conditional_accuracy']==.5
    assert holm([.01,.04])==[.02,.04]

def test_cluster_bootstrap_reproducible():
    a=summary([1,2,3,4],[1,-1,2,0],['a','a','b','b'],B=100)
    assert a==summary([1,2,3,4],[1,-1,2,0],['a','a','b','b'],B=100)

def test_cell_contrast_and_permutation(tmp_path):
    p=pd.DataFrame([dict(split='test',pair_id='p',locus='a',cell_a='A',cell_b='B')]);p.to_csv(tmp_path/'preselected_cell_pairs.tsv',sep='\t',index=False)
    f=pd.DataFrame([dict(pair_id='p',cell='A',beta=1,SE=.01,delta=.2,locus='a',observation_id='x'),dict(pair_id='p',cell='B',beta=-1,SE=.01,delta=-.2,locus='a',observation_id='y')])
    c=contrasts(tmp_path,'test',f);assert c.observed.iloc[0]==2 and c.predicted.iloc[0]==.4 and c.precision_subset.iloc[0]
    a=primary_permutation(f,c,R=30);assert a==primary_permutation(f,c,R=30);assert 0<a['p_A']<=1

def test_contrast_quantization_avoids_float_rank_artifacts(tmp_path):
    p=pd.DataFrame([dict(split='test',pair_id=k,locus=k,cell_a='A',cell_b='B') for k in ['p','q']]);p.to_csv(tmp_path/'preselected_cell_pairs.tsv',sep='\t',index=False)
    f=pd.DataFrame([dict(pair_id=k,cell=cell,beta=1 if cell=='A' else 0,SE=.1,delta=d,locus=k,observation_id=k+cell) for k,ds in [('p',[.3,.2]),('q',[.2,.1])] for cell,d in zip(['A','B'],ds)])
    c=contrasts(tmp_path,'test',f);assert c.predicted.iloc[0]==c.predicted.iloc[1]==.1

def test_runner_checkpoint_resume_and_geometry_failure(tmp_path,monkeypatch):
    from gxi_tenk10k.followup import runner
    from gxi_tenk10k.util import save
    from gxi_tenk10k.gi.client import cache_key
    frame=pd.DataFrame([dict(gene='g',chrom='chr1',tss1=26001,strand=s,pos1=26001,ref='A',alt='C',observation_id=s,role='candidate',locus='l',cell='B',pair_id='p',beta=1,SE=.1) for s in ['+','-']])
    frame.to_parquet(tmp_path/'cohort_development.parquet',index=False)
    save(tmp_path/'api_preflight.json',{'passed':True});save(tmp_path/'gi_openapi.json',{'info':{'version':'fixture'}})
    monkeypatch.setattr(runner,'freeze',lambda out:None)
    monkeypatch.setattr(runner,'registry',lambda out:pd.DataFrame([dict(split='development',gene='g',cell='B',arm=arm,description='human B cell') for arm in ['baseline','canonical_cell_only']]))
    class FA:
        def __init__(self,*args):pass
        def fetch(self,c,a,b):return 'A'*(b-a)
        def close(self):pass
    monkeypatch.setattr(runner.pysam,'FastaFile',FA)
    class Client:
        total=0;concurrency=2;calls=0;hits=0;rate_events=0;bad=False
        def __init__(self,*args):pass
        async def predict(self,seq,context,tss):
            type(self).total+=1
            return {'score':1+seq.count('C')*.1+seq.count('G')*.1,'window':[0,100] if self.bad else [0,40960],'cache_key':str(type(self).total)}
        def adapt(self):pass
        async def close(self):pass
    monkeypatch.setattr(runner,'AsyncGI',Client)
    a=asyncio.run(runner.score(tmp_path,'development',arms=['baseline']));assert len(a)==2 and np.allclose(a.delta,.1)
    calls=Client.total;b=asyncio.run(runner.score(tmp_path,'development',arms=['baseline']));assert len(b)==2 and Client.total==calls
    # An uncovered ALT/REF interval must fail without accepting a new job.
    Client.bad=True
    with pytest.raises(ValueError,match='scored geometry'):asyncio.run(runner.score(tmp_path,'development',arms=['canonical_cell_only']))
