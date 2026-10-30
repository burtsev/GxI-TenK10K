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
    monkeypatch.setattr(runner,'validate_reference',lambda out:None)
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

def test_cell_permutation_is_shared_across_locus_despite_row_order(tmp_path):
    p=pd.DataFrame([dict(split='test',pair_id=k,locus='L',cell_a='A',cell_b='B') for k in ['p','q']]);p.to_csv(tmp_path/'preselected_cell_pairs.tsv',sep='\t',index=False)
    rows=[]
    for pair,order in [('p',['A','B']),('q',['B','A'])]:
        for cell in order:rows.append(dict(pair_id=pair,cell=cell,beta=1 if cell=='A' else -1,SE=.1,delta=.2 if cell=='A' else -.2,locus='L',observation_id=pair+cell))
    f=pd.DataFrame(rows);c=contrasts(tmp_path,'test',f);m=primary_permutation(f,c,R=50)
    assert {x[0] for x in m['null']}=={0.,1.}

def test_paired_context_comparison_identical_models_has_zero_gain():
    from gxi_tenk10k.followup.analysis import paired_comparison
    f=pd.DataFrame({'observation_id':list('abcd'),'locus':['L1','L1','L2','L2'],'beta':[1,2,-1,-2],'delta':[.1,.2,-.1,-.2]})
    c=pd.DataFrame({'locus':f.locus,'observed':f.beta,'predicted':f.delta})
    m=paired_comparison(f,f,c,c,B=50,R=30)
    assert m['variant']['gain']==0 and m['variant']['rho_gain']==0
    assert m['cross_cell']['ci95']==[0.,0.] and m['cross_cell']['rho_gain_ci95']==[0.,0.]
    assert m['holm_p']==[1.,1.]

def test_selection_binding_is_immutable_on_resume(tmp_path):
    from gxi_tenk10k.followup.binding import validate_selection
    from gxi_tenk10k.followup.analysis import dev_select
    from gxi_tenk10k.util import save,file_hash
    source=tmp_path/'analysis_fixture.py';source.write_text('unchanged')
    save(tmp_path/'code_manifest.json',{str(source):file_hash(source)})
    lock={'selected_arm':'baseline','file_hashes':{'code_manifest.json':file_hash(tmp_path/'code_manifest.json')}};save(tmp_path/'context_selection.json',lock)
    assert dev_select(tmp_path)==lock
    before=(tmp_path/'code_manifest.json').read_bytes();source.write_text('changed')
    with pytest.raises(ValueError,match='Frozen implementation'):dev_select(tmp_path)
    assert (tmp_path/'code_manifest.json').read_bytes()==before
    source.write_text('unchanged');(tmp_path/'code_manifest.json').write_text('{}')
    with pytest.raises(ValueError,match='Frozen selection input'):validate_selection(tmp_path)

def test_reference_binding_detects_sequence_and_coordinate_changes(tmp_path):
    from gxi_tenk10k.followup.binding import validate_reference
    from gxi_tenk10k.util import save,file_hash
    fa=tmp_path/'fixture.fa';fa.write_text('>chr1\nACGT\n');index=tmp_path/'fixture.fa.fai';index.write_text('fixture-index')
    save(tmp_path/'data_verification.json',{'reference':{'sha256':file_hash(fa)},'reference_index_sha256':file_hash(index)})
    validate_reference(tmp_path,fa);index.write_text('different-offset')
    with pytest.raises(ValueError,match='coordinate index'):validate_reference(tmp_path,fa)
    fa.write_text('>chr1\nTGCA\n')
    with pytest.raises(ValueError,match='FASTA changed'):validate_reference(tmp_path,fa)

def test_context_choice_obeys_gene_free_coverage_and_zero_rho_tie_break():
    from gxi_tenk10k.followup.analysis import select_gene_free
    def metric(acc,rho,coverage=100):return {'cross_cell':{'all_pair_accuracy':acc,'spearman':rho,'n_all':100,'n_scored':coverage},'variant':{'all_pair_accuracy':.6}}
    m={'baseline':metric(.4,.1),'canonical_cell_only':metric(.44,0.),'structured_verified_assay':metric(.44,-.05),'baseline_plus_gene_symbol':metric(.99,.9),'duplicated_baseline':metric(.45,.2,96),'historical_repeated_paraphrase':metric(.42,.5)}
    selected,eligible=select_gene_free(m);assert selected=='canonical_cell_only';assert set(eligible)=={'canonical_cell_only','structured_verified_assay'}
    assert select_gene_free({'baseline':m['baseline'],'historical_repeated_paraphrase':m['historical_repeated_paraphrase']})==('baseline',[])
