import json
import httpx
import pytest
from gxi_tenk10k.gi.client import Client, Blocked, cache_key, validate_response
from gxi_tenk10k.pipeline.score import score_pair

CONTRACT = {'info':{'version':'fixture-v1'}}
SEQ = 'ACGT'*2300
CONTEXT = 'Human B naive 10x RNA'


def response(body, score=2):
    return {'data':{'task':'expression','model':body['model'],
            'input':{'description':body['options']['description'], 'tss_index':body['tss_index']},
            'prediction':{'expression_log_tpm':score,'unit':'log(TPM + 1)'}},
            'meta':{'sequence_length':len(body['sequence']),
                    'task_specific_counts':{'tss_index':body['tss_index'], 'scored_window':[100,9100]}}}


def client(tmp_path, monkeypatch, handler, **kwargs):
    monkeypatch.setenv('GI_API_KEY','offline-fixture')
    return Client(tmp_path, CONTRACT,'model',transport=httpx.MockTransport(handler),sleep=lambda t:None,**kwargs)


def test_retry_cache_resume(tmp_path,monkeypatch):
    calls=[]
    def handler(r):
        calls.append(r)
        if len(calls) == 1: return httpx.Response(429,headers={'retry-after':'0'})
        return httpx.Response(200,json=response(json.loads(r.content)))
    c=client(tmp_path,monkeypatch,handler)
    assert c.predict(SEQ,CONTEXT,4600)['score'] == 2
    assert c.calls == 2
    c2=client(tmp_path,monkeypatch,handler)
    assert c2.predict(SEQ,CONTEXT,4600)['score'] == 2
    assert len(calls)==2 and c2.spent==2
    assert 'offline-fixture' not in ''.join(p.read_text() for p in tmp_path.iterdir())


@pytest.mark.parametrize('code',[401,403,422,500,429])
def test_failures_bounded(tmp_path,monkeypatch,code):
    c=client(tmp_path,monkeypatch,lambda r:httpx.Response(code),retries=1,max_requests=2)
    with pytest.raises(Blocked):c.predict(SEQ,CONTEXT,4600)
    assert c.calls <=2
    assert not list(tmp_path.glob('[a-f0-9]'*64+'.json'))


def test_budget_cumulative(tmp_path,monkeypatch):
    c=client(tmp_path,monkeypatch,lambda r:httpx.Response(503),max_requests=1)
    with pytest.raises(Blocked): c.predict(SEQ,CONTEXT,4600)
    c=client(tmp_path,monkeypatch,lambda r:pytest.fail('No budget'),max_requests=1)
    with pytest.raises(Blocked): c.predict(SEQ,CONTEXT,4600)
    assert c.calls==0


def test_missing_key(tmp_path,monkeypatch):
    monkeypatch.delenv('GI_API_KEY',raising=False)
    c=Client(tmp_path,CONTRACT,'model')
    with pytest.raises(Blocked):c.predict(SEQ,CONTEXT,4600)
    assert c.calls==0


def test_hash_inputs():
    key=lambda seq=SEQ,ctx=CONTEXT,model='m',tss=4600:cache_key(seq,ctx,model,tss,'v','h','https://api')
    assert key()==key()
    assert len({key(),key(ctx='other'),key(model='m2'),key(seq=SEQ[:-1]+'A'),key(tss=4601)})==5


@pytest.mark.parametrize('alter', ['unit','window','model','context','tss','length','nan'])
def test_response_validation(alter):
    b={'model':'m','sequence':SEQ,'tss_index':4600,'options':{'description':CONTEXT}}
    d=response(b)
    if alter=='unit':d['data']['prediction']['unit']='TPM'
    if alter=='window':d['meta']['task_specific_counts']['scored_window']=None
    if alter=='model':d['data']['model']='wrong'
    if alter=='context':d['data']['input']['description']='wrong'
    if alter=='tss':d['data']['input']['tss_index']=4599
    if alter=='length':d['meta']['sequence_length']=3
    if alter=='nan':d['data']['prediction']['expression_log_tpm']=float('nan')
    with pytest.raises(ValueError):validate_response(d,SEQ,CONTEXT,'m',4600)


def test_offline_pair_integration(tmp_path,monkeypatch):
    import subprocess
    path=tmp_path/'ref.fa';path.write_text('>1\n'+'A'*20000+'\n')
    subprocess.run(['samtools','faidx',str(path)],check=True)
    calls=[]
    def handler(r):
        b=json.loads(r.content);calls.append(b)
        score=2+b['sequence'].count('C')
        return httpx.Response(200,json=response(b,score))
    c=client(tmp_path/'cache',monkeypatch,handler)
    row={'chrom':'1','tss1':10000,'strand':'+','pos1':10001,'ref':'A','alt':'C'}
    result=score_pair(c,path,row,4600,CONTEXT)
    assert result['delta']==1
    assert calls[0]['options']==calls[1]['options']
    score_pair(c,path,row,4600,CONTEXT)
    assert len(calls)==2
    row['pos1']=5400
    with pytest.raises(ValueError,match='outside effective'):score_pair(c,path,row,4600,CONTEXT)
