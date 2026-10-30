"""Unrestricted-total, rate-paced async GI client with durable checkpoints."""
import asyncio,json,os,time
from pathlib import Path
from datetime import datetime,timezone
from email.utils import parsedate_to_datetime
import httpx
from ..gi.client import cache_key,validate_response,Blocked
from ..util import digest,save


class AsyncGI:
    def __init__(self,root,contract,model='g0-expression-8192',timeout=120,retries=5):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.contract=contract;self.model=model;self.contract_hash=digest(contract)
        self.base=os.getenv('GI_BASE_URL','https://api.genomicintelligence.ai').rstrip('/')
        if not self.base.startswith('https://'):raise ValueError('HTTPS required')
        self.key=os.getenv('GI_API_KEY')
        if not self.key:raise Blocked('GI_API_KEY absent')
        self.http=httpx.AsyncClient(base_url=self.base,timeout=timeout)
        self.retries=retries;self.locks={};self.pace_lock=asyncio.Lock();self.io_lock=asyncio.Lock()
        self.next_request=0.;self.interval=.1;self.cooldown=0.;self.calls=0;self.hits=0;self.rate_events=0
        self.concurrency=2;self.maximum_concurrency=32;self.batch_successes=0
        self.ledger=self.root/'requests.jsonl';self.responses=self.root/'responses.jsonl'
    def key_for(self,seq,context,tss):
        return cache_key(seq,context,self.model,tss,self.contract['info']['version'],self.contract_hash,self.base)
    async def log(self,path,row):
        async with self.io_lock:
            with path.open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
    async def pace(self):
        async with self.pace_lock:
            delay=max(self.next_request,self.cooldown)-time.monotonic()
            if delay>0:await asyncio.sleep(delay)
            self.next_request=time.monotonic()+self.interval
    async def predict(self,seq,context,tss,force=False):
        if not 9198<=len(seq)<=500000 or not 4599<=tss<=len(seq)-4599:raise ValueError('Invalid input geometry')
        if set(seq)-set('ACGTN') or not context.strip():raise ValueError('Invalid sequence/context')
        key=self.key_for(seq,context,tss);path=self.root/(key+'.json')
        lock=self.locks.setdefault(key,asyncio.Lock())
        async with lock:
            if path.exists() and not force:
                self.hits+=1;return dict(validate_response(json.load(path.open()),seq,context,self.model,tss),cache_key=key)
            fail=self.root/(key+'.failure.json')
            if fail.exists() and not force:raise Blocked('Persisted request failure: '+json.load(fail.open())['reason'])
            body={'sequence':seq,'model':self.model,'tss_index':tss,'options':{'description':context}}
            for attempt in range(self.retries+1):
                await self.pace();self.calls+=1
                await self.log(self.ledger,{'cache_key':key,'attempt':attempt,'model':self.model,'timestamp':datetime.now(timezone.utc).isoformat(),'sequence_hash':digest(seq),'context_hash':digest(context),'input_bp':len(seq),'tss_index':tss})
                try:r=await self.http.post('/v1/tasks/expression/predict',json=body,headers={'Authorization':'Bearer '+self.key})
                except httpx.TransportError as e:
                    await self.log(self.responses,{'cache_key':key,'status':'transport_error','type':type(e).__name__})
                    if attempt==self.retries:save(fail,{'reason':'Transport retries exhausted','attempts':attempt+1});raise Blocked('Transport retries exhausted') from None
                    await asyncio.sleep(min(2**attempt,30));continue
                rate={k:v for k,v in r.headers.items() if k.startswith('ratelimit') or k in ('retry-after','x-request-id')}
                await self.log(self.responses,{'cache_key':key,'status':r.status_code,'rate_headers':rate})
                # Official contract: burst capacity=rpm/6; sustained rps=burst/10.
                if 'ratelimit-limit' in rate:
                    try:self.interval=max(.125,1/(float(rate['ratelimit-limit'])/10*.75))
                    except (ValueError,ZeroDivisionError):pass
                if r.status_code==429:
                    self.rate_events+=1;self.concurrency=max(2,self.concurrency//2)
                if r.status_code in (429,500,502,503,504) and attempt<self.retries:
                    value=r.headers.get('retry-after',str(min(2**attempt,30)))
                    try:seconds=float(value)
                    except ValueError:seconds=(parsedate_to_datetime(value)-datetime.now(timezone.utc)).total_seconds()
                    if seconds>300:save(fail,{'reason':'Long Retry-After; resume after service recovery'});raise Blocked('Retry-After operational stop')
                    self.cooldown=max(self.cooldown,time.monotonic()+max(seconds,0));await asyncio.sleep(max(seconds,0));continue
                if r.status_code!=200:
                    reason=f'HTTP {r.status_code}';save(fail,{'reason':reason,'attempts':attempt+1})
                    raise Blocked(reason)
                data=r.json()
                try:result=validate_response(data,seq,context,self.model,tss)
                except (ValueError,KeyError) as e:
                    save(fail,{'reason':'Response integrity: '+str(e)});raise
                save(path,data);self.batch_successes+=1
                return dict(result,cache_key=key)
            raise Blocked('Bounded retries exhausted')
    def adapt(self,had_errors=False):
        if had_errors:self.concurrency=max(2,self.concurrency//2)
        elif self.batch_successes>=32:
            self.concurrency=min(self.maximum_concurrency,self.concurrency*2);self.batch_successes=0
    async def close(self):await self.http.aclose()
