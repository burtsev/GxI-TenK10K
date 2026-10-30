"""Live contract discovery plus historical-locus smoke, never holdout tuning."""
import asyncio,json,os,platform,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import httpx,pandas as pd
from ..util import save,file_hash
from ..genome.reference import prepare
from ..genome.variants import covered
from .client import AsyncGI
async def audit(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    import shutil
    if not (out/'source_manifest.json').exists():shutil.copy2('docs/runs/followup-20261030/preregistration/source_manifest.json',out/'source_manifest.json')
    metadata=httpx.get('https://zenodo.org/api/records/18221260',timeout=60);metadata.raise_for_status();save(out/'zenodo_live.json',metadata.json())
    source=json.load(open(out/'source_manifest.json'));ontology=Path('data/raw/cell-ontology.obo')
    if not ontology.exists():
        r=httpx.get(source['ontology'].get('pinned_url',source['ontology']['url']),timeout=120,follow_redirects=True);r.raise_for_status();ontology.write_bytes(r.content)
    if file_hash(ontology)!=source['ontology']['sha256']:raise ValueError('Ontology version mismatch')
    sources={}
    for name,url in [('gi_openapi.json','https://api.genomicintelligence.ai/v1/openapi.json'),('gi_models.json','https://api.genomicintelligence.ai/v1/tasks/expression/models'),('gi_docs.md','https://docs.genomicintelligence.ai/rest-api.md'),('gi_limits.md','https://docs.genomicintelligence.ai/reference/limits.md')]:
        r=httpx.get(url,follow_redirects=True,timeout=60);r.raise_for_status()
        p=out/name
        if p.exists() and p.read_bytes()!=r.content:raise ValueError('Pinned API contract changed; explicit amendment required')
        p.write_bytes(r.content);sources[name]={'url':url,'sha256':file_hash(p)}
    contract=json.load(open(out/'gi_openapi.json'));catalog=json.load(open(out/'gi_models.json'))
    if not any(m['id']=='g0-expression-8192' for m in catalog['models']):raise ValueError('Model not live')
    rows=pd.read_csv('configs/followup_historical_smoke.tsv',sep='\t')
    client=AsyncGI('data/cache/followup-smoke',contract);checks=[]
    try:
        for row in rows.to_dict('records'):
            ref,alt,tss,index=prepare('data/raw/gencode44/GRCh38.primary_assembly.genome.fa',row,20480)
            context=f"Human peripheral blood {row['cell']}, unstimulated, 10x 3-prime single-cell RNA-seq"
            a=await client.predict(ref,context,tss);b=await client.predict(alt,context,tss)
            repeat=await client.predict(ref,context,tss,force=True)
            passed=a['window']==b['window'] and covered(index,a['window'],40960) and tss==20480
            if not passed:raise ValueError('Smoke scored geometry failure')
            checks.append(dict(gene=row['gene'],SNP=row['SNP'],ref_score=a['score'],alt_score=b['score'],delta=b['score']-a['score'],scored_window=a['window'],tss_index=tss,variant_index=index,exact_repeat_equal=a['score']==repeat['score'],passed=passed))
    finally:await client.close()
    import importlib.metadata
    save(out/'environment.json',dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,platform=platform.platform(),git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()},api_key_present=bool(os.getenv('GI_API_KEY'))))
    result={'passed':all(x['passed'] for x in checks),'checks':checks,'sources':sources,'total_call_cap':None,'unrestricted_total_usage':'user API developer authorization','requests':client.calls,'rate_events':client.rate_events,'rate_interval_seconds':client.interval}
    save(out/'api_preflight.json',result);return result
