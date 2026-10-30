"""Account for every historical request, including unpaired reference diagnostics."""
import json
from pathlib import Path
import pandas as pd
from ..util import save,file_hash,digest
from ..gi.client import cache_key
from ..genome.reference import fasta_sequence,window,prepare
from ..genome.variants import reverse_complement


def recover(out: Path) -> dict:
    out.mkdir(parents=True,exist_ok=True)
    cache=Path('data/cache');ledger=[json.loads(l) for l in (cache/'requests.jsonl').read_text().splitlines()]
    keys={r['key'] for r in ledger};mapped={};rows=[];sources=[]
    for path in sorted(Path('results').rglob('*.parquet')):
        if out in path.parents or path.name in ('association_pool.parquet','observations.parquet','expanded_prescreen.parquet'):continue
        f=pd.read_parquet(path)
        if not {'gene','chrom','tss1'}<=set(f):continue
        sources.append({'path':str(path),'sha256':file_hash(path),'rows':len(f)})
        for r in f.to_dict('records'):
            row={k:r.get(k) for k in ['gene','chrom','tss1','strand','SNP','cell']};row['source']=str(path);rows.append(row)
            for name in ['ref_cache_key','alt_cache_key']:
                if isinstance(r.get(name),str):mapped[r[name]]=row
    contract=json.load(open('results/live-20261029/gi_openapi.json'))
    fasta='data/raw/gencode44/GRCh38.primary_assembly.genome.fa'
    diagnostics=Path('results/live-20261029/smoke_geometry_40960.json')
    for path in [diagnostics,Path('results/live-20261029/smoke_geometry_diagnostics.json'),Path('results/live-20261029/smoke_geometry_20480.json')]:
        if not path.exists():continue
        for r in json.load(open(path)):
            gene='ENSG00000119535' if r['chrom']=='chr1' else 'ENSG00000087237'
            row=dict(gene=gene,chrom=r['chrom'],tss1=r['tss1'],strand=r['strand'],SNP=f"{r['chrom'].removeprefix('chr')}:{r['pos1']}:{r['ref']}:{r['alt']}",cell='CD14_Mono',source=str(path));rows.append(row)
            ref,alt,tss,index=prepare(fasta,r,r['input_length']//2)
            context='Human peripheral blood CD14 classical monocyte, unstimulated, 10x 3-prime single-cell RNA-seq'
            for seq in [ref,alt]:
                k=cache_key(seq,context,'g0-expression-8192',tss,contract['info']['version'],digest(contract),'https://api.genomicintelligence.ai')
                if k in keys:mapped[k]=row
    # Any unaccounted keys must be recognized HBB geometry diagnostics, never
    # guessed to be irrelevant. Cache response gives context and submitted size.
    hbb={'gene':'ENSG00000244734','chrom':'chr11','tss1':5227071,'strand':'-','SNP':None,'cell':'K562','source':'official HBB positive-control diagnostics'}
    for key in keys-set(mapped):
        data=json.load(open(cache/(key+'.json')));meta=data['meta'];entry=data['data'];ctx=entry['input']['description'];length=meta['sequence_length']
        if 'K562' not in ctx:continue
        a,b=window(hbb['tss1'],'-',length//2);seq=reverse_complement(fasta_sequence(fasta,'chr11',a,b))
        expected=cache_key(seq,ctx,entry['model'],entry['input']['tss_index'],contract['info']['version'],digest(contract),'https://api.genomicintelligence.ai')
        if expected==key:mapped[key]=hbb;rows.append(hbb)
    unmatched=sorted(keys-set(mapped))
    f=pd.DataFrame(rows).dropna(subset=['gene','chrom','tss1']);f.tss1=f.tss1.astype(int)
    f=f.drop_duplicates(['gene','chrom','tss1','SNP','cell','source'])
    f.to_csv(out/'legacy_exclusions.tsv',sep='\t',index=False,na_rep='NA')
    requestmap=[dict(cache_key=k,**mapped[k]) for k in sorted(keys&set(mapped))]
    pd.DataFrame(requestmap).to_parquet(out/'historical_request_registry.parquet',index=False)
    audit={'historical_ledger_requests':len(ledger),'unique_keys':len(keys),'mapped_keys':len(keys&set(mapped)),'unmapped_keys':unmatched,'historical_genes':f.gene.nunique(),'historical_tss':len(f[['chrom','tss1']].drop_duplicates()),'preledger_diagnostics':{'count':1,'gene':'ENSG00000244734','record':'results/pilot-20261029/positive_control_response.json'},'complete':not unmatched,'sources':sources}
    # The preledger call is represented even if no cached HBB call was present.
    if hbb['gene'] not in set(f.gene):pd.concat([f,pd.DataFrame([hbb])]).to_csv(out/'legacy_exclusions.tsv',sep='\t',index=False,na_rep='NA')
    save(out/'recovery_audit.json',audit)
    return audit
