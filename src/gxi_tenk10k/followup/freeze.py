"""Freeze paired cohorts, cells and contrasts using observations/input QA only."""
import json,zipfile
from pathlib import Path
from functools import lru_cache
import numpy as np
import pandas as pd
import pysam
from ..data.parse import read_table,variant
from ..genome.reference import window
from ..util import save,file_hash,digest
from .data import CELLS


def freeze(out):
    out=Path(out)
    if (out/'frozen_protocol.json').exists():
        frozen=json.load(open(out/'frozen_protocol.json'))
        for name,h in frozen['file_hashes'].items():
            if file_hash(out/name)!=h:raise ValueError('Frozen file changed: '+name)
        return frozen
    recovery=json.load(open(out/'recovery/recovery_audit.json'))
    if not recovery['complete']:raise ValueError('Historical registry incomplete; cannot freeze confirmation')
    if len(list((out/'summary_pool').glob('*.json')))!=28:raise ValueError('All source cells must finish audit')
    candidates=pd.read_parquet(out/'candidate_reserve.parquet')
    pool=pd.concat([pd.read_parquet(p) for p in sorted((out/'summary_pool').glob('*.parquet'))],ignore_index=True)
    groups=pool.groupby(['gene','cell'],sort=False).indices
    gene_groups=pool.groupby('gene',sort=False).indices
    high_ids=set(candidates.SNP);pip={};genes=set(candidates.gene)
    with zipfile.ZipFile('data/raw/susie_summary.zip') as a:
        for name in a.namelist():
            if not name.endswith('.tsv'):continue
            with a.open(name) as b:r,_=read_table(b);f=pd.DataFrame(r)
            f=f[f.gene.isin(genes)]
            for r in f.itertuples():pip[(r.gene,r.SNP,r.celltype)]=max(pip.get((r.gene,r.SNP,r.celltype),0),float(r.PIP))
    fa=pysam.FastaFile('data/raw/gencode44/GRCh38.primary_assembly.genome.fa')
    @lru_cache(maxsize=4096)
    def reference(gene,chrom,tss,strand):
        start,end=window(int(tss),strand,20480)
        if end>fa.get_reference_length(chrom):raise ValueError('Contig boundary')
        s=fa.fetch(chrom,start,end).upper()
        if len(s)!=40960 or s.count('N')/len(s)>.01:raise ValueError('Reference missing/ambiguous sequence')
        p=np.zeros(len(s)+1,dtype=np.int32);p[1:]=np.cumsum(np.frombuffer(s.encode(),dtype='S1').astype('U1')=='G')+np.cumsum(np.frombuffer(s.encode(),dtype='S1').astype('U1')=='C')
        return start,s,p
    def gc(pos,start,s,p):
        i=int(pos)-1-start
        if i<100 or i+101>len(s) or 'N' in s[i-100:i+101]:return None
        return float((p[i+101]-p[i-100])/201)
    def valid(r):
        return np.isfinite(r['beta']) and np.isfinite(r['SE']) and r['SE']>0 and np.isfinite(r['maf']) and .01-1e-6<=r['maf']<=.5
    excluded=[];matched=[];used=set()
    for c in candidates.to_dict('records'):
        try:
            g=pool.iloc[groups.get((c['gene'],c['discovery_cell']),np.array([],dtype=int))]
            rows=g[g.SNP==c['SNP']]
            if len(rows)!=1:raise ValueError('Missing or duplicate discovery estimate')
            r=rows.iloc[0].to_dict()
            if not valid(r):raise ValueError('Invalid BETA/SE/MAF')
            start,s,p=reference(r['gene'],r['chrom'],r['tss1'],r['strand'])
            if variant(r['SNP'])!=(r['chrom'].removeprefix('chr'),int(r['pos1']),r['ref'],r['alt']):raise ValueError('Source allele conflict')
            if s[int(r['pos1'])-1-start]!=r['ref']:raise ValueError('Candidate REF mismatch')
            rgc=gc(r['pos1'],start,s,p)
            if rgc is None:raise ValueError('Candidate GC unavailable')
            controls=g[(abs(g.maf-r['maf'])<=.05)&(abs(g.distance-r['distance'])<=5000)&~g.SNP.isin(high_ids|used)].copy()
            controls=controls[controls.SNP.str.fullmatch(r'(?:[0-9]+|X|Y):[1-9][0-9]*:[ACGT]:[ACGT]')]
            eligible=[]
            for x in controls.to_dict('records'):
                xp=pip.get((x['gene'],x['SNP'],x['cell']))
                if xp is not None and xp>=.5:continue
                if not valid(x) or x['ref']==x['alt']:continue
                i=int(x['pos1'])-1-start
                if not 0<=i<len(s) or s[i]!=x['ref']:continue
                xgc=gc(x['pos1'],start,s,p)
                if xgc is None or abs(xgc-rgc)>.1:continue
                x['gc']=xgc;x['PIP']=xp
                x['match_cost']=abs(x['maf']-r['maf'])/.05+abs(x['distance']-r['distance'])/5000+abs(xgc-rgc)/.1
                x['match_order']=digest([20261030,r['gene'],x['SNP'],r['SNP']]);eligible.append(x)
            if not eligible:raise ValueError('No fixed-covariate control')
            control=min(eligible,key=lambda x:(x['match_cost'],x['match_order']))
            # Preselect cells with complete estimates for both candidate/control.
            allgene=pool.iloc[gene_groups[r['gene']]]
            ca=allgene[allgene.SNP==r['SNP']];co=allgene[allgene.SNP==control['SNP']]
            ca={x['cell']:x for x in ca.to_dict('records') if valid(x)};co={x['cell']:x for x in co.to_dict('records') if valid(x)}
            cells={cell for cell in set(ca)&set(co) if abs(ca[cell]['maf']-co[cell]['maf'])<=.05 and (pip.get((r['gene'],control['SNP'],cell)) is None or pip[(r['gene'],control['SNP'],cell)]<.5)}
            if len(cells)<2:raise ValueError('Fewer than two complete measured cells')
            family=CELLS[r['cell']][1]
            related=sorted((x for x in cells if x!=r['cell'] and CELLS[x][1]==family),key=lambda x:digest([20261030,r['gene'],x]))
            distant=sorted((x for x in cells if x!=r['cell'] and CELLS[x][1]!=family),key=lambda x:(CELLS[x][1],digest([20261030,r['gene'],x])))
            # Discovery + up to two related + round-robin distinct distant families.
            selected=[r['cell']]+related[:2];familyqueues={}
            for cell in distant:familyqueues.setdefault(CELLS[cell][1],[]).append(cell)
            while len(selected)<6 and any(familyqueues.values()):
                for fam in sorted(familyqueues):
                    if familyqueues[fam] and len(selected)<6:selected.append(familyqueues[fam].pop(0))
            for cell in related[2:]:
                if len(selected)<6:selected.append(cell)
            r.update(PIP=float(c['PIP']),gc=rgc,locus=c['locus'],split=c['split'],discovery_cell=r['cell'])
            matched.append({'candidate':r,'control':control,'cells':selected,'candidate_rows':ca,'control_rows':co})
            used.add(control['SNP'])
        except (ValueError,KeyError) as e:
            excluded.append(dict(gene=c['gene'],SNP=c['SNP'],cell=c['discovery_cell'],split=c['split'],reason=str(e)))
    pd.DataFrame(excluded).to_csv(out/'selection_exclusions.tsv',sep='\t',index=False)
    dev=[x for x in matched if x['candidate']['split']=='development'];test=[x for x in matched if x['candidate']['split']=='test']
    def ordered(items,limit,unique):
        selected=[];seen=set()
        # First cover distinct clusters, then remaining pairs if permitted.
        for x in items:
            locus=x['candidate']['locus']
            if locus in seen:continue
            selected.append(x);seen.add(locus)
            if len(selected)==limit:return selected
        if not unique:
            for x in items:
                if any(y['candidate']['gene']==x['candidate']['gene'] and y['candidate']['SNP']==x['candidate']['SNP'] for y in selected):continue
                selected.append(x)
                if len(selected)==limit:break
        return selected
    dev=ordered(dev,250,True)
    target=2000 if len(test)>=2000 else 1000
    test=ordered(test,target,False)
    outputs={};contrasts=[]
    for split,items in [('development',dev),('test',test)]:
        rows=[]
        for item in items:
            c=item['candidate'];pair_id=digest([c['gene'],c['SNP']])[:20]
            for role,key in [('candidate','candidate_rows'),('control','control_rows')]:
                for cell in item['cells']:
                    r=dict(item[key][cell]);r.update(pair_id=pair_id,role=role,split=split,locus=c['locus'],discovery_cell=c['discovery_cell'],PIP=pip.get((r['gene'],r['SNP'],cell)),gc=c['gc'] if role=='candidate' else item['control']['gc'])
                    r['observation_id']=digest([split,pair_id,role,cell])
                    rows.append(r)
            from itertools import combinations
            for a,b in combinations(sorted(item['cells']),2):contrasts.append(dict(split=split,pair_id=pair_id,gene=c['gene'],SNP=c['SNP'],locus=c['locus'],cell_a=a,cell_b=b))
        f=pd.DataFrame(rows);f.to_parquet(out/('cohort_'+split+'.parquet'),index=False)
        f[['observation_id','pair_id','gene','SNP','cell','role','locus','discovery_cell','PIP','symbol','tss1','strand','ref','alt','gc','maf','distance']].to_csv(out/('cohort_'+split+'.tsv'),sep='\t',index=False,na_rep='NA')
        outputs[split]={'pairs':len(items),'loci':len({x['candidate']['locus'] for x in items}),'cells':f.cell.nunique(),'observations':len(f),'cell_counts':f.cell.value_counts().to_dict()}
    cf=pd.DataFrame(contrasts);cf.to_csv(out/'preselected_cell_pairs.tsv',sep='\t',index=False)
    development=pd.read_parquet(out/'cohort_development.parquet');testing=pd.read_parquet(out/'cohort_test.parquet')
    assert not set(development.gene)&set(testing.gene)
    assert not set(development.locus)&set(testing.locus)
    legacy=pd.read_csv(out/'recovery/legacy_exclusions.tsv',sep='\t');assert not set(legacy.gene)&set(testing.gene)
    for f in [development,testing]:
        for r in f[['gene','chrom','tss1']].drop_duplicates().itertuples():
            q=legacy[legacy.chrom==r.chrom]
            assert (abs(q.tss1-r.tss1)>200000).all()
    files=['cohort_development.parquet','cohort_test.parquet','cohort_development.tsv','cohort_test.tsv','preselected_cell_pairs.tsv','cell_mapping.tsv','recovery/legacy_exclusions.tsv']
    protocol={'schema_version':1,'status':'cohorts_frozen_before_followup_inference','counts':outputs,'file_hashes':{name:file_hash(out/name) for name in files},'criteria':'PIP>=.9, per-cell BH gene q<=.05, SNV, explicit ALT beta, reference match, >=2 measured cells, fixed uncertain control match, historical union clusters excluded.','test_target_met':outputs['test']['pairs']>=1000 and outputs['test']['loci']>=600 and outputs['test']['cells']>=12,'target_expansion_2000_available':target==2000,'physical_not_ld_independence':True,'selection_time':'before any followup cohort GI inference'}
    save(out/'cohort_hashes.json',protocol['file_hashes']);save(out/'frozen_protocol.json',protocol)
    flow=pd.DataFrame([dict(stage='candidate_reserve',pairs=len(candidates)),dict(stage='matched_complete',pairs=len(matched)),dict(stage='development_frozen',pairs=len(dev)),dict(stage='test_frozen',pairs=len(test))]);flow.to_csv(out/'cohort_flow.tsv',sep='\t',index=False)
    print(json.dumps(outputs,indent=2),flush=True);return protocol
