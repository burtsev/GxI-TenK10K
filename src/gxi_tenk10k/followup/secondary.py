"""Separate geometry diagnostics and descriptive covariate analyses."""
import asyncio,json
from pathlib import Path
import numpy as np,pandas as pd,pysam
from ..util import save,digest,file_hash
from ..genome.reference import window
from ..genome.variants import oriented_pair,covered
from .client import AsyncGI
from .runner import FASTA

def freeze_geometry(out):
    out=Path(out);f=pd.read_parquet(out/'cohort_development.parquet');f=f[f.role=='candidate'].copy();f['order']=[digest([20261030,r.locus,r.pair_id,r.cell]) for r in f.itertuples()]
    f=f.sort_values('order').drop_duplicates('locus').head(150);target=out/'geometry_cohort.parquet'
    if target.exists():
        if pd.read_parquet(target).observation_id.tolist()!=f.observation_id.tolist():raise ValueError('Geometry subset changed')
    else:f.to_parquet(target,index=False)
    return f
async def geometry(out):
    out=Path(out);f=freeze_geometry(out);reg=pd.read_parquet(out/'prompt_registry.parquet');contexts={(r.gene,r.cell):r.description for r in reg[(reg.split=='development')&(reg.arm=='baseline')].itertuples()}
    client=AsyncGI('data/cache/followup-geometry',json.load(open(out/'gi_openapi.json')));fa=pysam.FastaFile(FASTA)
    path=out/'geometry_results.jsonl';rows=[json.loads(l) for l in path.open()] if path.exists() else [];done={r['observation_id'] for r in rows}
    try:
        for row in f.to_dict('records'):
            if row['observation_id'] in done:continue
            result=dict(observation_id=row['observation_id'],locus=row['locus'],gene=row['gene'],cell=row['cell'],SNP=row['SNP'],beta=row['beta'],status='complete')
            start,end=window(int(row['tss1']),row['strand'],40960)
            if start<0 or end>fa.get_reference_length(row['chrom']):result.update(status='contig_boundary')
            else:
                s=fa.fetch(row['chrom'],start,end).upper();ref,alt,tss,index=oriented_pair(s,start,int(row['pos1']),row['ref'],row['alt'],int(row['tss1']),row['strand'])
                a=await client.predict(ref,contexts[row['gene'],row['cell']],tss);b=await client.predict(alt,contexts[row['gene'],row['cell']],tss)
                result.update(ref_window=a['window'],alt_window=b['window'],variant_index=index,tss_index=tss,delta=b['score']-a['score'],valid_comparison=a['window']==b['window'] and covered(index,a['window'],len(ref)),ref_cache_key=a['cache_key'],alt_cache_key=b['cache_key'])
                if not result['valid_comparison']:result['status']='unequal_or_uncovered_effective_window'
            with path.open('a') as h:h.write(json.dumps(result)+'\n')
            rows.append(result)
            if len(rows)%10==0:print('geometry',len(rows),'of',len(f),flush=True)
    finally:await client.close();fa.close()
    from .analysis import observations,summary
    short=observations(out,'development','baseline');g=pd.DataFrame(rows).merge(short[['observation_id','delta']].rename(columns={'delta':'delta_40960'}),on='observation_id',validate='one_to_one')
    g.to_parquet(out/'geometry_results.parquet',index=False);valid=g[g.valid_comparison.fillna(False)]
    m={'n_preselected':len(f),'n_attempted':len(g),'status_counts':g.status.value_counts().to_dict(),'n_exact_effective_spans':int(len(valid)),'primary_40960':summary(g.beta,g.delta_40960,g.locus,B=3000),'sensitivity_81920_valid_only':summary(valid.beta,valid.delta,valid.locus,B=3000),'not_pooled':True}
    save(out/'geometry_metrics.json',m);return m

def secondary(out):
    from .analysis import observations,summary,contrasts,rho
    from .data import CELLS
    out=Path(out);lock=json.load(open(out/'context_selection.json'));results={}
    for arm in dict.fromkeys(['baseline',lock['selected_arm']]):
        f=observations(out,'test',arm);candidate=f[f.role=='candidate'].copy();candidate['family']=candidate.cell.map(lambda c:CELLS[c][1]);candidate['abs_distance']=abs(candidate.distance)
        strata={}
        for field in ['maf','gc','abs_distance','PIP','SE']:
            # Fixed domain edges or descriptive quantiles only; no retuning primary.
            x=candidate[field];bins=pd.qcut(x,4,duplicates='drop') if x.nunique()>1 else pd.Series('constant',index=x.index)
            strata[field]={str(k):summary(g.beta,g.delta,g.locus,B=2000) for k,g in candidate.groupby(bins,observed=True)}
        family={k:summary(g.beta,g.delta,g.locus,B=2000) for k,g in candidate.groupby('family')}
        percell={k:summary(g.beta,g.delta,g.locus,B=2000) for k,g in candidate.groupby('cell')}
        heuristics={}
        for field,sign in [('abs_distance',-1),('maf',1),('gc',1)]:heuristics[field]={'magnitude_spearman_observed':rho(abs(candidate.beta.to_numpy()),sign*candidate[field].to_numpy()),'allele_direction':'unavailable without allele features'}
        # Descriptive precision weighting uses observed SE, never input selection.
        ok=np.isfinite(candidate.delta);weights=1/(candidate.SE.to_numpy()**2);from .analysis import correct
        weighted=float(np.sum(weights*correct(candidate.beta,candidate.delta))/weights.sum())
        results[arm]={'strata':strata,'family':family,'per_cell_ci':percell,'heuristics':heuristics,'precision_weighted_all_pair_direction':weighted,'cell_agnostic_contrast':{'correct_fraction':0.,'nonzero_coverage':0.,'explanation':'A shared variant effect across cells implies exactly zero contrast.'}}
    results['external_baselines']={name:{'status':'not_run','reason':reason} for name,reason in [('PromoterAI','No authorized executable model or licensed matched precomputed predictions supplied; published median is an unmatched narrative comparison.'),('AlphaGenome','No authenticated executable external-model access configured for the exact frozen cohort.'),('Borzoi','No validated local checkpoint, output-to-cell mapping and matching executable workflow configured.'),('caQTL','No independent matched caQTL source included in the authorized Zenodo/source release.'),('regulatory_pathways','No prespecified verified pathway-membership source in the frozen allowed-source manifest; do not invent post hoc pathway labels.')]}
    save(out/'secondary_metrics.json',results);return results
