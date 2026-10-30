"""Final input/provenance acceptance checks and sanitized inference registry."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from ..util import save,file_hash
from .freeze import freeze
from .runner import export

def qa(out):
    out=Path(out);frozen=freeze(out);dev=pd.read_parquet(out/'cohort_development.parquet');test=pd.read_parquet(out/'cohort_test.parquet');legacy=pd.read_csv(out/'recovery/legacy_exclusions.tsv',sep='\t')
    checks={'historical_registry_complete':json.load(open(out/'recovery/recovery_audit.json'))['complete'],'no_dev_test_genes':not bool(set(dev.gene)&set(test.gene)),'no_dev_test_clusters':not bool(set(dev.locus)&set(test.locus)),'no_historical_genes':not bool(set(legacy.gene)&set(pd.concat([dev,test]).gene)),'frozen_targets_met':frozen['test_target_met']}
    for split,f in [('development',dev),('test',test)]:
        checks[split+'_unique_observations']=not f.observation_id.duplicated().any();checks[split+'_finite_effects']=bool(np.isfinite(f.beta).all() and np.isfinite(f.SE).all() and (f.SE>0).all());checks[split+'_at_most_six_cells']=bool(f.groupby(['pair_id','role']).cell.nunique().le(6).all())
        a=f[f.role=='candidate'];b=f[f.role=='control'];matched=a.merge(b,on=['pair_id','cell'],suffixes=('_candidate','_control'),validate='one_to_one')
        checks[split+'_matched_all_rows']=len(matched)==len(a)==len(b)
        checks[split+'_maf_match']=bool((abs(matched.maf_candidate-matched.maf_control)<=.05+1e-7).all())
        checks[split+'_distance_match']=bool((abs(matched.distance_candidate-matched.distance_control)<=5000).all())
        checks[split+'_gc_match']=bool((abs(matched.gc_candidate-matched.gc_control)<=.10+1e-7).all())
        discovery=a[a.cell==a.discovery_cell];checks[split+'_discovery_pip']=bool((discovery.PIP>=.9).all())
    alleles=pd.read_csv(out/'allele_checks.tsv',sep='\t');checks['reference_all_pass']=bool(alleles.ref_pass.all())
    unit=(out/'offline_tests.txt').read_text().strip().splitlines()[-1];checks['offline_tests_passed']='passed' in unit and 'failed' not in unit
    rows=[]
    for split in ['development','test']:
        f=export(out,split)
        for r in f[f.status=='complete'].itertuples():
            for allele in ['ref','alt']:
                key=getattr(r,allele+'_cache_key');rows.append(dict(cache_key=key,split=split,arm=r.arm,observation_id=r.observation_id,gene=r.gene,SNP=r.SNP,cell=r.cell,locus=r.locus,role=r.role,allele=allele))
    # Geometry responses are separate experiments but also historical exposure.
    geometry=pd.read_parquet(out/'geometry_results.parquet')
    for r in geometry.itertuples():
        for allele in ['ref','alt']:
            key=getattr(r,allele+'_cache_key',None)
            if isinstance(key,str):rows.append(dict(cache_key=key,split='development_geometry',arm='baseline',observation_id=r.observation_id,gene=r.gene,SNP=r.SNP,cell=r.cell,locus=r.locus,role='candidate',allele=allele))
    registry=pd.DataFrame(rows).drop_duplicates();registry.to_parquet(out/'historical_registry_for_future_studies.parquet',index=False)
    keys=sorted(registry.cache_key.unique());response_hashes={}
    for key in keys:
        candidates=[Path('data/cache')/name/(key+'.json') for name in ['followup-v1','followup-geometry','followup-smoke']]
        path=next(p for p in candidates if p.exists());response_hashes[key]=file_hash(path);save(out/'response_hashes.json',response_hashes)
    checks['all_inference_cache_files_present']=len(response_hashes)==len(keys)
    import shutil
    predictions=pd.concat([export(out,s) for s in ['development','test']],ignore_index=True);predictions.to_parquet(out/'predictions.parquet',index=False)
    predictions[predictions.status!='complete'].reindex(columns=['observation_id','split','arm','gene','SNP','cell','status','reason']).to_csv(out/'coverage_exclusions.tsv',sep='\t',index=False)
    pd.concat([dev,test],ignore_index=True).to_parquet(out/'observed_effects.parquet',index=False)
    shutil.copy2(out/'recovery/legacy_exclusions.tsv',out/'legacy_exclusions.tsv');shutil.copy2('scripts/reproduce_followup.sh',out/'reproduce.sh')
    result={'checks':checks,'passed':sum(bool(v) for v in checks.values()),'failed':sum(not bool(v) for v in checks.values()),'reference_variants_checked':len(alleles),'reference_failures':int((~alleles.ref_pass).sum()),'unit_tests':unit,'response_hashes_sha256':file_hash(out/'response_hashes.json'),'distinct_inference_responses':len(keys)}
    save(out/'qa.json',result)
    if result['failed']:raise ValueError('Scientific acceptance QA failed: '+str([k for k,v in checks.items() if not v]))
    return result
