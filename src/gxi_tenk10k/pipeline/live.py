"""Execute verified-summary, one-discovery-context pilot under a hard call cap.

The user authorized the benchmark; this caps attempts below 800 because no
account balance/pricing endpoint is published. All quota errors stop expansion.
Run with GI_API_KEY exported; never accepts credentials on the command line.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from gxi_tenk10k.pipeline.prepare import prepare, EFFECT_SOURCE
from gxi_tenk10k.pipeline.select import CELL_MAP
from gxi_tenk10k.pipeline.score import score_pair
from gxi_tenk10k.genome.reference import fasta_sequence, window
from gxi_tenk10k.genome.variants import apply_snv
from gxi_tenk10k.gi.client import Client, Blocked
from gxi_tenk10k.util import save, digest, file_hash
from gxi_tenk10k.analysis.metrics import evaluate
from gxi_tenk10k.pipeline.matching import control_pool
from gxi_tenk10k.data.parse import variant


def run_live(args):
    if not 1<=args.max_requests<=760:raise ValueError('Live workflow is capped at 760 ledger attempts plus one initial preflight POST')
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    config=yaml.safe_load(open(getattr(args,'config','configs/pilot.yaml')))
    root=Path('results/pilot-20261029')
    for name in ['gi_models.json','gi_openapi.json','reference_manifest.json']:
        (out/name).write_bytes((root/name).read_bytes())
    fasta=Path('data/raw/gencode44/GRCh38.primary_assembly.genome.fa')
    if (out/'association_pool.parquet').exists():
        candidates=pd.read_parquet(out/'expanded_prescreen.parquet');pool=pd.read_parquet(out/'association_pool.parquet')
    else:candidates,pool=prepare(config,out)
    if args.prepare_only:return
    catalog=json.loads((out/'gi_models.json').read_text());model=config['models']['primary']
    spec=next(m['bio_spec'] for m in catalog['models'] if m['id']==model)
    if spec['trained_window_bp'] is not None:raise ValueError('Actual model coverage cannot be established')
    amendment=yaml.safe_load(open('configs/live_geometry_amendment.yaml'))
    flank=amendment['geometry']['flank_bp']
    save(out/'execution_amendments.json',amendment)
    save(out/'geometry_amendment.json',{'original_input_bp':81920,'amended_input_bp':40960,'amended_flank_bp':flank,'reason':'Real 81920-bp REF/ALT smoke responses differed by 4-7 bp in reported scored-window endpoints; fixed shorter input preserves identical effective scoring windows.','timing':'Before the 20-variant trial; no efficacy threshold tuning.','interpretation':'Protocol deviation: report amended-geometry results separately from original primary geometry.'})
    client=Client('data/cache',json.loads((out/'gi_openapi.json').read_text()),model,max_requests=args.max_requests)
    biological={v:k for k,v in CELL_MAP.items()}
    stages={};predictions=[];exclusions=[];chosen=[]
    def score(row,category):
        row=dict(row);row.update(category=category,beta_alt=float(row['beta']),locus=row['gene'])
        context=config['models']['context_template'].format(cell_type=biological[row['cell']])
        try:
            chrom,pos,ref,alt=variant(row['SNP'])
            if (chrom,pos,ref,alt)!=(row['chrom'].removeprefix('chr'),int(row['pos1']),row['ref'],row['alt']):
                raise ValueError('Observed marker ID and effect-allele columns disagree')
            if not np.isfinite(row['beta']) or not np.isfinite(row['SE']) or row['SE']<=0:
                raise ValueError('Missing/nonfinite observed beta or invalid SE')
            scored=score_pair(client,fasta,row,flank,context)
            predictions.append(scored)
            pd.DataFrame(predictions).to_parquet(out/'predictions.parquet',index=False)
            return scored
        except (ValueError,KeyError,Blocked) as e:
            exclusions.append(dict(gene=row['gene'],SNP=row['SNP'],cell=row['cell'],category=category,reason=str(e)))
            pd.DataFrame(exclusions).to_csv(out/'exclusions.tsv',sep='\t',index=False)
            if isinstance(e,Blocked):raise
            return None
    # Validate both nominated real alleles against the same reference and source.
    for snp in ['1:36478697:C:A','16:56971389:C:T']:
        rows=pool[(pool.SNP==snp)&(pool.cell=='CD14_Mono')]
        if len(rows)!=1:raise ValueError('Missing/ambiguous nominated smoke observation '+snp)
        r=rows.iloc[0].to_dict()
        if score(r,'smoke') is None:
            stages['smoke']={'status':'failed','reason':'coverage/sequence/response mismatch'}
            save(out/'stages.json',stages);return
    stages['smoke']={'status':'complete','pairs':2,'calls_cumulative':client.spent}
    save(out/'stages.json',stages)
    # Candidate order is prespecified by seed, significance, PIP and nominal
    # coverage. Actual scored coverage is a mandatory exclusion, not tuning.
    for candidate in candidates.to_dict('records'):
        if len(chosen)>=100:break
        if abs(int(candidate['pos1'])-int(candidate['tss1']))>=flank-100:continue
        rows=pool[(pool.gene==candidate['gene'])&(pool.SNP==candidate['SNP'])&(pool.cell==candidate['cell'])]
        if len(rows)!=1:
            exclusions.append(dict(gene=candidate['gene'],SNP=candidate['SNP'],cell=candidate['cell'],category='candidate',reason='missing_or_ambiguous_observed_beta'));continue
        r=rows.iloc[0].to_dict()
        available,_=control_pool(r,pool,{r['SNP']},fasta,flank,config['study']['seed'])
        if available.empty:
            exclusions.append(dict(gene=r['gene'],SNP=r['SNP'],cell=r['cell'],category='candidate',reason='no_covariate_matched_control;outcome_independent_eligibility'))
            continue
        try:result=score(r,'candidate')
        except Blocked as e:
            stages['budget_stop']={'reason':str(e)};break
        if result is None:continue
        chosen.append(r)
        if len(chosen)==20:
            stages['trial']={'status':'complete' if len({x['cell'] for x in chosen})>=4 else 'incomplete','candidate_pairs':20,'cells':len({x['cell'] for x in chosen}),'calls_cumulative':client.spent}
            save(out/'stages.json',stages)
            # Save independent trial metrics before proceeding; no threshold fit.
            frame=pd.DataFrame([p for p in predictions if p['category']=='candidate'])
            save(out/'trial_metrics.json',evaluate(frame,replicates=2000,permutations=1000))
    # Do not scale controls until a 20-pair trial across four cells has passed.
    controls=[]
    if stages.get('trial',{}).get('status')=='complete':
        used={x['SNP'] for x in chosen}
        for r in chosen:
            g,rgc=control_pool(r,pool,used,fasta,flank,config['study']['seed'])
            matched=False
            for control in g.sort_values(['cost','order']).head(10).to_dict('records'):
                control.pop('order');control.pop('cost');control['control_category']='matched_low_or_unreported_PIP_uncertain'
                control['matched_candidate']=r['SNP'];control['matched_gc_difference']=control['gc']-rgc
                try:result=score(control,'control')
                except Blocked as e:stages['budget_stop']={'reason':str(e)};break
                if result is not None:
                    controls.append(control);used.add(control['SNP']);matched=True;break
            if not matched:exclusions.append(dict(gene=r['gene'],SNP=r['SNP'],cell=r['cell'],category='control',reason='no_covered_covariate_match_within_bounded_attempts'))
            if client.spent>=args.max_requests:break
    # Explicitly merge overlapping 200-kb gene windows into conservative loci.
    allrows=pd.DataFrame(predictions)
    if not allrows.empty:
        loci={};last_chr=None;last_end=-1;locus_id=0
        genes=allrows[['gene','chrom','tss1']].drop_duplicates().sort_values(['chrom','tss1'])
        for gene,chrom,tss in genes.itertuples(index=False,name=None):
            if chrom!=last_chr:
                locus_id+=1;last_end=tss+100000
            elif tss-100000>last_end:
                locus_id+=1;last_end=tss+100000
            else:last_end=max(last_end,tss+100000)
            loci[gene]=f'locus_{locus_id}';last_chr=chrom
        allrows['locus']=allrows.gene.map(loci)
        allrows.to_parquet(out/'predictions.parquet',index=False)
        trial=allrows[allrows.category=='candidate'].head(20)
        if len(trial)==20:save(out/'trial_metrics.json',evaluate(trial,replicates=2000,permutations=1000))
    pd.DataFrame(chosen).to_csv(out/'candidates.tsv',sep='\t',index=False)
    pd.DataFrame(controls).to_csv(out/'controls.tsv',sep='\t',index=False)
    pd.DataFrame(exclusions).to_csv(out/'exclusions.tsv',sep='\t',index=False)
    stages['pilot']={'status':'complete' if len(chosen)==len(controls)==100 else 'incomplete','candidate_pairs':len(chosen),'control_pairs':len(controls),'cells':len({x['cell'] for x in chosen}),'calls_cumulative':client.spent,'new_calls_this_invocation':client.calls}
    save(out/'stages.json',stages)
    a=config['analysis'];metrics={}
    for category in ['candidate','control','smoke']:
        f=allrows[allrows.category==category]
        metrics[category]=evaluate(f,a['near_zero_delta_epsilon'],a['bootstrap_replicates'],a['permutation_replicates'],config['study']['seed'])
    metrics['candidate_by_cell']={cell:evaluate(g,replicates=2000,permutations=1000) for cell,g in allrows[allrows.category=='candidate'].groupby('cell')}
    save(out/'metrics.json',metrics)
    save(out/'live_run_manifest.json',{'config':config,'reference':json.loads((out/'reference_manifest.json').read_text()),
        'model':model,'quota':'No account balance/pricing endpoint published. Key active; burst headers observed; total attempts capped below 800 including preflight. User authorized up to 1600; cap reduced to 760.',
        'preflight_paid_calls_outside_cache':1,'matching_amendment':'After initial 100-candidate/76-control pass, match availability was moved before inference; same deterministic source ordering, same covariate tolerances, no beta/delta-based inclusion. Initial unmatched pass retained separately.', 'context_count_per_pair':1,'ascertainment':'Discovery-cell candidates; PIP>=0.9 and significant gene; primary requires strict equal scored windows and REF match',
        'control':'MAF <=.05 difference, absolute distance <=5000bp difference, GC <=.1 difference, same gene/cell; low or unreported PIP; uncertain comparators'})
    print(json.dumps(stages,indent=2),flush=True)
