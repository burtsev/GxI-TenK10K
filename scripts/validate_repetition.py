"""Direct once-versus-twice checks on round-two fresh loci; no re-selection."""
import json
from pathlib import Path
import pandas as pd
from optimize_context import paired,NAMES
from gxi_tenk10k.gi.client import Client
from gxi_tenk10k.pipeline.score import score_pair
from gxi_tenk10k.analysis.metrics import evaluate
from gxi_tenk10k.util import save


def main():
    out=Path('results/context-repetition-20261030');base=Path('results/live-20261029')
    templates=json.load(open(out/'design.json'))['templates']
    ids=pd.read_csv(out/'fresh_panel.tsv',sep='\t')
    pool=pd.read_parquet(base/'association_pool.parquet');fresh=ids.merge(pool,on=['gene','SNP','cell','tss1'],validate='one_to_one')
    fresh['beta_alt']=fresh.beta;fresh['locus']=fresh.gene
    client=Client('data/cache',json.load(open(base/'gi_openapi.json')),'g0-expression-8192',max_requests=1600,retries=0)
    records=[];metrics={}
    winner=json.load(open(out/'frozen_selection.json'))['winner']
    selected_family=winner.removesuffix('_twice') if winner.endswith('_twice') and winner.removesuffix('_twice')+'_once' in templates else 'paraphrase_assay_first'
    families=['combined',selected_family]
    save(out/'repeat_validation_design.json',{'families':families,'reason':'Combined-description test and development-selected repeated paraphrase family; no winner re-selection using validation.','worst_case_new_requests':160,'cap':1600})
    for family in families:
        frames=[]
        for suffix in ['once','twice']:
            prompt=family+'_'+suffix;rows=[]
            for row in fresh.to_dict('records'):
                r=score_pair(client,'data/raw/gencode44/GRCh38.primary_assembly.genome.fa',row,20480,templates[prompt].format(cell=NAMES[row['cell']]))
                r['prompt']=prompt;rows.append(r);records.append(r)
                pd.DataFrame(records).to_parquet(out/'repeat_validation_predictions.parquet',index=False)
            frames.append(pd.DataFrame(rows))
        metrics[family]={'once':evaluate(frames[0]),'twice':evaluate(frames[1]),'paired':paired(*frames)}
        save(out/'repeat_validation_metrics.json',metrics);print(family,metrics[family]['paired'],flush=True)
    pd.DataFrame(records).to_csv(out/'repeat_validation_predictions.tsv',sep='\t',index=False,na_rep='NA')
    save(out/'repeat_validation_calls.json',{'new_calls':client.calls,'ledger_calls':client.spent,'total_calls_including_preledger':client.spent+1})

if __name__=='__main__':main()
