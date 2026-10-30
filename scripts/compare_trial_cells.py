"""Additional prespecified lineage contexts for the first 20 trial candidates."""
import json
from pathlib import Path
import pandas as pd
from gxi_tenk10k.pipeline.select import CELL_MAP
from gxi_tenk10k.pipeline.score import score_pair
from gxi_tenk10k.gi.client import Client,Blocked
from gxi_tenk10k.analysis.contrasts import cell_contrasts
from gxi_tenk10k.analysis.metrics import evaluate
from gxi_tenk10k.util import save


def main():
    out=Path('results/live-20261029');p=pd.read_parquet(out/'predictions.parquet')
    trial=p[p.category=='candidate'].head(20)
    pool=pd.read_parquet(out/'association_pool.parquet')
    client=Client('data/cache',json.loads((out/'gi_openapi.json').read_text()),'g0-expression-8192',max_requests=760)
    rows=[];exclusions=[];biological={v:k for k,v in CELL_MAP.items()}
    stop=False
    for candidate in trial.to_dict('records'):
        for cell in CELL_MAP.values():
            source=pool[(pool.gene==candidate['gene'])&(pool.SNP==candidate['SNP'])&(pool.cell==cell)]
            if len(source)!=1:
                exclusions.append(dict(gene=candidate['gene'],SNP=candidate['SNP'],cell=cell,reason='missing_comparison_beta'));continue
            row=source.iloc[0].to_dict();row.update(beta_alt=float(row['beta']),locus=candidate['locus'],category='trial_comparison')
            context=f'Human peripheral blood {biological[cell]}, unstimulated, 10x 3-prime single-cell RNA-seq'
            try:rows.append(score_pair(client,'data/raw/gencode44/GRCh38.primary_assembly.genome.fa',row,20480,context))
            except (ValueError,KeyError,Blocked) as e:
                exclusions.append(dict(gene=row['gene'],SNP=row['SNP'],cell=cell,reason=str(e)))
                if isinstance(e,Blocked):stop=True;break
            pd.DataFrame(rows).to_parquet(out/'trial_comparison_predictions.parquet',index=False)
        if stop:break
    frame=pd.DataFrame(rows);contrasts=cell_contrasts(frame)
    contrasts.to_csv(out/'trial_cell_contrasts.tsv',sep='\t',index=False)
    pd.DataFrame(exclusions).to_csv(out/'comparison_exclusions.tsv',sep='\t',index=False)
    save(out/'cell_contrast_metrics.json',evaluate(contrasts,replicates=2000,permutations=1000))
    save(out/'comparison_stage.json',{'pairs':20,'observations':len(rows),'cells':frame.cell.nunique(),'contrasts':len(contrasts),'exclusions':len(exclusions),'new_requests':client.calls,'cumulative_requests':client.spent})

if __name__=='__main__':main()
