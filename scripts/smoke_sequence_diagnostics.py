"""Live real-allele geometry checks, independent of observed effect estimates."""
import json
from pathlib import Path
from gxi_tenk10k.genome.annotations import load_tss
from gxi_tenk10k.genome.reference import prepare
from gxi_tenk10k.genome.variants import covered
from gxi_tenk10k.gi.client import Client
from gxi_tenk10k.util import save


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--flank",type=int,default=40960);args=parser.parse_args()
    out=Path('results/live-20261029')
    annotation=load_tss('data/raw/gencode44/gencode.v44.basic.annotation.gtf.gz')
    client=Client('data/cache',json.load(open(out/'gi_openapi.json')),'g0-expression-8192',max_requests=760)
    rows=[]
    for gene,pos,ref,alt in [('ENSG00000119535',36478697,'C','A'),('ENSG00000087237',56971389,'C','T')]:
        chrom,tss,strand=annotation[gene]
        row=dict(chrom=chrom,tss1=tss,strand=strand,pos1=pos,ref=ref,alt=alt)
        r,a,tss_index,index=prepare('data/raw/gencode44/GRCh38.primary_assembly.genome.fa',row,args.flank)
        context='Human peripheral blood CD14 classical monocyte, unstimulated, 10x 3-prime single-cell RNA-seq'
        reference=client.predict(r,context,tss_index);alternative=client.predict(a,context,tss_index)
        rows.append(dict(row,ref_prediction=reference,alt_prediction=alternative,
            delta=alternative['score']-reference['score'],input_length=len(r),variant_index=index,tss_index=tss_index,
            covered_in_both=covered(index,reference['window'],len(r)) and covered(index,alternative['window'],len(a)),
            equal_windows=reference['window']==alternative['window']))
    save(out/f'smoke_geometry_{args.flank}.json',rows)
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
