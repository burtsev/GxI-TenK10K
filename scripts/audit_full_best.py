"""Offline scientific QA on the full paired description repeat."""
from pathlib import Path
import numpy as np
import pandas as pd
from gxi_tenk10k.data.parse import variant
from gxi_tenk10k.util import save


def main():
    root=Path('results/full-best-description-20261030');f=pd.read_parquet(root/'predictions.parquet');o=pd.read_parquet(root/'observations.parquet')
    for r in o.itertuples():
        assert variant(r.SNP)==(r.chrom.removeprefix('chr'),int(r.pos1),r.ref,r.alt)
        assert r.effect_allele==r.alt and r.genome_build=='GRCh38'
    assert not f.duplicated(['gene','SNP','cell','description']).any()
    assert len(f)==2*len(o)
    assert (f.tss_index==20480).all() and (f.scored_start==0).all() and (f.scored_end==40960).all()
    assert ((f.variant_index>=0)&(f.variant_index<40960)).all()
    assert np.isfinite(f[['beta_alt','SE','ref_score','alt_score','delta']]).all().all()
    a=f[f.description=='baseline'];b=f[f.description=='best'];paired=a.merge(b,on=['gene','SNP','cell'],suffixes=('_baseline','_best'),validate='one_to_one')
    for c in ['beta_alt','SE','ref_sequence_hash','alt_sequence_hash','locus']:
        assert paired[c+'_baseline'].equals(paired[c+'_best']),c
    assert np.allclose(f.delta,f.alt_score-f.ref_score,rtol=0,atol=1e-12)
    save(root/'quality_audit.json',{'passed':True,'observation_rows':len(o),'prediction_rows':len(f),'paired_rows':len(paired),'allele_orientation':'All SNP IDs, genomic REF/ALT and effect alleles agree; minus-strand editing retained ALT-minus-REF.','geometry':'Every oriented TSS index 20480, effective interval [0,40960), variant covered.','sequence_control':'Exact REF/ALT sequence hashes equal across descriptions.','finite_scores':True,'duplicate_rows':0,'prediction_units':'log(TPM+1), checked by GI adapter on each response.'})
    print('Full scientific QA passed',len(paired),'paired rows')

if __name__=='__main__':main()
