import io
import zipfile
import numpy as np
import pandas as pd
import pytest
from gxi_tenk10k.analysis.metrics import evaluate, split_loci
from gxi_tenk10k.pipeline.select import screen, match_controls
from gxi_tenk10k.data.zenodo import members


def test_missing_and_cluster_denominators():
    df=pd.DataFrame({'beta_alt':[1,1,-1,np.nan,1], 'delta':[2,2,-1,1,0],
                     'locus':['a','a','b','c','d'],'gene':['g','g','h','i','j'],'cell':['B']*5})
    m=evaluate(df,replicates=30,permutations=30)
    assert m['n_total']==5 and m['n_missing']==1 and m['n_abstain']==1 and m['n_scored']==3
    assert m['n_loci']==2 and m['ci95']==[1,1]
    assert evaluate(df,replicates=30,permutations=30)==m
    s=split_loci(df,42)
    assert s[0]==s[1]


def test_zero_data():
    df=pd.DataFrame({k:pd.Series(dtype='float64' if k in ['delta','beta_alt'] else 'object') for k in ['delta','beta_alt','gene','locus','cell']})
    m=evaluate(df)
    assert m['n_scored']==0 and m['sign_concordance'] is None and m['ci95'] is None


def test_deterministic_screen():
    f=pd.DataFrame([dict(SNP=f'1:{i+1}:A:C',gene=f'g{i}',celltype=['B','T'][i%2],PIP=.95) for i in range(20)])
    a=screen(f,1,10)
    b=screen(f.sample(frac=1,random_state=3),1,10)
    assert a.SNP.tolist()==b.SNP.tolist()
    assert len(a)==a.gene.nunique()==10


def test_control_matching():
    c=pd.DataFrame([dict(SNP='x',locus='l',maf=.2,distance=100,gc=.5,PIP=.99)])
    pool=pd.DataFrame([dict(SNP='y',locus='l',maf=.21,distance=200,gc=.51,PIP=.1),dict(SNP='z',locus='other',maf=.2,distance=100,gc=.5,PIP=.01)])
    chosen,excluded=match_controls(c,pool)
    assert chosen.SNP.tolist()==['y'] and excluded.empty


@pytest.mark.parametrize('name',['../escape','/absolute'])
def test_zip_traversal(name):
    b=io.BytesIO()
    with zipfile.ZipFile(b,'w') as z:z.writestr(name,'abc')
    b.seek(0)
    with zipfile.ZipFile(b) as z:
        with pytest.raises(ValueError):list(members(z))


def test_cell_contrasts_keep_allele_orientation():
    from gxi_tenk10k.analysis.contrasts import cell_contrasts
    f=pd.DataFrame([dict(gene='g',SNP='1:1:A:C',cell='A',beta_alt=-1,delta=-.1,locus='l'),dict(gene='g',SNP='1:1:A:C',cell='B',beta_alt=-2,delta=-.3,locus='l')])
    c=cell_contrasts(f)
    assert c.iloc[0].beta_alt==1 and abs(c.iloc[0].delta-.2)<1e-12
    with pytest.raises(ValueError):cell_contrasts(pd.concat([f,f]))


def test_matching_uses_covariates_not_outcomes(tmp_path):
    import subprocess
    from gxi_tenk10k.pipeline.matching import control_pool
    fasta=tmp_path/'reference.fa';fasta.write_text('>1\n'+'A'*100000+'\n')
    subprocess.run(['samtools','faidx',str(fasta)],check=True)
    row=dict(gene='g',cell='B',SNP='1:50000:A:C',chrom='1',tss1=50000,strand='+',pos1=50000,maf=.2,distance=0)
    pool=pd.DataFrame([dict(gene='g',cell='B',SNP='1:50100:A:G',pos1=50100,maf=.21,distance=100,PIP=np.nan,beta=1000,delta=-99)])
    a,_=control_pool(row,pool,{row['SNP']},fasta)
    pool['beta']=-1;pool['delta']=99
    b,_=control_pool(row,pool,{row['SNP']},fasta)
    assert a.SNP.tolist()==b.SNP.tolist()==['1:50100:A:G']
    assert a.PIP.isna().all()


def test_magnitude_metrics_retain_zero_predictions():
    f=pd.DataFrame(dict(beta_alt=[1.,2.,3.],delta=[0.,2.,1.],gene=['a','b','c'],locus=['a','b','c'],cell=['B']*3))
    m=evaluate(f,replicates=20,permutations=20)
    assert m['n_scored']==2 and m['spearman_n']==3
    assert abs(m['spearman']-.5)<1e-12
    assert abs(m['spearman_nonabstaining']+1)<1e-12
