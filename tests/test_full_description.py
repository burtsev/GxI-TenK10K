import importlib.util
from pathlib import Path
import pandas as pd

spec=importlib.util.spec_from_file_location('full_description',Path(__file__).parents[1]/'scripts/full_best_description.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def test_paired_gain_clusters_repeated_cells():
    a=pd.DataFrame(dict(gene=['g','g','h'],SNP=['s','s','t'],cell=['B','T','B'],locus=['l1','l1','l2'],beta_alt=[1.,1.,1.],delta=[-1.,-1.,1.]))
    b=a.copy();b['delta']=[1.,1.,1.]
    result=module.paired(a,b)
    assert result['n_pairs']==3 and result['n_loci']==2
    assert result['baseline_correct_all_pairs']==1 and result['best_correct_all_pairs']==3
    assert abs(result['all_pair_gain']-2/3)<1e-12
    assert result['gain_ci95']==[0.,1.]
