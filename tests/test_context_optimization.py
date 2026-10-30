import importlib.util
from pathlib import Path
import pandas as pd

spec=importlib.util.spec_from_file_location('optimization',Path(__file__).parents[1]/'scripts/optimize_context.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def test_panel_ignores_outcomes_and_excludes_nearby_loci():
    rows=[dict(gene=f'g{i}',SNP=f'1:{1000000+i*500000}:A:C',cell='B',chrom='chr1',tss1=1000000+i*500000,beta_alt=1,delta=1) for i in range(10)]
    f=pd.DataFrame(rows)
    prior=[dict(chrom='chr1',tss1=1000000)]
    a=module.panel(f,3,prior)
    f['delta']=-100;f['beta_alt']=-1
    b=module.panel(f,3,prior)
    assert a.gene.tolist()==b.gene.tolist() and 'g0' not in set(a.gene)

def test_abstentions_cannot_inflate_search_utility():
    f=pd.DataFrame(dict(delta=[0.,1.,-1.],beta_alt=[1.,1.,1.]))
    assert module.utility(f)==1
