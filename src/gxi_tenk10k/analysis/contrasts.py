"""Descriptive cell contrasts; no independent-SE assumption across donor cells."""
from itertools import combinations
import pandas as pd


def cell_contrasts(frame):
    required={'gene','SNP','cell','beta_alt','delta','locus'}
    if not required<=set(frame):raise ValueError('Missing contrast columns')
    if frame.duplicated(['gene','SNP','cell']).any():raise ValueError('Ambiguous cell observations')
    rows=[]
    for (gene,snp),g in frame.groupby(['gene','SNP'],sort=True):
        for a,b in combinations(g.sort_values('cell').to_dict('records'),2):
            rows.append(dict(gene=gene,SNP=snp,locus=a['locus'],cell=a['cell']+' minus '+b['cell'],
                beta_alt=a['beta_alt']-b['beta_alt'],delta=a['delta']-b['delta'],
                note='Descriptive contrast; covariance of observed estimates unknown'))
    return pd.DataFrame(rows,columns=['gene','SNP','locus','cell','beta_alt','delta','note'])
