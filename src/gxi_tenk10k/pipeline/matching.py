"""Outcome-independent local matching; PIP-unreported controls remain uncertain."""
import numpy as np
from ..genome.reference import window,fasta_sequence
from ..util import digest


def control_pool(row, pool, used, fasta, flank=20480, seed=20261029):
    g=pool[(pool.gene==row['gene'])&(pool.cell==row['cell'])&~pool.SNP.isin(used)].copy()
    g=g[g.SNP.str.fullmatch(r'[0-9XY]+:[1-9][0-9]*:[ACGT]:[ACGT]')]
    g=g[~g.SNP.str.contains(r':(?:A:T|T:A|C:G|G:C)$')]
    g=g[(abs(g.maf-row['maf'])<=.05)&(abs(g.distance-row['distance'])<=5000)&(g.distance<flank-100)]
    g=g[g.PIP.isna()|(g.PIP<.5)]
    start,end=window(int(row['tss1']),row['strand'],flank)
    seq=fasta_sequence(fasta,row['chrom'],start,end)
    def gc(pos):
        i=int(pos)-1-start
        if i<100 or i+101>len(seq):return np.nan
        s=seq[i-100:i+101]
        if 'N' in s:return np.nan
        return (s.count('G')+s.count('C'))/len(s)
    rgc=gc(row['pos1']);g['gc']=[gc(p) for p in g.pos1]
    g=g[abs(g.gc-rgc)<=.1]
    g['order']=[digest([seed,row['gene'],snp]) for snp in g.SNP]
    g['cost']=abs(g.maf-row['maf'])/.05+abs(g.distance-row['distance'])/5000+abs(g.gc-rgc)/.1
    return g.sort_values(['cost','order']),rgc
