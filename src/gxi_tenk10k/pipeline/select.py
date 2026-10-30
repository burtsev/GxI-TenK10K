"""Prespecified PIP screening; these rows are provisional until all gates pass."""
import pandas as pd
from ..data.parse import variant
from ..util import digest

CELL_MAP = {'B naive': 'B_naive', 'CD4 naive T': 'CD4_Naive',
            'CD14 classical monocyte': 'CD14_Mono', 'CD16 nonclassical monocyte': 'CD16_Mono'}


def screen(susie, seed, target, minimum=.9):
    f = susie.copy()
    f['PIP'] = pd.to_numeric(f.PIP, errors='raise')
    if f.PIP.isna().any() or not f.PIP.between(0, 1).all():
        raise ValueError('Invalid PIP')
    f = f[f.PIP >= minimum].copy()
    parsed = []
    for _, r in f.iterrows():
        try:
            chrom, pos, ref, alt = variant(r.SNP)
        except ValueError:
            continue
        if len(ref) != 1 or len(alt) != 1 or {ref, alt} in ({'A', 'T'}, {'C', 'G'}):
            continue
        parsed.append(dict(r, chrom=chrom, pos1=pos, ref=ref, alt=alt))
    if not parsed:
        return pd.DataFrame()
    f = pd.DataFrame(parsed)
    f['order'] = [digest([seed, r.gene, r.SNP, r.celltype]) for r in f.itertuples()]
    f = f.sort_values('order').drop_duplicates(['gene', 'SNP'])
    # Round robin cells, favour distinct genes; never use GI outputs.
    chosen, seen = [], set()
    bycell = [iter(g.to_dict('records')) for _, g in f.groupby('celltype', sort=True)]
    while len(chosen) < target:
        added = False
        for rows in bycell:
            for r in rows:
                if r['gene'] not in seen:
                    chosen.append(r); seen.add(r['gene']); added = True; break
            if len(chosen) == target:
                break
        if not added:
            break
    result = pd.DataFrame(chosen).drop(columns='order')
    result['status'] = 'provisional_missing_beta_reference_coverage'
    return result


def match_controls(candidates, pool):
    """Within-locus matching with fixed tolerances; no replacement or outcome use.

    Requires precomputed MAF, absolute TSS distance and 200-bp genomic GC.
    Lower-PIP controls are uncertain comparators, not causal negatives.
    """
    required = {'locus', 'maf', 'distance', 'gc', 'SNP', 'PIP'}
    if not required <= set(pool) or not required <= set(candidates):
        raise ValueError('Missing control matching covariates')
    used = set(candidates.SNP)
    selected, excluded = [], []
    for r in candidates.sort_values(['locus', 'SNP']).to_dict('records'):
        g = pool[(pool.locus == r['locus']) & (pool.PIP < .5) & ~pool.SNP.isin(used)].copy()
        g = g[(abs(g.maf-r['maf']) <= .05) & (abs(g.distance-r['distance']) <= 5000) & (abs(g.gc-r['gc']) <= .1)]
        if g.empty:
            excluded.append(dict(SNP=r['SNP'], reason='no_covariate_matched_control'))
            continue
        g['cost'] = abs(g.maf-r['maf'])/.05 + abs(g.distance-r['distance'])/5000 + abs(g.gc-r['gc'])/.1
        control = g.sort_values(['cost', 'SNP']).iloc[0].drop('cost').to_dict()
        used.add(control['SNP']); selected.append(control)
    return pd.DataFrame(selected), pd.DataFrame(excluded)
