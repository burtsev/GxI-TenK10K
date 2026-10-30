"""Cluster resampling preserves dependence within loci and across cells."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, theilslopes


def evaluate(frame, epsilon=1e-6, replicates=2000, permutations=1000, seed=20261029):
    needed = ['beta_alt', 'delta', 'locus', 'gene', 'cell']
    for col in needed:
        if col not in frame:
            raise ValueError('Missing evaluation field '+col)
    finite = np.isfinite(frame.beta_alt) & np.isfinite(frame.delta) & frame[needed].notna().all(axis=1)
    valid = frame[finite & (frame.beta_alt != 0) & (abs(frame.delta) > epsilon)].copy()
    out = {'n_total': len(frame), 'n_missing': int((~finite).sum()),
           'n_abstain': int((finite & ((frame.beta_alt == 0) | (abs(frame.delta) <= epsilon))).sum()),
           'n_scored': len(valid), 'n_loci': valid.locus.nunique(), 'n_genes': valid.gene.nunique(),
           'n_cells': valid.cell.nunique(), 'sign_concordance': None, 'ci95': None,
           'spearman': None, 'permutation_p': None, 'epsilon': epsilon}
    magnitude = frame[finite]
    if len(magnitude) and magnitude.beta_alt.nunique() > 1 and magnitude.delta.nunique() > 1:
        out['spearman'] = float(spearmanr(magnitude.beta_alt, magnitude.delta).statistic)
        out['robust_slope'] = float(theilslopes(magnitude.beta_alt, magnitude.delta).slope)
    out['spearman_n'] = len(magnitude)
    if valid.empty:
        return out
    valid['agree'] = (np.sign(valid.beta_alt) == np.sign(valid.delta)).astype(float)
    out['sign_concordance'] = float(valid.agree.mean())
    groups = [g for _, g in valid.groupby('locus', sort=True)]
    sizes = np.array([len(g) for g in groups])
    successes = np.array([g.agree.sum() for g in groups])
    rng = np.random.default_rng(seed)
    if len(groups) > 1:
        boots = []
        for _ in range(replicates):
            sampled = rng.integers(len(groups), size=len(groups))
            boots.append(successes[sampled].sum()/sizes[sampled].sum())
        out['ci95'] = [float(x) for x in np.quantile(boots, [0.025, .975])]
        # Locus-level random sign-flips; cells in a locus move together.
        extreme = 0
        for _ in range(permutations):
            flips = rng.choice([-1, 1], size=len(groups))
            agreement = np.where(flips == 1, successes, sizes-successes).sum()/sizes.sum()
            extreme += agreement >= out['sign_concordance']
        out['permutation_p'] = (extreme+1)/(permutations+1)
    if valid.beta_alt.nunique() > 1 and valid.delta.nunique() > 1:
        out['spearman_nonabstaining'] = float(spearmanr(valid.beta_alt, valid.delta).statistic)
    out['chance_baseline'] = .5
    out['calibration_sensitivity'] = {}
    for e in [0, 1e-6, 1e-5, 1e-4]:
        subset = frame[finite & (frame.beta_alt != 0) & (abs(frame.delta) > e)]
        out['calibration_sensitivity'][str(e)] = {'retained': len(subset), 'sign_concordance': float((np.sign(subset.beta_alt) == np.sign(subset.delta)).mean()) if len(subset) else None}
    return out


def split_loci(frame, seed):
    loci = sorted(frame.locus.unique())
    rng = np.random.default_rng(seed)
    rng.shuffle(loci)
    holdout = set(loci[:max(1, len(loci)//5)])
    return frame.locus.map(lambda x: 'holdout' if x in holdout else 'development')
