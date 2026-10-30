"""Standalone figures from measured predictions; no synthetic scores."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def figures(frame, metrics, directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    candidates=frame[frame.category=='candidate']
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for cell,g in candidates.groupby('cell'):
        axes[0,0].scatter(g.beta_alt,g.delta,label=cell,s=16,alpha=.7)
    axes[0,0].axhline(0,color='gray',lw=.5);axes[0,0].axvline(0,color='gray',lw=.5)
    axes[0,0].set(xlabel='Observed ALT beta (SAIGE-QTL)',ylabel='Predicted ALT − REF log(TPM+1)',title='Amended geometry: candidate effects')
    axes[0,0].legend(fontsize=7)
    labels=[];values=[];lower=[];upper=[]
    for cell,m in metrics['candidate_by_cell'].items():
        if m['sign_concordance'] is None:continue
        labels.append(cell);values.append(m['sign_concordance'])
        ci=m['ci95'] or [m['sign_concordance']]*2
        lower.append(max(0,m['sign_concordance']-ci[0]));upper.append(max(0,ci[1]-m['sign_concordance']))
    axes[0,1].errorbar(values,range(len(labels)),xerr=[lower,upper],fmt='o')
    axes[0,1].set_yticks(range(len(labels)),labels,fontsize=8)
    axes[0,1].axvline(.5,color='gray',ls='--')
    axes[0,1].set(xlim=(0,1),xlabel='Sign agreement; locus-bootstrap 95% CI',title='Non-abstaining candidates by cell')
    for category,g in frame[frame.category.isin(['candidate','control'])].groupby('category'):
        axes[1,0].hist(g.delta,bins=30,histtype='step',label=category)
    axes[1,0].legend();axes[1,0].set(xlabel='Predicted ALT − REF',ylabel='Pairs',title='Measured effect distribution')
    axes[1,1].scatter(candidates.ref_score,candidates.alt_score,s=15,alpha=.6)
    limits=[float(min(candidates.ref_score.min(),candidates.alt_score.min())),float(max(candidates.ref_score.max(),candidates.alt_score.max()))]
    axes[1,1].plot(limits,limits,color='gray',ls='--')
    axes[1,1].set(xlabel='Reference log(TPM+1)',ylabel='Alternative log(TPM+1)',title='Prediction scale and sensitivity')
    fig.savefig(directory/'benchmark.png',dpi=140);plt.close(fig)

    valid=candidates[(abs(candidates.delta)>1e-6)&(candidates.beta_alt!=0)].copy()
    valid['agree']=(np.sign(valid.delta)==np.sign(valid.beta_alt)).astype(int)
    groups=valid.groupby('locus').agree.agg(['sum','count'])
    rng=np.random.default_rng(20261029)
    flips=rng.choice([False,True],size=(1000,len(groups)))
    null=np.where(flips,groups['sum'].to_numpy(),(groups['count']-groups['sum']).to_numpy()).sum(axis=1)/len(valid)
    fig,axes=plt.subplots(1,2,figsize=(9,3),layout='constrained')
    axes[0].hist(null,bins=25,color='steelblue',alpha=.7)
    axes[0].axvline(valid.agree.mean(),color='black',label='Observed')
    axes[0].legend();axes[0].set(xlabel='Locus sign-flip concordance',ylabel='Monte Carlo replicates',title='Null distribution (display seed fixed)')
    pairs=frame[frame.category=='control'].merge(candidates[['gene','delta']],on='gene',suffixes=('_control','_candidate'))
    axes[1].hist(abs(pairs.delta_candidate)-abs(pairs.delta_control),bins=25)
    axes[1].axvline(0,color='gray',ls='--')
    axes[1].set(xlabel='|Candidate delta| − |Control delta|',ylabel='Matched pairs',title='Uncertain-comparator priority proxy')
    fig.savefig(directory/'diagnostics.png',dpi=140);plt.close(fig)
