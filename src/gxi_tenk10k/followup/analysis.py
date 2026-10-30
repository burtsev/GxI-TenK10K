"""Frozen all-pair estimands and independent-locus resampling."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata
from ..util import save,file_hash,digest
from .prompts import ARMS,GENE_ARMS
from .runner import export,freeze_analysis
SEED=20261030

def holm(p):
    a=np.array(p,float);order=np.argsort(a);v=np.maximum.accumulate(a[order]*np.arange(len(a),0,-1));out=np.empty(len(a));out[order]=np.minimum(v,1);return out.tolist()
def correct(y,x,eps=1e-6):
    y=np.asarray(y,float);x=np.asarray(x,float)
    return (np.isfinite(y)&np.isfinite(x)&(abs(x)>eps)&(y*x>0)).astype(float)
def rho(y,x):
    ok=np.isfinite(y)&np.isfinite(x)
    if ok.sum()<3 or np.ptp(y[ok])==0 or np.ptp(x[ok])==0:return None
    return float(spearmanr(y[ok],x[ok]).statistic)
def observations(out,split,arm):
    out=Path(out);f=pd.read_parquet(out/f'cohort_{split}.parquet')
    p=export(out,split);p=p[p.arm==arm]
    if p.observation_id.duplicated().any():raise ValueError('Duplicate frozen jobs')
    return f.merge(p[['observation_id','delta','status']],on='observation_id',how='left',validate='one_to_one')
def contrasts(out,split,f):
    p=pd.read_csv(Path(out)/'preselected_cell_pairs.tsv',sep='\t');p=p[(p.split==split)&p.pair_id.isin(f.pair_id)]
    # Candidate primary; controls have their own exploratory contrast population.
    g=f[['pair_id','cell','beta','SE','delta','locus','observation_id']].copy()
    a=p.merge(g.rename(columns={'cell':'cell_a','beta':'beta_a','SE':'SE_a','delta':'delta_a','observation_id':'id_a'}),on=['pair_id','locus','cell_a'],how='left',validate='many_to_one')
    b=a.merge(g.rename(columns={'cell':'cell_b','beta':'beta_b','SE':'SE_b','delta':'delta_b','observation_id':'id_b'}),on=['pair_id','locus','cell_b'],how='left',validate='many_to_one')
    b['observed']=b.beta_a-b.beta_b;b['predicted']=np.round(b.delta_a-b.delta_b,10)
    b['precision_subset']=abs(b.observed)>2*(b.SE_a+b.SE_b)
    return b

def summary(y,x,loci,eps=1e-6,B=0):
    y=np.asarray(y,float);x=np.asarray(x,float);loci=np.asarray(loci);c=correct(y,x,eps);ok=np.isfinite(y)&np.isfinite(x);nz=ok&(abs(x)>eps)
    result={'n_all':len(y),'n_loci':len(set(loci)),'n_scored':int(ok.sum()),'n_missing':int((~ok).sum()),'n_zero':int((ok&~nz).sum()),'n_nonzero':int(nz.sum()),'n_correct':int(c.sum()),'all_pair_accuracy':float(c.mean()) if len(c) else None,'conditional_accuracy':float(c.sum()/nz.sum()) if nz.sum() else None,'spearman':rho(y,x)}
    if B and len(y):
        unique,index=np.unique(loci,return_inverse=True);groups=[np.flatnonzero(index==i) for i in range(len(unique))];rng=np.random.default_rng(SEED);acc=[];rs=[]
        for _ in range(B):
            idx=np.concatenate([groups[i] for i in rng.integers(0,len(unique),len(unique))]);acc.append(c[idx].mean());r=rho(y[idx],x[idx]);rs.append(np.nan if r is None else r)
        result['accuracy_ci95']=np.quantile(acc,[.025,.975]).tolist();result['spearman_ci95']=np.nanquantile(rs,[.025,.975]).tolist() if np.isfinite(rs).any() else None
    return result

def dev_select(out):
    out=Path(out);freeze_analysis(out);metrics={}
    raw=export(out,'development');expected=len(pd.read_parquet(out/'cohort_development.parquet'))
    if any(len(raw[raw.arm==arm])!=expected for arm in ARMS):raise ValueError('Development arm family incomplete; selection prohibited')
    for arm in ARMS:
        f=observations(out,'development',arm);f=f[f.role=='candidate'];c=contrasts(out,'development',f)
        metrics[arm]={'gene_aware':arm in GENE_ARMS,'variant':summary(f.beta,f.delta,f.locus),'cross_cell':summary(c.observed,c.predicted,c.locus)}
    base=metrics['baseline'];eligible=[]
    for arm,m in metrics.items():
        if arm=='baseline' or arm in GENE_ARMS:continue
        gain=m['cross_cell']['all_pair_accuracy']-base['cross_cell']['all_pair_accuracy'];coverage=m['cross_cell']['n_scored']/m['cross_cell']['n_all']-base['cross_cell']['n_scored']/base['cross_cell']['n_all']
        if gain>=.03-1e-12 and coverage>=-.03-1e-12:eligible.append(arm)
    selected=max(eligible,key=lambda a:(metrics[a]['cross_cell']['all_pair_accuracy'],metrics[a]['cross_cell']['spearman'] or -2,metrics[a]['variant']['all_pair_accuracy'])) if eligible else 'baseline'
    save(out/'development_arm_metrics.json',metrics)
    code={str(p):file_hash(p) for root in ['src','tests','configs','scripts'] for p in sorted(Path(root).rglob('*')) if p.is_file() and p.suffix in ('.py','.yaml','.tsv','.json')}
    save(out/'code_manifest.json',code)
    files=['code_manifest.json','frozen_protocol.json','prompt_registry.parquet','analysis_plan.json','gi_openapi.json','development_arm_metrics.json']
    lock={'selected_arm':selected,'eligible_alternatives':eligible,'gene_free_only':True,'selection_rule':'3pp cross-cell all-pair gain; <=3pp pair-coverage loss; rho then variant direction tie-break','git_commit':__import__('subprocess').check_output(['git','rev-parse','HEAD'],text=True).strip(),'frozen_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'file_hashes':{n:file_hash(out/n) for n in files}}
    target=out/'context_selection.json'
    if target.exists():
        previous=json.load(target.open());lock['frozen_utc']=previous['frozen_utc'];lock['git_commit']=previous['git_commit']
        if previous!=lock:raise ValueError('Frozen context selection differs')
    save(target,lock);return lock

def primary_permutation(f,c,R=2000):
    """Common random cell priorities within locus induce within-pair permutations.

    Each pair retains its frozen measured-cell set; pairs in the same locus share
    random priorities. Missing values and numerical zeros are permuted too.
    """
    f=f.reset_index(drop=True);y=f.beta.to_numpy(float);x=f.delta.to_numpy(float);ids={x:i for i,x in enumerate(f.observation_id)}
    ia=np.array([ids[x] for x in c.id_a]);ib=np.array([ids[x] for x in c.id_b]);cy=c.observed.to_numpy(float)
    contrast=np.round(x[ia]-x[ib],10);obsA=correct(cy,contrast).mean();obsR=rho(cy,contrast);obsB=correct(y,x).mean()
    loci,li=np.unique(f.locus,return_inverse=True);cellnames=sorted(f.cell.unique());ci=np.array([cellnames.index(v) for v in f.cell]);groups=[g.index.to_numpy() for _,g in f.groupby('pair_id',sort=False)]
    bysize={k:np.stack([g for g in groups if len(g)==k]) for k in sorted(set(map(len,groups)))}
    rng=np.random.default_rng(SEED);na=nr=nb=0;null=[]
    for _ in range(R):
        priorities=rng.random((len(loci),len(cellnames)));xp=x.copy()
        for idx in bysize.values():
            order=np.argsort(priorities[li[idx],ci[idx]],axis=1);xp[idx]=x[np.take_along_axis(idx,order,axis=1)]
        ca=np.round(xp[ia]-xp[ib],10);a=correct(cy,ca).mean();r=rho(cy,ca);flips=rng.choice([-1,1],len(loci));b=correct(y,x*flips[li]).mean()
        na+=a>=obsA-1e-12;nr+=obsR is None or (r is not None and r>=obsR-1e-12);nb+=b>=obsB-1e-12;null.append([a,r,b])
    pdir=(na+1)/(R+1);prho=(nr+1)/(R+1);pb=(nb+1)/(R+1);pa=max(pdir,prho);pa_raw=pa;pb_raw=pb
    # Null randomization is not a 50% test when zeros are present. Require the
    # prespecified >50% success floor before claiming direction skill.
    if obsA<=.5:pa=1.
    if obsB<=.5:pb=1.
    return {'direction_decision_floor':.5,'raw_p_A':pa_raw,'raw_p_variant_signflip':pb_raw,'subfloor_primary_p_rule':'set primary p to 1 when all-pair direction<=0.5','p_cell_direction':pdir,'p_cell_rho':prho,'p_A':pa,'p_B':pb,'holm_A_B':holm([pa,pb]),'replicates':R,'null':null}

def paired_comparison(a,b,ca,cb,B=3000,R=2000):
    # Frozen same row order and cell pairs, including missing and zero predictions.
    assert a.observation_id.tolist()==b.observation_id.tolist()
    def endpoint(y,xa,xb,loci):
        d=correct(y,xb)-correct(y,xa);labels,ix=np.unique(loci,return_inverse=True);s=np.bincount(ix,weights=d);n=np.bincount(ix);rng=np.random.default_rng(SEED);boot=[]
        for _ in range(B):
            sample=rng.integers(0,len(s),len(s));boot.append(s[sample].sum()/n[sample].sum())
        obs=d.mean();count=0
        for _ in range(R):count+=(s*rng.choice([-1,1],len(s))).sum()/n.sum()>=obs-1e-12
        return dict(gain=float(obs),ci95=np.quantile(boot,[.025,.975]).tolist(),p=(count+1)/(R+1),n=len(d))
    v=endpoint(a.beta,a.delta,b.delta,a.locus);c=endpoint(ca.observed,ca.predicted,cb.predicted,ca.locus)
    return {'variant':v,'cross_cell':c,'holm_p':holm([c['p'],v['p']])}

def evaluate(out):
    out=Path(out);plan=freeze_analysis(out);lock=json.load(open(out/'context_selection.json'));arms=list(dict.fromkeys(['baseline',lock['selected_arm']]))
    raw=export(out,'test');expected=len(pd.read_parquet(out/'cohort_test.parquet'))
    if any(len(raw[raw.arm==arm])!=expected for arm in arms):raise ValueError('Test incomplete; complete scoring or document an explicit partial-study amendment')
    result={'test':{},'development_all_arms':json.load(open(out/'development_arm_metrics.json')),'selection':lock,'analysis_plan':plan}
    frames={};cs={}
    for arm in arms:
        f=observations(out,'test',arm);candidate=f[f.role=='candidate'].copy();c=contrasts(out,'test',candidate);frames[arm]=candidate;cs[arm]=c
        print('Evaluating',arm,flush=True)
        primary=primary_permutation(candidate,c,plan['permutation_replicates']);null=primary.pop('null');pd.DataFrame(null,columns=['cell_accuracy','cell_rho','variant_accuracy']).to_csv(out/f'null_{arm}.tsv',sep='\t',index=False)
        m={'variant':summary(candidate.beta,candidate.delta,candidate.locus,B=3000),'cross_cell':summary(c.observed,c.predicted,c.locus,B=3000),'primary_tests':primary,'precision_contrasts':summary(c[c.precision_subset].observed,c[c.precision_subset].predicted,c[c.precision_subset].locus,B=3000),'per_cell':{cell:summary(g.beta,g.delta,g.locus) for cell,g in candidate.groupby('cell')},'epsilon_sensitivity':{str(e):{'variant':summary(candidate.beta,candidate.delta,candidate.locus,eps=e),'cross_cell':summary(c.observed,c.predicted,c.locus,eps=e)} for e in [1e-4,1e-3]}}
        control=f[f.role=='control'];m['uncertain_controls']=summary(control.beta,control.delta,control.locus,B=3000)
        dev=observations(out,'development',arm);dev=dev[dev.role=='candidate'];ok=np.isfinite(dev.delta)&np.isfinite(dev.beta);coef=np.linalg.lstsq(np.column_stack([np.ones(ok.sum()),dev.delta[ok]]),dev.beta[ok],rcond=None)[0]
        ok=np.isfinite(candidate.delta);pred=coef[0]+coef[1]*candidate.delta[ok];m['development_fit_calibration']={'intercept':coef[0],'slope':coef[1],'n_fit':int(np.isfinite(dev.delta).sum()),'test_n':int(ok.sum()),'test_MAE':float(abs(pred-candidate.beta[ok]).mean()),'test_RMSE':float(np.sqrt(np.mean((pred-candidate.beta[ok])**2))),'scale_warning':'SAIGE Poisson coefficient is not log(TPM+1); empirical calibration only'}
        paired=f[['pair_id','cell','role','delta','PIP']].pivot(index=['pair_id','cell'],columns='role',values=['delta','PIP']);a=abs(paired['delta']['candidate']);b=abs(paired['delta']['control']);valid=a.notna()&b.notna();m['candidate_control_ranking']={'n_matched_scored':int(valid.sum()),'candidate_larger_abs_delta':int((a[valid]>b[valid]).sum()),'ties':int((a[valid]==b[valid]).sum()),'control_pip_missing':int(paired['PIP']['control'].isna().sum()),'control_interpretation':'uncertain; missing PIP not zero; no definitive negatives'}
        candidate.to_parquet(out/f'evaluated_test_{arm}.parquet',index=False);c.to_parquet(out/f'contrasts_test_{arm}.parquet',index=False);result['test'][arm]=m
    if len(arms)==2:result['paired_alternative']=paired_comparison(frames['baseline'],frames[arms[1]],cs['baseline'],cs[arms[1]])
    diagnostics=export(out,'development');result['negative_diagnostics']={}
    for arm in ['swapped_gene_label','nonsense_gene_label','shuffled_cell_label']:
        if arm in set(diagnostics.arm):
            f=observations(out,'development',arm);f=f[f.observation_id.isin(diagnostics[diagnostics.arm==arm].observation_id)&(f.role=='candidate')];result['negative_diagnostics'][arm]=summary(f.beta,f.delta,f.locus,B=3000)
    save(out/'metrics.json',result);return result

def development_inference(out):
    """Report full arm family after selection; never feed p-values into choice."""
    out=Path(out);base=observations(out,'development','baseline');base=base[base.role=='candidate'];bc=contrasts(out,'development',base);rows={};p=[];names=[]
    for arm in ARMS:
        f=observations(out,'development',arm);f=f[f.role=='candidate'];c=contrasts(out,'development',f)
        m={'variant':summary(f.beta,f.delta,f.locus,B=3000),'cross_cell':summary(c.observed,c.predicted,c.locus,B=3000),'gene_aware':arm in GENE_ARMS}
        if arm!='baseline':
            m['paired_vs_baseline']=paired_comparison(base,f,bc,c);p.append(m['paired_vs_baseline']['cross_cell']['p']);names.append(arm)
        rows[arm]=m
    for arm,adjusted in zip(names,holm(p)):rows[arm]['cross_cell_improvement_holm_9_arms']=adjusted
    diagnostic=export(out,'development');neg={}
    for arm in ['swapped_gene_label','nonsense_gene_label','shuffled_cell_label']:
        ids=set(diagnostic[diagnostic.arm==arm].observation_id);f=observations(out,'development',arm);f=f[(f.role=='candidate')&f.observation_id.isin(ids)];b=base[base.observation_id.isin(ids)]
        if len(f):
            c=contrasts(out,'development',f);cbase=contrasts(out,'development',b)
            neg[arm]={'variant':summary(f.beta,f.delta,f.locus,B=3000),'cross_cell':summary(c.observed,c.predicted,c.locus,B=3000),'paired_vs_baseline':paired_comparison(b,f,cbase,c)}
    save(out/'development_inference.json',{'arms':rows,'diagnostics':neg});return rows
