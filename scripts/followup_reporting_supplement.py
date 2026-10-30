"""Post-analysis reporting supplement; frozen primary analysis remains unchanged.

Adds uncertainty to the prespecified matched-control ranking summaries and makes
editorial distinctions between the published decision gate and raw randomization
results. No cohorts, contexts, predictions or primary metrics are changed.
"""
import json,shutil,re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from gxi_tenk10k.followup.analysis import observations,summary
from gxi_tenk10k.followup.binding import validate_selection
from gxi_tenk10k.followup.report import build_pdf,pct,number
from gxi_tenk10k.util import save,file_hash

def main(out=Path('results/followup-20261030')):
    out=Path(out);validate_selection(out);quality=json.load(open(out/'qa.json'));assert quality['failed']==0
    original=json.load(open(out/'report_manifest.json'));assert file_hash(original['report_pdf'])==original['sha256']
    m=json.load(open(out/'metrics.json'));b=m['test']['baseline'];v=b['variant'];c=b['cross_cell'];counts=quality['cohort_counts']['test'];plan=m['analysis_plan'];f=observations(out,'test','baseline')
    keys=['pair_id','cell','locus'];cols=keys+['delta','PIP']
    paired=f[f.role=='candidate'][cols].merge(f[f.role=='control'][cols],on=keys,suffixes=('_candidate','_control'),validate='one_to_one')
    d=abs(paired.delta_candidate)-abs(paired.delta_control);valid=np.isfinite(d);non_tied=valid&(abs(d)>plan['epsilon'])
    allrank=summary(np.ones(valid.sum()),d[valid],paired.locus[valid],B=plan['bootstrap_replicates'])
    conditional=summary(np.ones(non_tied.sum()),d[non_tied],paired.locus[non_tied],B=plan['bootstrap_replicates'])
    known=paired.PIP_candidate.notna()&paired.PIP_control.notna();pipdiff=paired.PIP_candidate-paired.PIP_control
    ordered=known&(pipdiff!=0)&valid
    piprank=summary(np.sign(pipdiff[ordered]),d[ordered],paired.locus[ordered],B=plan['bootstrap_replicates'])
    supplement={'status':'reporting supplement completed after frozen primary analysis; no primary parameters or data changed','bootstrap_replicates':plan['bootstrap_replicates'],'seed':plan['seed'],'ranking':{'n_all':int(valid.sum()),'wins':int((d[valid]>plan['epsilon']).sum()),'ties':int((abs(d[valid])<=plan['epsilon']).sum()),'losses':int((d[valid]<-plan['epsilon']).sum()),'all_pair_win_fraction':allrank['all_pair_accuracy'],'all_pair_ci95':allrank['accuracy_ci95'],'non_tied_n':int(non_tied.sum()),'non_tied_win_fraction':conditional['all_pair_accuracy'],'non_tied_ci95':conditional['accuracy_ci95']},'reported_pip_ordering':{'n_both_pip_reported':int(known.sum()),'n_equal_pip':int((known&(pipdiff==0)).sum()),'metrics':piprank,'interpretation':'Reported posterior-probability ordering only; unknown PIP omitted explicitly, never imputed zero. Not experimental causal truth.'},'observed_expression_stratification':{'status':'not_run','reason':'The authorized association summary tables do not provide observed gene-expression means; donor-level expression is not reconstructed.'},'timing_utc':datetime.now(timezone.utc).isoformat(),'frozen_primary_metrics_sha256':file_hash(out/'metrics.json'),'supplement_source_sha256':file_hash(__file__)}
    save(out/'reporting_supplement.json',supplement)
    text=(out/'report.md').read_text();preamble,*parts=text.split('\n## ');sections=[];tail=[]
    for part in parts:
        name,body=part.split('\n\n',1)
        if name in ['Entire development family','Legacy and independent study comparison']:tail.append((name,body));continue
        if name in ['Secondary matched ranking uncertainty','Discussion: signal, abstention and the decision gate']:continue
        sections.append((name,body.strip()))
    abstract=(f"We assessed GI predictions against TenK10K summary effects using 250 development loci and a frozen test of {counts['variant_gene_pairs']:,} variant–gene pairs across {counts['genes']:,} genes and {counts['locus_clusters']:,} physical locus clusters, with matched uncertain controls and 28 immune populations. Previously scored loci were excluded using the recovered request registry. No gene-free description alternative qualified in development, so baseline alone was tested. Variant direction was correct in {v['n_correct']:,}/{v['n_all']:,} observations: {pct(v['all_pair_accuracy'])} (95% locus CI {pct(v['accuracy_ci95'][0])}–{pct(v['accuracy_ci95'][1])}) with zeros counted incorrect, and {pct(v['conditional_accuracy'])} among {v['n_nonzero']:,} nonzero predictions. Cross-cell direction was correct in {c['n_correct']:,}/{c['n_all']:,} contrasts: {pct(c['all_pair_accuracy'])} (CI {pct(c['accuracy_ci95'][0])}–{pct(c['accuracy_ci95'][1])}), or {pct(c['conditional_accuracy'])} among {c['n_nonzero']:,} nonzero contrasts. Cross-cell Spearman rho was {number(c['spearman'])} (CI {number(c['spearman_ci95'][0])}–{number(c['spearman_ci95'][1])}). Raw clustered randomizations and the published all-pair decision gate are distinguished below. The assessment does not establish independence from GI training data or provide a matched head-to-head comparison with PromoterAI.")
    dev=m['development_all_arms'];base_accuracy=dev['baseline']['cross_cell']['all_pair_accuracy']
    dev_statement=("None of the other development descriptions exceeded baseline's point cross-cell accuracy." if all(r['cross_cell']['all_pair_accuracy']<=base_accuracy for arm,r in dev.items() if arm!='baseline') else "No gene-free alternative met the prespecified development gain/coverage rule.")
    interpretation=(f"The distinction between statistical signal and the prespecified utility gate is material. Variant-effect rho was {number(v['spearman'])} (95% locus CI {number(v['spearman_ci95'][0])}–{number(v['spearman_ci95'][1])}); the cross-cell estimate and its interval are reported above. Raw zero-preserving randomization p-values were {b['primary_tests']['raw_p_variant_signflip']:.4g} for variant direction and {b['primary_tests']['raw_p_A']:.4g} for the joint cell direction/rho statistic. The gated co-primary decision values were {b['primary_tests']['p_A']} and {b['primary_tests']['p_B']}, using the rule published before test inference. Small raw null p-values indicate orientation information even when the all-pair gate fails; a gate failure must not be interpreted as proof that predictions contain no information. Zero predictions accounted for {v['n_zero']:,}/{v['n_all']:,} variant observations and {c['n_zero']:,}/{c['n_all']:,} cell contrasts. The conditional percentages expose these abstentions rather than replacing the primary denominator. {dev_statement} This is not a held-out comparison of the unselected alternatives, because baseline alone qualified for test.")
    rank=supplement['ranking'];ranking=(f"Candidate |delta| exceeded its matched uncertain control in {rank['wins']:,}/{rank['n_all']:,} comparisons, with {rank['ties']:,} ties and {rank['losses']:,} losses. Strict all-pair success was {pct(rank['all_pair_win_fraction'])} (95% locus CI {pct(rank['all_pair_ci95'][0])}–{pct(rank['all_pair_ci95'][1])}); conditional on a non-tie, it was {pct(rank['non_tied_win_fraction'])} among {rank['non_tied_n']:,} comparisons (CI {pct(rank['non_tied_ci95'][0])}–{pct(rank['non_tied_ci95'][1])}). This is a prioritization proxy, not causal classification. Continuous reported-PIP ordering was assessable only where both probabilities were released; its full counts and locus intervals are supplied in reporting_supplement.json. Missing PIP was never assigned zero.")
    pips=supplement['reported_pip_ordering'];pm=pips['metrics']
    ranking+=f" Both PIPs were reported for {pips['n_both_pip_reported']:,} comparisons across {pm['n_loci']:,} clusters. Predicted magnitude ordering agreed in {pm['n_correct']:,}/{pm['n_all']:,}: {pct(pm['all_pair_accuracy'])} (95% locus CI {pct(pm['accuracy_ci95'][0])}–{pct(pm['accuracy_ci95'][1])}), with {pm['n_zero']:,} predicted ties. This compares released posterior probabilities, not experimental causal labels."
    revised=[]
    for name,body in sections:
        if name=='Abstract':body=abstract
        if name=='Methods: context selection and statistics':
            body=body.replace('We used 3,000 physical-locus bootstrap resamples and 2,000 clustered randomizations.','Co-primary and development-family intervals used 3,000 physical-locus bootstrap resamples; exploratory per-cell, family and covariate/pathway strata used 2,000. Clustered randomizations used 2,000 replicates.')
        if name=='Results: controls and diagnostics':
            body=re.sub(r' Matched magnitude ranking counts: \{.*?\}\.',f" Paired ranking used {b['candidate_control_ranking']['n_matched_scored']:,} scored comparisons; control PIP was unreported in {b['candidate_control_ranking']['control_pip_missing']:,} comparisons.",body)
        if name=='Results: independent test':
            body=body.replace('Holm-adjusted co-primary p-values were A=1.000 and B=1.000.','The preregistered all-pair gate was not crossed, giving Holm-adjusted gated decision values A=1 and B=1. Raw randomization evidence is distinguished from that gate below.')
        revised.append((name,body))
        if name=='Results: controls and diagnostics':revised.append(('Secondary matched ranking uncertainty',ranking))
        if name=='Discussion: relation to the Nature study':revised.append(('Discussion: signal, abstention and the decision gate',interpretation))
        if name=='Discussion: limitations':
            n,t=revised[-1];t=t.replace(' Observed gene-expression means were unavailable in the permitted association summaries, so observed-expression-level stratification was not run and no donor-level expression was inferred.','');revised[-1]=(n,t+' Observed gene-expression means were unavailable in the permitted association summaries, so observed-expression-level stratification was not run and no donor-level expression was inferred.')
    text=preamble.strip()+'\n\n'+'\n\n'.join('## '+n+'\n\n'+t for n,t in revised+tail)+'\n'
    (out/'report.md').write_text(text)
    pdf=final_pdf(out,preamble.splitlines()[0].removeprefix('# '),datetime.now(timezone.utc).strftime('%d %B %Y'),revised)
    original.update(sha256=file_hash(pdf),bytes=pdf.stat().st_size,editorial_supplement='reporting_supplement.json',editorial_note='Prose clarifies raw randomization versus gated decisions and secondary ranking uncertainty. Frozen primary metrics unchanged.',editorial_source_sha256=file_hash(__file__))
    save(out/'report_manifest.json',original)
    dest=Path('docs/runs/followup-20261030')
    for source,target in [('report.md','summary.md'),(pdf.name,pdf.name),('reporting_supplement.json','reporting_supplement.json'),('report_manifest.json','report_manifest.json'),('offline_tests.txt','offline_tests.txt')]:shutil.copy2(out/source,dest/target)
    shutil.copy2(out/'figures/flow.png',dest/'figures/flow.png')
    print(json.dumps({'report':str(pdf),'ranking':rank,'reported_PIP_pairs':int(known.sum()),'sha256':original['sha256']},indent=2))

def final_pdf(out,title,date,sections):
    """Final manuscript layout; headings stay with text, flow branches correctly."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    from reportlab.platypus import BaseDocTemplate,PageTemplate,Frame,Paragraph,Spacer,PageBreak,NextPageTemplate,Image
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.enums import TA_JUSTIFY,TA_LEFT
    from xml.sax.saxutils import escape
    counts=json.load(open(out/'cohort_description.json'));flow=pd.read_csv(out/'cohort_flow.tsv',sep='\t').set_index('stage').pairs
    fig,ax=plt.subplots(figsize=(10,3.6));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    def box(x,y,w,h,text):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015',facecolor='#edf3f7',edgecolor='#476b86'));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=10)
    box(.02,.32,.23,.32,f"Candidate reserve\n{flow['candidate_reserve']:,} pairs")
    box(.34,.32,.23,.32,f"Matched complete\n{flow['matched_complete']:,} pairs")
    box(.70,.62,.27,.28,f"Development\n{counts['development']['variant_gene_pairs']:,} pairs / {counts['development']['locus_clusters']:,} clusters")
    box(.70,.15,.27,.32,f"Frozen test\n{counts['test']['variant_gene_pairs']:,} pairs / {counts['test']['locus_clusters']:,} clusters\n{counts['test']['genes']:,} genes; 28 populations")
    for a,b in [((.26,.48),(.32,.48)),((.59,.53),(.68,.76)),((.59,.43),(.68,.31))]:ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#476b86','lw':1.5})
    ax.text(.5,.98,'Historical 167 genes and connected physical TSS clusters excluded before sampling',ha='center',va='top',fontsize=10)
    unused=int(flow['matched_complete']-counts['development']['variant_gene_pairs']-counts['test']['variant_gene_pairs'])
    ax.text(.5,.03,f'{unused:,} matched reserve pairs not selected; one uncertain control per frozen candidate',ha='center',fontsize=9)
    fig.tight_layout();fig.savefig(out/'figures/flow.png',dpi=180);plt.close(fig)
    W,H=A4;M=42;gap=18;cw=(W-2*M-gap)/2
    body=ParagraphStyle('body',fontName='Times-Roman',fontSize=9.3,leading=11.4,alignment=TA_JUSTIFY,spaceAfter=6)
    head=ParagraphStyle('head',parent=body,fontName='Times-Bold',fontSize=10.5,leading=12.2,spaceBefore=6,keepWithNext=True)
    small=ParagraphStyle('small',parent=body,fontSize=8.5,leading=10.4,alignment=TA_LEFT)
    title_style=ParagraphStyle('title',fontName='Times-Bold',fontSize=20,leading=23,spaceAfter=10)
    def para(text,style=body):return Paragraph(escape(text).replace('\n','<br/>').replace('≥','&gt;=').replace('≤','&lt;='),style)
    def footer(can,doc):
        can.saveState();can.setFont('Times-Roman',7.5);can.drawString(M,22,'GI × TenK10K | Independent follow-up | Research manuscript, not peer reviewed');can.drawRightString(W-M,22,str(doc.page));can.restoreState()
    pdf=out/'GI_TenK10K_followup_manuscript.pdf';doc=BaseDocTemplate(str(pdf),pagesize=A4,title=title,author='GI × TenK10K benchmark project',pageCompression=1)
    full=Frame(M,42,W-2*M,H-84,id='full',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
    columns=[Frame(M,42,cw,H-84,id='l',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0),Frame(M+cw+gap,42,cw,H-84,id='r',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)]
    doc.addPageTemplates([PageTemplate(id='front',frames=[full],onPage=footer),PageTemplate(id='text',frames=columns,onPage=footer),PageTemplate(id='figures',frames=[full],onPage=footer)])
    story=[para(title,title_style),para('GI × TenK10K benchmark project | '+date),para('Abstract',head),para(sections[0][1]),Spacer(1,8),Image(str(out/'figures/flow.png'),width=W-2*M,height=(W-2*M)*.36)]
    rest=list(sections[1:])
    # Fill the front with introductory text when it fits, without forcing an overflow.
    for name in ['Introduction','Methods: data and isolation']:
        match=next(((n,t) for n,t in rest if n==name),None)
        if not match:continue
        proposed=story+[para(match[0],head),para(match[1])]
        height=sum(item.wrap(W-2*M,10000)[1]+getattr(item,'spaceBefore',0)+getattr(item,'spaceAfter',0) for item in proposed)
        if height<700:story=proposed;rest.remove(match)
    story.extend([NextPageTemplate('text'),PageBreak()])
    for name,text in rest:
        story.append(para(name,head));style=small if name in ('Reproducibility and availability','References') else body
        story.extend(para(t,style) for t in text.split('\n\n'))
    story.extend([NextPageTemplate('figures'),PageBreak(),para('Figure 1. Independent test results',head),Image(str(out/'figures/primary.png'),width=W-2*M,height=W-2*M),para('A: fixed candidate cell contrasts; units differ across axes. B: all-pair variant direction by population. C: all-pair and conditional success with separate zero rate. D: within-locus shuffled-cell null; red line is the observed statistic. All model zeros remain in the all-pair denominator.'),PageBreak(),para('Figure 2. Development, calibration and geometry',head),Image(str(out/'figures/secondary.png'),width=W-2*M,height=(W-2*M)*.7),para('A: complete development family, including the separate gene-aware track, with locus intervals. B: test evaluation of development-fitted calibration. C: family sample sizes. D: longer-window apparent deltas; red points fail the identical-effective-span rule and are not valid biological effects.')])
    doc.build(story);return pdf

if __name__=='__main__':main()
