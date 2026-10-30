"""Create a reproducible full-scoring report without new GI requests."""
import json,shutil,subprocess
from pathlib import Path
from xml.sax.saxutils import escape
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from gxi_tenk10k.util import save,file_hash


def main():
    root=Path('results/full-best-description-20261030');out=Path('docs/runs/full-best-description-20261030');out.mkdir(parents=True,exist_ok=True)
    m=json.load(open(root/'metrics.json'));gate=json.load(open(root/'discovery_gate.json'));run=json.load(open(root/'execution.json'));design=json.load(open(root/'design.json'))
    f=pd.read_parquet(root/'predictions.parquet')
    from gxi_tenk10k.analysis.metrics import evaluate
    discovery={d:{category:evaluate(g) for category,g in x.groupby('category')} for d,x in f[f.is_discovery].groupby('description')}
    save(root/'discovery_metrics.json',discovery)
    styles=getSampleStyleSheet();styles['Normal'].fontName='Times-Roman';styles['Normal'].fontSize=10;styles['Normal'].leading=12
    def p(t,style='Normal'):return Paragraph(t,styles[style])
    def pct(x):return 'NA' if x is None else f'{100*x:.1f}%'
    def ci(x):return 'NA' if x is None else f'{100*x[0]:.1f}–{100*x[1]:.1f}%'
    def table(data):
        t=Table([[p(escape(str(v))) for v in row] for row in data],colWidths=[125,65,70,120,125]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eef4')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]));return t
    b=m['baseline']['candidate'];n=m['best']['candidate'];q=m['paired']['candidate']
    story=[p('Full matched study with optimized context text','Title'),p('Repeated sentence description | 30 October 2026'),p('Summary','Heading2'),p(f"The original matched cohort was fixed at 100 candidate and 100 control variant–gene pairs. Both descriptions were scored across all four blood-cell types wherever an observed beta was available: {len(f)//2} observations per description, including {b['n_total']} candidate-derived and {m['baseline']['control']['n_total']} control-derived observations. There were {design['missing_observed_betas_per_description']} unavailable variant–gene–cell estimates per description out of 800 possible; these remain exclusions, not zero effects."),p(f"Across all available candidate-derived observations, direction agreement changed from {pct(b['sign_concordance'])} (95% locus CI {ci(b['ci95'])}) to {pct(n['sign_concordance'])} (CI {ci(n['ci95'])}). Pooled Spearman rho changed from {b['spearman']:.3f} to {n['spearman']:.3f}. Correct directions per all available observations changed from {q['baseline_correct_all_pairs']}/{q['n_pairs']} to {q['best_correct_all_pairs']}/{q['n_pairs']}; paired gain {100*q['all_pair_gain']:.1f} pp (95% CI {100*q['gain_ci95'][0]:.1f} to {100*q['gain_ci95'][1]:.1f} pp). Repeated cells were clustered by locus. These are prompt-tuned, observational results, not an independent zero-shot replication."),p('Original discovery-cell pilot repeated','Heading2')]
    rows=[['Condition','All','Nonzero','Sign agreement / 95% CI','Pooled Spearman rho']]
    for d in ['baseline','best']:
        for category in ['candidate','control']:
            x=discovery[d][category];rows.append([d+' '+category,x['n_total'],x['n_scored'],pct(x['sign_concordance'])+' / '+ci(x['ci95']),f"{x['spearman']:.3f}"])
    story += [table(rows),p('Full scoring across available cells','Heading2')]
    rows=[['Condition','All','Nonzero','Sign agreement / 95% CI','Pooled Spearman rho']]
    for d in ['baseline','best']:
        for category in ['candidate','control']:
            x=m[d][category];rows.append([d+' '+category,x['n_total'],x['n_scored'],pct(x['sign_concordance'])+' / '+ci(x['ci95']),f"{x['spearman']:.3f}"])
    story += [table(rows),PageBreak(),p('Methods, scope and interpretation','Heading1'),p('The expansion gate required both discovery-candidate correct/all-pair count and Spearman correlation to improve. The frozen repeated description passed this gate. No further prompt or variant selection used the full-scoring results. The source cohorts, control matches, allele orientation, TSSs, model and scored windows were identical across descriptions.'),p('Chosen description','Heading2'),p(escape(design['best_template'])),p('Cell specificity','Heading2')]
    for d in ['baseline','best']:
        x=m[d]['cell_contrasts'];story.append(p(f"<b>{d}:</b> {x['n_total']} pairwise contrasts, {x['n_scored']} non-abstaining; sign agreement {pct(x['sign_concordance'])} (95% CI {ci(x['ci95'])}); Spearman rho {x['spearman']:.3f}. Contrasts share variants, genes and donors; covariance of observed cell effects remains unknown, so these are descriptive model comparisons, not tests of biological heterogeneity."))
    contrast_pair=m['paired']['cell_contrasts']
    story.append(p(f"The paired cell-contrast gain in correct predictions per all contrasts was {100*contrast_pair['all_pair_gain']:.1f} pp (95% locus CI {100*contrast_pair['gain_ci95'][0]:.1f} to {100*contrast_pair['gain_ci95'][1]:.1f} pp). Individual above-chance concordance does not by itself demonstrate improvement over the baseline."))
    story += [p('Limitations and comparison','Heading2'),p('The text was selected using development outcomes in earlier prompt searches. That original-cohort discovery gain was partly in-sample. The prior new-locus check gave 10/20 correct for both original and selected text; thus independent improvement was not established before this expansion. Full scoring adds cell contexts to the same selected loci, rather than adding an independent cohort. Candidate-derived comparison-cell associations need not have PIP≥0.9 or be significant outside the discovery cell. Controls remain uncertain comparators; absent PIP is not zero.'),p('The model g0-expression-8192 used 40,960-bp oriented inputs and exact matched REF/ALT effective windows, as documented in the earlier geometry amendment. Minus-strand edits preserve ALT-minus-REF effect orientation. The initial 81,920-bp geometry was not restored. Scores are log(TPM+1), while SAIGE-QTL betas are count-model effects. GI training overlap and an immutable deployed weight revision remain unresolved. Gene-window clusters approximate LD. No causal validation is implied.'),p('Cuomo et al., Nature (2026), doi:10.1038/s41586-026-10721-2, reports PromoterAI median correlation 0.66 on lead eQTLs. Our pooled Spearman coefficients use different variants, selection and aggregation; they cannot rank the models without running them on matched inputs and observations.'),p('Execution and reproducibility','Heading2'),p(f"Completed all available cells: {run['completed_all_available_cells']}. Stored prediction records: {run['completed_records']}. New requests this repeat/expansion: {run['new_requests']}; cumulative prediction attempts including the initial preledger control: {run['total_requests_including_initial_preledger']}. The API developer/user confirmed quota lifted for this key. Requests remained serial with a 4,000-attempt software cap, no automatic retries and logged rate headers. No monetary price was available. Small measured tables, source hash, all missing-data exclusions, both descriptions and metrics accompany this report. Original benchmark outputs remain unchanged."),p('Reproduce with GI_API_KEY exported: python scripts/full_best_description.py; python scripts/report_full_best.py. Model calls resume from the content-addressed cache. Data: Zenodo 18221260; GRCh38/GENCODE v44. Offline testing result is in test_results.txt.')]
    def footer(c,doc):c.setFont('Times-Roman',8);c.drawString(44,20,'GI × TenK10K | Full paired scoring');c.drawRightString(A4[0]-44,20,str(doc.page))
    path=out/'GI_full_best_description_report.pdf';SimpleDocTemplate(str(path),pagesize=A4,leftMargin=44,rightMargin=44,topMargin=36,bottomMargin=36).build(story,onFirstPage=footer,onLaterPages=footer)
    # Keep a small, sequence-free table; raw response cache and large source pool stay ignored.
    columns=['gene','SNP','cell','category','locus','discovery_cell','is_discovery','beta_alt','SE','PIP','delta','ref_score','alt_score','description','model','scored_start','scored_end','ref_sequence_hash','alt_sequence_hash','source_file']
    f[columns].to_csv(out/'measured_predictions.tsv',sep='\t',index=False,na_rep='NA')
    for name in ['design.json','execution.json','api_preflight.json','discovery_gate.json','discovery_metrics.json','metrics.json','quality_audit.json','exclusions.tsv','baseline_cell_contrasts.tsv','best_cell_contrasts.tsv']:shutil.copyfile(root/name,out/name)
    shutil.copyfile('/tmp/gxi-full-best-tests.txt',out/'test_results.txt')
    save(out/'manifest.json',{'files':{x.name:file_hash(x) for x in out.iterdir() if x.is_file() and x.name!='manifest.json'},'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'code_hashes':{name:file_hash(name) for name in ['scripts/full_best_description.py','scripts/report_full_best.py']},'scope':'Paired prompt-tuned repeat and all-available-cell scoring; no new zero-shot or independent validation claim.'})
    print(path)

if __name__=='__main__':main()
