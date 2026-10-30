"""Report exact-repeat and paraphrase tests, including unsuccessful variants."""
import json,shutil
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from gxi_tenk10k.util import save,file_hash


def main():
    root=Path('results/context-repetition-20261030');out=Path('docs/runs/context-repetition-20261030');out.mkdir(parents=True,exist_ok=True)
    s=json.load(open(root/'summary.json'));d=json.load(open(root/'design.json'));m=json.load(open(root/'development_metrics.json'))
    extra=json.load(open(root/'repeat_validation_calls.json')) if (root/'repeat_validation_calls.json').exists() else {}
    round_calls=s['new_calls']+extra.get('new_calls',0)
    total_calls=extra.get('total_calls_including_preledger',s['total_calls_including_preledger'])
    styles=getSampleStyleSheet();styles['Normal'].fontName='Times-Roman';styles['Normal'].fontSize=10;styles['Normal'].leading=12
    def p(t,sty='Normal'):return Paragraph(t,styles[sty])
    def percent(x):return 'NA' if x is None else f'{100*x:.1f}%'
    story=[p('Combined, paraphrased and repeated GI descriptions','Title'),p('Exploratory follow-up | 30 October 2026'),Spacer(1,10),p('Design','Heading2'),p('We tested concatenation of all four earlier descriptions, that concatenation twice, compact text twice, and three new paraphrases each once or exactly twice. Biological facts were held constant: human peripheral-blood cell identity, unstimulated state and 10x 3-prime single-cell RNA-seq. Model, DNA, TSS, allele orientation and 40,960-bp effective geometry were unchanged. REF and ALT used identical text.'),p('The same 16 development pairs were used for selection. Correct directions divided by all 16 pairs was the objective; zero effects count as incorrect. The incumbent compact description is retained on primary ties. Spearman correlation breaks ties only among strictly better alternatives. A new 20-locus validation set excludes all previously tested loci, including the first optimization holdout, using nonoverlapping ±100-kb TSS windows. No fresh predictions were evaluated before selection froze.')]
    rows=[['Description','Correct / 16','Nonzero','Spearman rho']]+[[key,str(v['correct_all_pairs']),str(v['n_scored']),f"{v['spearman']:.3f}" if v['spearman'] is not None else 'NA'] for key,v in m.items()]
    t=Table([[p(escape(x)) for x in row] for row in rows],colWidths=[218,87,75,125]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eef4')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story += [p('Development results','Heading2'),t,p(f"Frozen winner: <b>{escape(s['winner'])}</b>. Searching nine alternatives on sixteen pairs is highly susceptible to development-set overfitting; the complete table includes unsuccessful descriptions."),PageBreak(),p('Does exact repetition help?','Heading1')]
    for name in ['combined','compact','paraphrase_assay_first','paraphrase_cell_first','paraphrase_sentence']:
        a='compact' if name=='compact' else name+'_once';b=name+'_twice'
        story.append(p(f"<b>{name}:</b> correct directions {m[a]['correct_all_pairs']} → {m[b]['correct_all_pairs']} /16; Spearman rho {m[a]['spearman']:.3f} → {m[b]['spearman']:.3f}."))
    if (root/'repeat_validation_metrics.json').exists():
        direct=json.load(open(root/'repeat_validation_metrics.json'))
        story.append(p('Direct once-versus-twice checks on new loci','Heading2'))
        for family,result in direct.items():
            q=result['paired'];story.append(p(f"<b>{family}:</b> {q['baseline_correct_all_pairs']}/20 correct once versus {q['selected_correct_all_pairs']}/20 twice; paired all-pair change {100*q['difference']:.1f} pp (95% CI {100*q['ci95'][0]:.1f} to {100*q['ci95'][1]:.1f}). These diagnostics do not re-select the frozen prompt."))
    fresh=s['fresh'];q=fresh['paired'];a=s['baseline_original'];b=s['selected_original']
    story += [p('Frozen-choice validation and original cohort','Heading2'),p(f"At the 20 new loci, baseline versus selected correct directions were {q['baseline_correct_all_pairs']}/20 versus {q['selected_correct_all_pairs']}/20. The all-pair gain was {100*q['difference']:.1f} percentage points (95% paired locus-bootstrap interval {100*q['ci95'][0]:.1f} to {100*q['ci95'][1]:.1f}). Conditional direction agreement was {percent(fresh['baseline']['sign_concordance'])} versus {percent(fresh['selected']['sign_concordance'])}; abstentions were {fresh['baseline']['n_abstain']} versus {fresh['selected']['n_abstain']}. These small samples do not establish general improvement when the paired interval includes zero."),p(f"On the already-reported 100-pair cohort, pooled Spearman rho changed from {a['spearman']:.3f} to {b['spearman']:.3f}; conditional direction agreement from {percent(a['sign_concordance'])} to {percent(b['sign_concordance'])}. All-pair correct directions changed from {s['original_paired']['baseline_correct_all_pairs']}/100 to {s['original_paired']['selected_correct_all_pairs']}/100. This is post-tuning description, not an untouched validation set."),p('Next strategies (proposed, not executed)','Heading2'),p('Use a controlled ablation matrix: gene symbol; Ensembl gene ID; both; symbol-only repetition; canonical cell name; verified Cell Ontology ID; and explicit assay fields. Add one field at a time and keep REF/ALT text identical. Gene identity may enable learned gene lookup and must be analyzed separately from sequence-only conditioning. Freeze choices on development data before using another new validation cohort.'),p('Selected description','Heading2'),p(escape(d['templates'][s['winner']])),p('Interpretation and audit','Heading2'),p('Repetition is a testable conditioning change, not a guaranteed way to improve prediction. A development gain should be distinguished from fresh-locus gain, and higher correlation from better directional accuracy. Multiple prompt searches, small validation sets, rounding/abstentions, observational eQTLs and unknown training overlap limit interpretation. The original zero-shot benchmark remains unchanged. Exact strings for all descriptions and measured results accompany this PDF.'),p(f"API preflight: the user identified themselves as the API developer and confirmed this key’s quota was lifted. The software retained a bounded cumulative ceiling of 1,600 attempts and serial requests; headers were logged, with no automatic retries. This round made {round_calls} new attempts; cumulative attempts including the initial preledger call were {total_calls}. Monetary pricing was not published and no dollar amount is claimed. Reproduce: python scripts/optimize_context_round2.py; python scripts/build_repetition_pdf.py (GI_API_KEY in environment). Data: doi:10.5281/zenodo.18221260.")]
    def footer(c,doc):c.setFont('Times-Roman',8);c.drawString(44,20,'GI × TenK10K | Exploratory prompt repetition');c.drawRightString(A4[0]-44,20,str(doc.page))
    pdf=out/'GI_repetition_paraphrasing_report.pdf';SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=44,rightMargin=44,topMargin=36,bottomMargin=36).build(story,onFirstPage=footer,onLaterPages=footer)
    for name in ['design.json','development_metrics.json','development_panel.tsv','fresh_panel.tsv','frozen_selection.json','fresh_metrics.json','summary.json','measured_predictions.tsv']:shutil.copyfile(root/name,out/name)
    for name in ['repeat_validation_design.json','repeat_validation_metrics.json','repeat_validation_calls.json','repeat_validation_predictions.tsv']:
        if (root/name).exists():shutil.copyfile(root/name,out/name)
    save(out/'manifest.json',{'files':{p.name:file_hash(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'},'scope':'Exploratory outcome-tuned description search; fresh second validation cohort; user confirmed quota lifted.'})
    print(pdf)

if __name__=='__main__':main()
