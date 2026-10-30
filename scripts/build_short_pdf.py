"""Build a three-page manuscript-style PDF from committed benchmark outputs."""
import json
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.pdfgen import canvas
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, PageBreak, NextPageTemplate, Image, Table, TableStyle, FrameBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/runs/live-20261029/short_report'
DATA=ROOT/'docs/runs/live-20261029'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    m=json.loads((DATA/'metrics.json').read_text());c=m['candidate'];ctrl=m['control']
    h=json.loads((DATA/'cell_contrast_metrics.json').read_text())
    style=ParagraphStyle('body',fontName='Times-Roman',fontSize=10,leading=12.3,alignment=TA_JUSTIFY,spaceAfter=6)
    head=ParagraphStyle('head',fontName='Times-Bold',fontSize=11,leading=13,spaceBefore=6,spaceAfter=5)
    small=ParagraphStyle('small',fontName='Times-Roman',fontSize=8.4,leading=10.4,spaceAfter=5)
    title=ParagraphStyle('title',fontName='Times-Bold',fontSize=20,leading=23,spaceAfter=8)
    abstract=ParagraphStyle('abstract',parent=style,fontSize=9.5,leading=11.6)
    W,H=A4;M=44;gap=18;cw=(W-2*M-gap)/2
    def p(t,s=style):return Paragraph(t,s)
    pdf=OUT/'GI_TenK10K_short_report.pdf'
    abstract_text=("<b>Abstract.</b> We assessed zero-shot Genomic Intelligence predictions against TenK10K immune-cell cis-eQTL summaries. "
        "A matched pilot comprised 100 high-PIP candidate variant–gene pairs and 100 uncertain controls across four blood-cell types. "
        "Candidate direction agreement was <b>67.5%</b> (95% locus-bootstrap CI, <b>56.4–77.9%</b>) among 77 non-abstaining effects; 23 had near-zero predicted changes. "
        "Pooled effect correlation was Spearman rho=0.370. Controls showed 51.7% direction agreement. Cross-cell contrasts did not show directional agreement above chance. "
        "The published Nature study reports a PromoterAI median correlation of 0.66, an unmatched comparison. "
        "Geometry and matching adjustments make these results exploratory; the unchanged original primary protocol was not completed.")
    def footer(can,doc):
        can.saveState();can.setStrokeColor(colors.HexColor('#b1b7bf'));can.line(M,32,W-M,32)
        can.setFont('Times-Roman',8);can.setFillColor(colors.HexColor('#555555'))
        can.drawString(M,20,'GI × TenK10K | Research report | 29 October 2026')
        can.drawRightString(W-M,20,str(doc.page));can.restoreState()
    def first(can,doc):
        footer(can,doc)
        text='Zero-shot prediction of immune-cell eQTL effects with Genomic Intelligence'
        a=p(text,title);_,height=a.wrap(W-2*M,200);a.drawOn(can,M,H-M-height)
        y=H-M-height-17;can.setFont('Times-Italic',9);can.drawString(M,y,'A reproducible pilot benchmark with a comparison to the TenK10K Nature study')
        y-=17;can.setStrokeColor(colors.HexColor('#374b65'));can.line(M,y,W-M,y)
        a=p(abstract_text,abstract);_,height=a.wrap(W-2*M,200);a.drawOn(can,M,y-height-9)
    first_top=H-44-190
    first_frames=[Frame(M,44,cw,first_top-44,id='left1',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0),Frame(M+cw+gap,44,cw,first_top-44,id='right1',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)]
    regular=[Frame(M,44,cw,H-88,id='left',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0),Frame(M+cw+gap,44,cw,H-88,id='right',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)]
    full=[Frame(M,44,W-2*M,H-88,id='full',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)]
    doc=BaseDocTemplate(str(pdf),pagesize=A4,title='GI × TenK10K: zero-shot eQTL benchmark',author='GI × TenK10K computational benchmark',subject='Executed pilot results and comparison with Cuomo et al., Nature (2026)',pageCompression=1)
    doc.addPageTemplates([PageTemplate(id='first',frames=first_frames,onPage=first),PageTemplate(id='text',frames=regular,onPage=footer),PageTemplate(id='figure',frames=full,onPage=footer)])
    story=[]
    def section(name,texts):
        story.append(p(name,head));story.extend(p(t) for t in texts)
    section('Introduction',[
        'Sequence models can be tested against naturally occurring allelic effects without fitting to the target associations. TenK10K provides such a testbed. We used its released summaries, rather than donor-level expression or genotype data, to assess effect direction and cross-cell differences. The source preprint [1] and its Nature version of record [2] are treated as versions of the same study, not independent replications.'
    ])
    section('Methods',[
        'We audited Zenodo record 18221260 [3], verified supplied MD5 checksums and independent SHA256 hashes, and scanned 63,986,931 full-summary rows in B-naive, CD4-naive T, CD14 classical-monocyte and CD16 nonclassical-monocyte populations. Candidate eligibility required gene-level ACAT Benjamini–Hochberg q≤0.05, SuSiE PIP≥0.9, a nonpalindromic biallelic SNV and an unambiguous canonical GENCODE v44 TSS. Distinct genes were sampled reproducibly using seed 20261029.',
        'GRCh38 REF alleles were checked against the release FASTA. Edits occurred on genomic-plus sequence before reverse complementation for minus-strand genes. Observed effects were SAIGE-QTL BETA for Allele2=ALT. The GI model g0-expression-8192 received identical biological-context descriptions for REF and ALT; the predicted effect was ALT minus REF in log(TPM+1).',
        'Controls shared gene and cell, with MAF difference≤0.05, absolute TSS-distance difference≤5 kb and 201-bp reference-GC difference≤0.1. Their PIP was low (&lt;0.5) or unreported. Missing PIP was not replaced by zero. Directional analysis abstained at |delta|≤10<super>−6</super>. We used 2,000 bootstrap resamples and 1,000 sign-flip null replicates clustered by overlapping ±100-kb TSS windows. These are locus proxies, not measured LD blocks.'
    ])
    section('Results',[
        'A two-variant smoke and a 20-candidate trial preceded pilot expansion. The final pilot contained 100 candidates and 100 matched controls. Candidate agreement was 52/77 (67.5%; 95% CI 56.4–77.9%); 23/100 predictions abstained, with no missing observed effects. Across all 100 candidates, including zero predictions, Spearman rho was 0.370. The locus sign-flip P value was 0.001 (Monte Carlo resolution).',
        'Controls agreed in 31/60 non-abstaining comparisons (51.7%; 95% CI 39.3–64.4%); 40/100 abstained. Their pooled rho was −0.006. Candidate |delta| exceeded matched-control |delta| in 60/76 non-tied pairs (78.9%; 95% CI 69.6–87.8%); 24 pairs tied within epsilon. Because 98 controls had unreported PIP, this is an uncertain-comparator priority proxy, not causal classification accuracy.'
    ])
    story.insert(5,FrameBreak())
    story.extend([NextPageTemplate('text'),PageBreak()])
    section('Cell specificity and protocol adjustments',[
        'For the first 20 matched candidates, comparison-cell data were available for 63/80 possible variant–gene–cell observations; all 17 missing comparisons were recorded. These yielded 81 pairwise cell contrasts, of which 71 were non-abstaining across 17 locus clusters. Agreement was 32/71 (45.1%; 95% CI 27.8–59.0%); rho across all contrasts was −0.040. This does not demonstrate prediction of cell-specific allelic differences.',
        'The original 81,920-bp inputs produced REF/ALT effective-window endpoints differing by 4–7 bp at the real smoke loci. Before the trial, we switched to 40,960-bp inputs, which were scored in full and identically for both alleles. CETP changed predicted effect sign between geometries. Therefore these results apply to the amended input, not an unchanged confirmatory primary endpoint.',
        'An initial pass produced 100 candidates but only 76 matched controls. We then checked control availability before inference, retaining the same source order and matching tolerances without beta- or delta-based eligibility. The earlier cohort remains available as a sensitivity analysis: 73.1% agreement among 78 non-abstaining candidates (95% CI 62.2–83.3%). The adjustment was outcome-independent in code but followed an unblinded first pass.'
    ])
    section('Comparison with the Nature paper',[
        'Cuomo et al. [2] analyse 1,925 individuals and 28 immune-cell types, identifying 17,674 common eGenes. They report median correlation 0.66 between PromoterAI scores and each gene’s lead-eQTL effects across cell types. Their main text does not label this coefficient as Pearson or Spearman. Our pooled Spearman rho=0.370 uses a different model, variant selection and aggregation. The numerical difference cannot establish relative model performance: PromoterAI was not rerun on this matched cohort.',
        'The Nature study finds about half of effects limited to one or a few cell types and uses LD- and power-aware sharedness criteria. Our small, noisy contrast analysis evaluates model predictions rather than population-level discovery. Its near-chance result neither reproduces that biological specificity nor contradicts it. A shared-data, shared-metric comparison is needed.'
    ])
    story.append(FrameBreak())
    section('Discussion',[
        'The directional signal among non-abstaining candidates is encouraging, but abstention and ascertainment matter: 23% of candidate predictions and 40% of controls were effectively zero. Per-cell concordance ranged from 41.7% in B-naive cells to 82.6% in CD16 monocytes, with small denominators and broad intervals. These strata are descriptive and do not establish cell-specific superiority.',
        'Matching limits trivial frequency, distance and GC explanations, yet summary associations remain observational and LD-confounded. SAIGE count-model coefficients are not calibrated to GI log(TPM+1) differences. The shorter input excludes distal mechanisms; neither rare variants nor indels were tested. Sequence-context sensitivity, the weak HBB/K562 qualitative control and the unblinded matching adjustment warrant replication under a frozen protocol.',
        'GI training includes GTEx, ENCODE and cellxgene; TenK10K donor or training overlap cannot be excluded. Zero-shot use therefore does not prove an independent training holdout. The runbook’s proximal CETP rsID label also needs reconciliation; coordinates and alleles defined this pilot. A defensible next study should freeze model weights and input geometry, pre-match before scoring, quantify abstention, use the same pairs for PromoterAI/GI comparisons and validate selected effects experimentally.'
    ])
    section('Availability and references',[
        'The run used 444 prediction requests; a complete CLI resume added none. The offline suite passed 49 tests. Code, sequence-free predictions, exclusions, hashes and the full report are in github.com/burtsev/GxI-TenK10K (benchmark results commit b7ecbed). This PDF introduces no new model inference.',
    ])
    refs=[
        '[1] Cuomo ASE et al. Impact of rare and common genetic variation on cell type-specific gene expression in human blood. medRxiv v2. <link href="https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2" color="#233b58">doi:10.1101/2025.03.20.25324352</link>.',
        '[2] Cuomo ASE et al. Rare and common genetic variants shape single-cell immune gene expression. Nature (28 October 2026). <link href="https://www.nature.com/articles/s41586-026-10721-2" color="#233b58">doi:10.1038/s41586-026-10721-2</link>. Published version of [1].',
        '[3] Cuomo ASE. TenK10K Phase 1 association summary statistics. Zenodo, record 18221260. <link href="https://doi.org/10.5281/zenodo.18221260" color="#233b58">doi:10.5281/zenodo.18221260</link>. CC-BY-4.0.'
    ]
    story.extend(p(t,small) for t in refs)
    story.extend([NextPageTemplate('figure'),PageBreak(),p('Figure 1 | Measured pilot performance',head)])
    fig=DATA/'figures/benchmark.png';story.append(Image(str(fig),width=W-2*M,height=(W-2*M)*.7))
    story.append(p('Candidate ALT-oriented effects (upper left); per-cell direction agreement with locus-bootstrap 95% intervals (upper right); measured candidate/control delta distributions (lower left); and REF versus ALT prediction scale (lower right). All panels use the amended 40,960-bp geometry. Zero-effect predictions remain in distributions and correlations but abstain from direction agreement. These are our measured data; no published figure is reproduced.',small))
    story.append(Spacer(1,10))
    table=[['Endpoint','All / sign-eligible','Estimate (95% locus CI)'],['Candidates: direction','100 / 77','67.5% (56.4–77.9%)'],['Controls: direction','100 / 60','51.7% (39.3–64.4%)'],['Cell contrasts: direction','81 / 71','45.1% (27.8–59.0%)'],['Candidate magnitude rank','100 / all','Spearman rho = 0.370'],['Nature: PromoterAI [2]','Different lead-eQTL sets','Median correlation = 0.66']]
    t=Table([[p(escape(v),small) for v in row] for row in table],colWidths=[(W-2*M)*.39,(W-2*M)*.25,(W-2*M)*.36],hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eef3')),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#718298')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),3)]));story.append(t)
    story.append(Spacer(1,8));story.append(p('Published and pilot correlations are not a matched benchmark. The source study’s reported PromoterAI coefficient and our pooled Spearman statistic differ in sampling and aggregation; sign accuracy must not be interpreted as AUROC. All subgroup results and protocol adjustments are exploratory.',small))
    doc.build(story)
    print(pdf)

if __name__=='__main__':main()
