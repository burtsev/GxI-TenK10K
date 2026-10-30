"""Test combined wording, paraphrases and exact repetition; fresh second holdout."""
from pathlib import Path
import pandas as pd
import optimize_context as experiment


def main():
    old=experiment.TEMPLATES.copy()
    canonical='Homo sapiens peripheral-blood {cell}, unstimulated. Single-cell RNA sequencing with 10x Genomics Chromium 3-prime poly(A) capture; source material: peripheral blood mononuclear cells (PBMCs).'
    variants={
        'paraphrase_assay_first':'10x Genomics Chromium 3-prime single-cell RNA sequencing of unstimulated {cell} from human peripheral blood mononuclear cells (PBMCs). The library captures polyadenylated RNA.',
        'paraphrase_cell_first':'Cell identity: {cell}. Organism: Homo sapiens. Tissue: peripheral blood, mononuclear-cell fraction. State: unstimulated. Measurement: 10x Genomics Chromium single-cell RNA-seq with 3-prime poly(A) capture.',
        'paraphrase_sentence':canonical}
    combined=' '.join(old.values())
    experiment.TEMPLATES={'baseline':old['baseline'],'compact':old['compact'],'combined_once':combined,'combined_twice':combined+' '+combined,'compact_twice':old['compact']+' '+old['compact']}
    for name,text in variants.items():
        experiment.TEMPLATES[name+'_once']=text
        experiment.TEMPLATES[name+'_twice']=text+' '+text
    experiment.OUT=Path('results/context-repetition-20261030')
    fresh=pd.read_csv('results/context-optimization-20261030/fresh_panel.tsv',sep='\t')
    pool=pd.read_parquet(experiment.BASE/'association_pool.parquet')
    prior=fresh.merge(pool,on=['gene','SNP','cell','tss1'])[['gene','chrom','tss1']].drop_duplicates()
    experiment.main(extra_prior=prior,incumbent='compact',cap=1600,planned=536,
        preflight={'quota':'User identifies as API developer and confirms quota lifted for this key.','authorization':'User requested combined descriptions, paraphrasing and exact repetition.','rate':'Sequential, rate headers logged, no automatic retries.','cost':'Price not published; no monetary estimate claimed. Existing configured 1600-request ceiling retained.'})

if __name__=='__main__':main()
