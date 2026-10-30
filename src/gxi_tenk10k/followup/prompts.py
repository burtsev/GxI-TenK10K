"""Immutable input-only description registry; no association outcomes used."""
import json
from pathlib import Path
import pandas as pd
from ..util import save,digest,file_hash
ARMS=['baseline','canonical_cell_only','cell_plus_verified_cl_id','verified_cl_id_only','structured_verified_assay','baseline_plus_gene_symbol','baseline_plus_ensembl_id','baseline_plus_symbol_and_ensembl','historical_repeated_paraphrase','duplicated_baseline']
GENE_ARMS=set(ARMS[5:8])
def descriptions(row,mapping):
    cell=mapping[row['cell']];name=cell['canonical'];cl=cell.get('cl_id')
    baseline=f'Human peripheral blood {name}, unstimulated, 10x 3-prime single-cell RNA-seq'
    historical=f'Homo sapiens peripheral-blood {name}, unstimulated. Single-cell RNA sequencing with 10x Genomics Chromium 3-prime poly(A) capture; source material: peripheral blood mononuclear cells (PBMCs). '
    eligible=bool(cell['ontology_eligible']) and pd.notna(cl)
    return dict(zip(ARMS,[baseline,name,f'{name}; {cl}' if eligible else None,cl if eligible else None,f'Organism: Homo sapiens. Tissue: peripheral blood mononuclear cells. Cell: {name}. State: unstimulated. Assay: 10x Genomics 3-prime single-cell RNA sequencing.',baseline+f". Gene: {row['symbol']}",baseline+f". Gene: {row['gene']}",baseline+f". Gene: {row['symbol']} ({row['gene']})",historical*2,baseline+' '+baseline]))
def registry(out):
    out=Path(out);mapping=pd.read_csv(out/'cell_mapping.tsv',sep='\t').set_index('cell').to_dict('index')
    # Actual source header is explicitly checked, never silently substituted.
    rows=[]
    for split in ['development','test']:
        f=pd.read_parquet(out/f'cohort_{split}.parquet')
        for r in f[['gene','symbol','cell']].drop_duplicates().to_dict('records'):
            for arm,context in descriptions(r,mapping).items():rows.append(dict(split=split,**r,arm=arm,description=context,context_hash=digest(context) if context else None,gene_aware=arm in GENE_ARMS))
    target=out/'prompt_registry.parquet';frame=pd.DataFrame(rows)
    if target.exists():
        old=pd.read_parquet(target)
        if digest(old.to_json())!=digest(frame.to_json()):raise ValueError('Prompt registry changed')
    else:frame.to_parquet(target,index=False)
    return frame
