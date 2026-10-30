"""Normalize the now-audited full SAIGE release with explicit allele semantics."""
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from .select import screen, CELL_MAP
from ..data.parse import read_table, variant
from ..genome.annotations import load_tss
from ..genome.reference import window, fasta_sequence
from ..util import save, file_hash, digest

EFFECT_SOURCE = 'https://weizhou0.github.io/SAIGE-QTL-doc/docs/single_step2.html (Allele2 effect allele; BETA effect of Allele2)'
CELLS = list(CELL_MAP.values())


def prepare(config, out, target=300):
    out=Path(out);raw=Path(config['paths']['raw'])
    annotation=raw/'gencode44/gencode.v44.basic.annotation.gtf.gz'
    fasta=raw/'gencode44/GRCh38.primary_assembly.genome.fa'
    mapping=load_tss(annotation)
    # Start with the previously audited significant genes; PIP/variant screening
    # is expanded by the same seed and rule to provide coverage replacements.
    significant=set()
    import scipy.stats
    with zipfile.ZipFile(raw/'common_variant_gene_level_results_with_annotated_fails.zip') as a:
        for name in a.namelist():
            cell=Path(name).name.split('_all_cis')[0]
            if cell not in CELLS:continue
            with a.open(name) as b:r,_=read_table(b);f=pd.DataFrame(r)
            p=pd.to_numeric(f.ACAT_p,errors='coerce');valid=p.between(0,1)
            q=scipy.stats.false_discovery_control(p[valid]);significant.update((g,cell) for g in f.loc[valid].loc[q<=.05,'gene'])
    frames=[]
    with zipfile.ZipFile(raw/'susie_summary.zip') as a:
        for name in a.namelist():
            if not name.endswith('.tsv'):continue
            cell=Path(name).name.removesuffix('_all_credible_snps.tsv')
            if cell not in CELLS:continue
            with a.open(name) as b:r,_=read_table(b);f=pd.DataFrame(r)
            f['PIP']=pd.to_numeric(f.PIP)
            f['source_file']='susie_summary.zip:'+name
            frames.append(f)
    susie=pd.concat(frames,ignore_index=True)
    high=susie[(susie.PIP>=.9)&susie.gene.isin(mapping)].copy()
    high=high[[ (g,c) in significant for g,c in zip(high.gene,high.celltype) ]]
    keep=[]
    for row in high.to_dict('records'):
        chrom,tss,strand=mapping[row['gene']]
        try:vc,pos,ref,alt=variant(row['SNP'])
        except ValueError:
            keep.append(False);continue
        try:start,end=window(tss,strand,40960)
        except ValueError:
            keep.append(False);continue
        keep.append(chrom.removeprefix('chr')==vc and start<=pos-1<end)
    high=high[keep]
    candidates=screen(high,config['study']['seed'],target,.9)
    # No gene boundary surrogate: all TSSs come from validated GENCODE v44.
    candidates['tss1']=[mapping[g][1] for g in candidates.gene]
    candidates['strand']=[mapping[g][2] for g in candidates.gene]
    candidates['chrom']=[mapping[g][0] for g in candidates.gene]
    candidates['cell']=candidates.celltype
    candidates.to_parquet(out/'expanded_prescreen.parquet',index=False)
    genes=set(candidates.gene)|{'ENSG00000119535','ENSG00000087237'}
    source=raw/'common_all_cis_pvalues_100kb.zip'
    sourcehash=file_hash(source)
    usecols=['CHR','POS','MarkerID','Allele1','Allele2','AF_Allele2','BETA','SE','p.value','N','gene']
    counts=[];parts=[]
    with zipfile.ZipFile(source) as a:
        for cell in CELLS:
            matches=[i for i in a.infolist() if i.filename==f'{cell}_common_all_cis_raw_pvalues.tsv']
            if len(matches)!=1:raise ValueError('Unexpected full-summary inventory')
            info=matches[0]
            if info.file_size>4_000_000_000 or info.compress_size<1 or info.file_size/info.compress_size>100:
                raise ValueError('Unsafe full summary expansion')
            n=0;kept=0
            with a.open(info) as b:
                for chunk in pd.read_csv(b,sep='\t',usecols=usecols,chunksize=200000,dtype={'CHR':str,'gene':str,'MarkerID':str,'Allele1':str,'Allele2':str}):
                    n+=len(chunk);f=chunk[chunk.gene.isin(genes)].copy();kept+=len(f)
                    f['cell']=cell;f['source_file']=source.name+':'+info.filename
                    parts.append(f)
            counts.append({'cell':cell,'rows_scanned':n,'rows_retained':kept,'member_crc32':f'{info.CRC:08x}'})
            print('Scanned',cell,n,'retained',kept,flush=True)
    full=pd.concat(parts,ignore_index=True)
    full=full.rename(columns={'MarkerID':'SNP','POS':'pos1','Allele1':'ref','Allele2':'alt','BETA':'beta','AF_Allele2':'af'})
    full['chrom']='chr'+full.CHR.astype(str)
    full['effect_allele']=full.alt;full['effect_allele_source']=EFFECT_SOURCE
    full['genome_build']='GRCh38';full['maf']=np.minimum(full.af,1-full.af)
    full['tss1']=[mapping[g][1] for g in full.gene];full['strand']=[mapping[g][2] for g in full.gene]
    full['distance']=abs(full.pos1-full.tss1)
    # Merge all credible-set membership; absent PIP remains unknown, never zero.
    pips=susie.groupby(['gene','SNP','celltype'],as_index=False).PIP.max().rename(columns={'celltype':'cell'})
    full=full.merge(pips,on=['gene','SNP','cell'],how='left',validate='many_to_one')
    full.to_parquet(out/'association_pool.parquet',index=False)
    save(out/'full_summary_audit.json',{'source_sha256':sourcehash,'source_md5':file_hash(source,'md5'),
        'effect_allele_source':EFFECT_SOURCE,'scanned_members':counts,'normalization':'SAIGE-QTL count-model beta; not log(TPM+1)',
        'coordinate':'GRCh38 one-based POS, Allele1=REF, Allele2=ALT/effect', 'genes_retained':len(genes)})
    return candidates,full
