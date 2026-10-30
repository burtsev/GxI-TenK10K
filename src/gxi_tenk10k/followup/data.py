"""Deterministic source filtering and bounded full-summary extraction."""
import gzip,json,re,zipfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.stats import false_discovery_control
from ..data.parse import read_table
from ..genome.annotations import load_tss
from ..util import save,file_hash,digest

CELLS={
 'ASDC':('AXL-positive SIGLEC6-positive dendritic cell','dendritic'),
 'B_intermediate':('intermediate B cell','B'),'B_memory':('memory B cell','B'),'B_naive':('B naive','B'),
 'CD14_Mono':('CD14 classical monocyte','monocyte'),'CD16_Mono':('CD16 nonclassical monocyte','monocyte'),
 'CD4_CTL':('CD4 cytotoxic T cell','CD4'),'CD4_Naive':('CD4 naive T','CD4'),'CD4_Proliferating':('proliferating CD4 T cell','CD4'),'CD4_TCM':('central memory CD4 T cell','CD4'),'CD4_TEM':('effector memory CD4 T cell','CD4'),
 'CD8_Naive':('naive CD8 T cell','CD8'),'CD8_Proliferating':('proliferating CD8 T cell','CD8'),'CD8_TCM':('central memory CD8 T cell','CD8'),'CD8_TEM':('effector memory CD8 T cell','CD8'),
 'HSPC':('hematopoietic stem and progenitor cell','progenitor'),'ILC':('innate lymphoid cell','innate'),
 'MAIT':('mucosal-associated invariant T cell','T'),'NK':('natural killer cell','NK'),'NK_CD56bright':('CD56-bright natural killer cell','NK'),'NK_Proliferating':('proliferating natural killer cell','NK'),
 'Plasmablast':('plasmablast','B'),'Treg':('regulatory T cell','CD4'),'cDC1':('conventional dendritic cell type 1','dendritic'),'cDC2':('conventional dendritic cell type 2','dendritic'),'dnT':('double-negative T cell','T'),'gdT':('gamma-delta T cell','T'),'pDC':('plasmacytoid dendritic cell','dendritic')}
# Require exact names in the downloaded official CL release; absence = ineligible.
CL_NAMES={'B_naive':'naive B cell','B_memory':'memory B cell','CD14_Mono':'classical monocyte','CD16_Mono':'non-classical monocyte','CD4_Naive':'naive thymus-derived CD4-positive, alpha-beta T cell','CD8_Naive':'naive thymus-derived CD8-positive, alpha-beta T cell','CD4_TCM':'central memory CD4-positive, alpha-beta T cell','CD8_TCM':'central memory CD8-positive, alpha-beta T cell','CD4_TEM':'effector memory CD4-positive, alpha-beta T cell','CD8_TEM':'effector memory CD8-positive, alpha-beta T cell','NK':'natural killer cell','NK_CD56bright':'CD16-negative, CD56-bright natural killer cell, human','Plasmablast':'plasmablast','Treg':'regulatory T cell','gdT':'gamma-delta T cell','pDC':'plasmacytoid dendritic cell','cDC1':'CD141-positive myeloid dendritic cell','cDC2':'CD1c-positive myeloid dendritic cell','ILC':'innate lymphoid cell'}


def annotations():
    p=Path('data/raw/gencode44/gencode.v44.basic.annotation.gtf.gz');tss=load_tss(p);names={}
    with gzip.open(p,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            c=line.rstrip().split('\t')
            if c[2]!='gene':continue
            a=dict(re.findall(r'(\w+) "([^"]+)"',c[8]));names[a['gene_id'].split('.')[0]]=a['gene_name']
    return pd.DataFrame([dict(gene=g,chrom=v[0],tss1=v[1],strand=v[2],symbol=names.get(g,g)) for g,v in tss.items() if re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)',v[0])])


def cell_mapping(out):
    p=Path('data/raw/cell-ontology.obo');text=p.read_text();lookup={}
    for b in text.split('[Term]'):
        lines=b.splitlines();identifier=next((x[4:] for x in lines if x.startswith('id: ')),'');name=next((x[6:] for x in lines if x.startswith('name: ')),'')
        if identifier.startswith('CL:') and 'is_obsolete: true' not in lines:lookup.setdefault(name,[]).append(identifier)
    rows=[]
    for cell,(canonical,family) in CELLS.items():
        name=CL_NAMES.get(cell,'');ids=lookup.get(name,[])
        rows.append(dict(cell=cell,canonical=canonical,family=family,ontology_label=name if len(ids)==1 else '',cl_id=ids[0] if len(ids)==1 else '',ontology_eligible=len(ids)==1,ontology_source='https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl.obo',ontology_sha256=file_hash(p)))
    f=pd.DataFrame(rows);f.to_csv(out/'cell_mapping.tsv',sep='\t',index=False,na_rep='NA');return f


def cluster_genes(ann,legacy):
    x=pd.concat([ann[['gene','chrom','tss1']],legacy[['gene','chrom','tss1']]]).drop_duplicates(['gene','chrom','tss1']).sort_values(['chrom','tss1','gene'])
    groups={};previous=None;end=-1;index=0
    for r in x.itertuples(index=False):
        if r.chrom!=previous or r.tss1-100000>end:index+=1;end=r.tss1+100000
        else:end=max(end,r.tss1+100000)
        previous=r.chrom;groups[r.gene]=f'cluster_{index:05d}'
    old={groups[g] for g in legacy.gene if g in groups}
    return groups,old


def prescreen(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    ann=annotations();significant=set();raw=Path('data/raw')
    with zipfile.ZipFile(raw/'common_variant_gene_level_results_with_annotated_fails.zip') as a:
        for name in a.namelist():
            if not name.endswith('.tsv'):continue
            cell=Path(name).name.split('_all_cis')[0]
            if cell not in CELLS:raise ValueError('Unknown cell '+cell)
            with a.open(name) as b:r,_=read_table(b);f=pd.DataFrame(r)
            p=pd.to_numeric(f.ACAT_p,errors='coerce');ok=p.between(0,1);q=false_discovery_control(p[ok])
            significant.update((g,cell) for g in f.loc[ok].loc[q<=.05,'gene'])
    frames=[]
    with zipfile.ZipFile(raw/'susie_summary.zip') as a:
        for name in a.namelist():
            if not name.endswith('.tsv'):continue
            with a.open(name) as b:r,_=read_table(b);f=pd.DataFrame(r)
            f['PIP']=pd.to_numeric(f.PIP);f=f[f.PIP>=.9]
            f=f[[(g,c) in significant for g,c in zip(f.gene,f.celltype)]]
            frames.append(f)
    high=pd.concat(frames,ignore_index=True).rename(columns={'celltype':'discovery_cell'})
    high=high[high.SNP.str.fullmatch(r'(?:[0-9]+|X|Y):[1-9][0-9]*:[ACGT]:[ACGT]')].copy()
    parts=high.SNP.str.split(':',expand=True);high['pos1']=parts[1].astype(int);high['ref']=parts[2];high['alt']=parts[3];high['variant_chrom']='chr'+parts[0]
    high=high[high.ref!=high.alt].merge(ann,on='gene',validate='many_to_one')
    high=high[(high.chrom==high.variant_chrom)&(abs(high.pos1-high.tss1)<20380)&(high.tss1>20481)]
    legacy=pd.read_csv(out/'recovery/legacy_exclusions.tsv',sep='\t')
    groups,old=cluster_genes(high[['gene','chrom','tss1']].drop_duplicates(),legacy)
    high['locus']=high.gene.map(groups);before=len(high);high=high[~high.locus.isin(old)]
    high['order']=[digest([20261030,r.locus,r.gene,r.SNP,r.discovery_cell]) for r in high.itertuples()]
    high=high.sort_values('order').drop_duplicates(['gene','SNP'])
    # Choose fixed disjoint locus pools before viewing any predictions. Reserve
    # surplus candidate loci for input/matching attrition; never replace after scoring.
    clusters=high[['locus']].drop_duplicates().copy();clusters['order']=[digest([20261030,l]) for l in clusters.locus];clusters=clusters.sort_values('order')
    dev=set(clusters.head(450).locus);high['split']=np.where(high.locus.isin(dev),'development','test')
    high.to_parquet(out/'candidate_reserve.parquet',index=False)
    cell_mapping(out)
    save(out/'prescreen.json',{'high_pip_gene_cell_rows_before_historical_exclusion':before,'reserved_pairs':len(high),'reserved_genes':high.gene.nunique(),'reserved_loci':high.locus.nunique(),'historical_cluster_count':len(old),'development_reserve_clusters':len(dev),'gene_fdr_method':'BH q<=0.05 on source ACAT per cell; prespecified conservative gene eligibility','primary_pip':.9,'source_cell_count':len(CELLS)})
    print('Prescreen',len(high),'pairs',high.gene.nunique(),'genes',high.locus.nunique(),'clusters',flush=True)
    return high


def extract(out):
    out=Path(out);reserve=pd.read_parquet(out/'candidate_reserve.parquet');ann=annotations().set_index('gene');genes=set(reserve.gene)
    folder=out/'summary_pool';folder.mkdir(exist_ok=True)
    usecols=['CHR','POS','MarkerID','Allele1','Allele2','AF_Allele2','BETA','SE','N','gene']
    def one(cell):
        target=folder/(cell+'.parquet');done=folder/(cell+'.json')
        if target.exists() and done.exists():return json.load(open(done))
        count=0;kept=0;writer=None;tmp=target.with_suffix('.partial')
        with zipfile.ZipFile('data/raw/common_all_cis_pvalues_100kb.zip') as a:
            name=cell+'_common_all_cis_raw_pvalues.tsv';info=a.getinfo(name)
            if info.file_size>4_000_000_000 or info.file_size/max(1,info.compress_size)>100:raise ValueError('Unsafe source expansion')
            with a.open(info) as b:
                for chunk in pd.read_csv(b,sep='\t',usecols=usecols,chunksize=200000,dtype={'CHR':str,'gene':str,'MarkerID':str,'Allele1':str,'Allele2':str}):
                    count+=len(chunk);f=chunk[chunk.gene.isin(genes)].copy()
                    if f.empty:continue
                    f['tss1']=f.gene.map(ann.tss1);f=f[abs(f.POS-f.tss1)<20380]
                    if f.empty:continue
                    f=f.rename(columns={'POS':'pos1','MarkerID':'SNP','Allele1':'ref','Allele2':'alt','BETA':'beta','AF_Allele2':'af'})
                    f['cell']=cell;f['chrom']='chr'+f.CHR;f['strand']=f.gene.map(ann.strand);f['symbol']=f.gene.map(ann.symbol)
                    f['maf']=np.minimum(f.af,1-f.af);f['distance']=abs(f.pos1-f.tss1)
                    f['beta_alt']=f.beta;f['effect_allele']=f.alt;f['genome_build']='GRCh38';f['source_file']='common_all_cis_pvalues_100kb.zip:'+name
                    table=pa.Table.from_pandas(f,preserve_index=False)
                    if writer is None:writer=pq.ParquetWriter(tmp,table.schema,compression='zstd')
                    writer.write_table(table);kept+=len(f)
        if writer is None:raise ValueError('Empty source cell '+cell)
        writer.close();tmp.replace(target)
        item={'cell':cell,'source_rows':count,'retained_rows':kept,'source_crc32':f'{info.CRC:08x}','sha256':file_hash(target)};save(done,item)
        print('Extracted',cell,count,kept,flush=True);return item
    with ThreadPoolExecutor(max_workers=2) as ex:profiles=list(ex.map(one,sorted(CELLS)))
    save(out/'full_source_audit.json',profiles)
