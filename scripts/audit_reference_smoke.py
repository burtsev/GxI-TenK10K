"""Independent real reference checks; not an authenticated GI smoke test.

Fetches official Ensembl GRCh38 gene records and plus-strand bases only.
Runtime reference snapshots remain in ignored data/; no sequences in Git.
"""
import json
from pathlib import Path
import httpx
from gxi_tenk10k.util import file_hash, save


def main():
    out=Path('results/pilot-20261029')
    raw=Path('data/raw/ensembl');raw.mkdir(parents=True,exist_ok=True)
    checks=[]
    with httpx.Client(base_url='https://rest.ensembl.org',timeout=60,headers={'Content-Type':'application/json'}) as c:
        for gene,chrom,pos,ref,alt in [('ENSG00000119535','1',36478697,'C','A'),('ENSG00000087237','16',56971389,'C','T'),('ENSG00000087237','16',57020104,'G','T')]:
            try:
                r=c.get('/lookup/id/'+gene,params={'content-type':'application/json'});r.raise_for_status()
                annotation=r.json();save(raw/(gene+'.json'),annotation)
                if annotation['assembly_name']!='GRCh38':raise ValueError('Unexpected assembly')
                r=c.get(f'/sequence/region/human/{chrom}:{pos}..{pos}:1',params={'content-type':'application/json'});r.raise_for_status()
                reference=r.json();save(raw/(f'{chrom}-{pos}.json'),reference)
                tss=annotation['start'] if annotation['strand']==1 else annotation['end']
                checks.append(dict(gene=gene,variant=f'{chrom}:{pos}:{ref}:{alt}',assembly=annotation['assembly_name'],
                    genomic_ref_matches=reference['seq'].upper()==ref,ensembl_gene_boundary_tss=tss,
                    strand=annotation['strand'],distance_from_gene_boundary_tss=pos-tss,
                    annotation_sha256=file_hash(raw/(gene+'.json')),sequence_record_sha256=file_hash(raw/(f'{chrom}-{pos}.json')),
                    note='Ensembl current gene boundary only; not a validated GENCODE v44 canonical transcript TSS'))
            except (httpx.HTTPError,ValueError,KeyError) as e:
                checks.append(dict(gene=gene,variant=f'{chrom}:{pos}:{ref}:{alt}',error=type(e).__name__))
    save(out/'smoke_reference_audit.json',checks)
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
