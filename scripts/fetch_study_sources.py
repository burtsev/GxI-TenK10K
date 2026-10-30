"""Fetch official study/reference sources; full association ZIP is explicit opt-in."""
import argparse,gzip,shutil
from pathlib import Path
import httpx
from gxi_tenk10k.data.zenodo import metadata,download
from gxi_tenk10k.util import save,file_hash


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--include-full',action='store_true');args=p.parse_args()
    z=metadata();names={'susie_summary.zip','common_variant_gene_level_results_with_annotated_fails.zip'}
    if args.include_full:names.add('common_all_cis_pvalues_100kb.zip')
    manifest=[]
    for entry in z['files']:
        if entry['key'] not in names:continue
        path=download(entry,'data/raw',allow_large=args.include_full)
        manifest.append(dict(name=path.name,url=entry['links']['self'],provided_checksum=entry['checksum'],sha256=file_hash(path),size=path.stat().st_size))
    root=Path('data/raw/gencode44');root.mkdir(parents=True,exist_ok=True)
    url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/'
    for name in ['MD5SUMS','gencode.v44.basic.annotation.gtf.gz','GRCh38.primary_assembly.genome.fa.gz']:
        path=root/name
        if not path.exists():
            part=path.with_suffix(path.suffix+'.part');n=0
            with httpx.stream('GET',url+name,timeout=120,follow_redirects=True) as r:
                r.raise_for_status()
                with part.open('wb') as f:
                    for b in r.iter_bytes(1024*1024):
                        n+=len(b)
                        if n>1_000_000_000:raise ValueError('Reference download exceeds cap')
                        f.write(b)
            part.replace(path)
    expected={l.split()[1].lstrip('*'):l.split()[0] for l in (root/'MD5SUMS').read_text().splitlines() if len(l.split())==2}
    references=[]
    for name in ['gencode.v44.basic.annotation.gtf.gz','GRCh38.primary_assembly.genome.fa.gz']:
        path=root/name;actual=file_hash(path,'md5')
        if actual!=expected[name]:raise ValueError('GENCODE checksum failed '+name)
        references.append(dict(name=name,url=url+name,md5=actual,sha256=file_hash(path)))
    fasta=root/'GRCh38.primary_assembly.genome.fa'
    if not fasta.exists():
        temp=fasta.with_suffix('.fa.part');n=0
        with gzip.open(root/'GRCh38.primary_assembly.genome.fa.gz','rb') as r,temp.open('wb') as f:
            while b:=r.read(1024*1024):
                n+=len(b)
                if n>3_500_000_000:raise ValueError('Reference expansion exceeds cap')
                f.write(b)
        temp.replace(fasta)
    import subprocess
    subprocess.run(['samtools','faidx',str(fasta)],check=True)
    out=Path('results/pilot-20261029');out.mkdir(exist_ok=True,parents=True)
    save(out/'reference_manifest.json',{'assembly':'GRCh38','annotation_release':'GENCODE v44','sha256':file_hash(fasta),'sources':references})
    save(out/'complete_source_manifest.json',{'doi':z['metadata']['doi'],'license':z['metadata']['license'],'files':manifest,'reference':references,'full_archive_opt_in':args.include_full})
    print('Verified',len(manifest),'study sources and GENCODE v44 reference')

if __name__=='__main__':main()
