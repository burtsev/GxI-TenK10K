"""Recheck cached official sources and audit retained-data missingness."""
import json,zipfile
from pathlib import Path
import pandas as pd,numpy as np
from ..util import file_hash,save

def verify(out):
    out=Path(out);m=json.load(open(out/'source_manifest.json'));live=json.load(open(out/'zenodo_live.json'));inventory={x['key']:x for x in live['files']};checks=[]
    for row in m['files']:
        path=Path('data/raw')/row['name'];entry=inventory[row['name']]
        assert entry['checksum']==row['provided_checksum'] and entry['size']==path.stat().st_size
        assert file_hash(path)==row['sha256'];alg,expected=entry['checksum'].split(':');assert file_hash(path,alg)==expected
        checks.append(dict(name=row['name'],sha256=row['sha256'],provided_checksum_verified=True,live_metadata_verified=True))
    root=Path('data/raw/gencode44');expected={l.split()[1].lstrip('*'):l.split()[0] for l in (root/'MD5SUMS').read_text().splitlines() if len(l.split())==2}
    for row in m['reference']:
        p=root/row['name'];assert file_hash(p)==row['sha256'];assert file_hash(p,'md5')==expected[row['name']]==row['md5']
    prior=json.load(open('docs/runs/live-20261029/reference_manifest.json'));fa=root/'GRCh38.primary_assembly.genome.fa';assert file_hash(fa)==prior['sha256']
    missing=[]
    for path in sorted((out/'summary_pool').glob('*.parquet')):
        assert file_hash(path)==json.load(open(path.with_suffix('.json')))['sha256']
        f=pd.read_parquet(path);counts={c:int(f[c].isna().sum()) for c in f.columns};missing.append(dict(cell=path.stem,n_retained=len(f),missing_by_field=counts,invalid_beta=int((~np.isfinite(f.beta)).sum()),invalid_se=int(((~np.isfinite(f.SE))|(f.SE<=0)).sum()),invalid_maf=int((~f.maf.between(0,.5)).sum()),min_N=float(f.N.min()),max_N=float(f.N.max())))
    schemas=[]
    with zipfile.ZipFile('data/raw/common_all_cis_pvalues_100kb.zip') as archive:
        from ..data.parse import read_table
        for name in archive.namelist():
            if not name.endswith('.tsv'):continue
            with archive.open(name) as raw:
                reader,delimiter=read_table(raw);first=next(reader);schemas.append({'member':name,'delimiter':delimiter,'columns':reader.fieldnames,'sample_fields':len(first)})
    save(out/'data_verification.json',{'archive_checks':checks,'reference':prior,'retained_data_missingness':missing,'prior_full_schema_audit':schemas,'all_full_source_members_read_to_eof_and_crc_checked':len(json.load(open(out/'full_source_audit.json')))==28,'missingness_scope':'retained source pool only; entire-archive missing-field counts not computed','live_doi':live['metadata']['doi']})
    return checks
