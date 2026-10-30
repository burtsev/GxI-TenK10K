"""Real source audit with bounded ZIP streaming, complete counts and manifests."""
import csv
import json
import os
import platform
import resource
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
import httpx
import pandas as pd
from .data.zenodo import metadata, download, members
from .data.parse import read_table, schema
from .pipeline.select import CELL_MAP
from .util import save, file_hash, digest


def audit(config, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    raw = Path(config['paths']['raw'])
    z = metadata()
    save(out/'zenodo.json', z)
    manifest = {'schema_version': 1, 'doi': z['metadata']['doi'], 'version': z['metadata'].get('relations'),
                'license': z['metadata']['license'], 'config_hash': digest(config),
                'files': [], 'total_compressed_bytes': sum(f['size'] for f in z['files'])}
    profiles = []
    for entry in z['files']:
        item = {'name': entry['key'], 'size': entry['size'], 'provided_checksum': entry['checksum'], 'url': entry['links']['self'], 'downloaded': False}
        if entry['key'] in config['sources']['initial_archives']:
            p = download(entry, raw)
            item.update(downloaded=True, sha256=file_hash(p), checksum_verified=True)
            with zipfile.ZipFile(p) as a:
                for info in members(a):
                    with a.open(info) as b:
                        reader, delimiter = read_table(b)
                        kind = schema(reader.fieldnames)
                        n, missing, pip90, pip50 = 0, {c: 0 for c in reader.fieldnames}, 0, 0
                        first, examples = None, []
                        for r in reader:
                            if None in r or any(v is None for v in r.values()):
                                raise ValueError('Malformed row: '+info.filename)
                            n += 1
                            if first is None:
                                first = r
                            for c, v in r.items():
                                missing[c] += str(v).strip().lower() in ('', 'na', 'nan', 'none', 'null')
                            if kind == 'susie':
                                pip = float(r['PIP'])
                                if not 0 <= pip <= 1:
                                    raise ValueError('Invalid source PIP')
                                pip90 += pip >= .9; pip50 += pip >= .5
                            marker = r.get('SNP', r.get('top_MarkerID', ''))
                            if marker.startswith(('1:36478697:', '16:56971389:', '16:57020104:')):
                                examples.append(r)
                    profiles.append(dict(archive=p.name, member=info.filename, uncompressed_bytes=info.file_size,
                                         delimiter=delimiter, columns=reader.fieldnames, kind=kind, rows=n,
                                         missing=missing, pip_ge_90=pip90, pip_ge_50=pip50, first_row=first,
                                         smoke_locus_rows=examples))
        manifest['files'].append(item)
    save(out/'manifest.json', manifest)
    save(out/'archive_profiles.json', profiles)
    # Endpoint schemas are discovered live, not inferred from model names.
    discovery = {}
    for name, url in [('gi_openapi', 'https://api.genomicintelligence.ai/v1/openapi.json'),
                      ('gi_models', 'https://api.genomicintelligence.ai/v1/tasks/expression/models'),
                      ('gi_docs', 'https://docs.genomicintelligence.ai/rest-api.md')]:
        try:
            r = httpx.get(url, timeout=60, follow_redirects=True)
            r.raise_for_status()
            (out/(name+('.md' if name == 'gi_docs' else '.json'))).write_bytes(r.content)
            discovery[name] = {'url': url, 'status': r.status_code, 'sha256': file_hash(out/(name+('.md' if name == 'gi_docs' else '.json')))}
        except (httpx.HTTPError, ValueError) as e:
            discovery[name] = {'url': url, 'error': type(e).__name__}
    save(out/'discovery.json', discovery)
    import importlib.metadata
    env = {'timestamp_utc': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
           'platform': platform.platform(), 'cpus': os.cpu_count(),
           'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
           'config_sha256': digest(config), 'GI_API_KEY_present': bool(os.getenv('GI_API_KEY')),
           'ulimit_nofile': list(resource.getrlimit(resource.RLIMIT_NOFILE)),
           'disk': subprocess.check_output(['df', '-h', '.'], text=True),
           'memory': subprocess.check_output(['free', '-h'], text=True),
           'packages': {p.metadata['Name']: p.version for p in importlib.metadata.distributions()}}
    env['implementation_sha256'] = {str(p): file_hash(p) for root in ['src', 'tests', 'scripts', 'configs'] for p in sorted(Path(root).rglob('*')) if p.is_file() and p.suffix in ('.py', '.yaml')}
    env['implementation_tree_hash'] = digest(env['implementation_sha256'])
    env['git_dirty_at_execution'] = bool(subprocess.check_output(['git', 'status', '--porcelain'], text=True).strip())
    save(out/'environment.json', env)
    counts = {k: sum(p['rows'] for p in profiles if p['kind'] == k) for k in ['gene', 'susie']}
    cells = sorted({p['first_row']['celltype'] for p in profiles if p['kind'] == 'susie'})
    text = f'''# Zenodo data audit

Record: {manifest['doi']}; licence {manifest['license']['id']}. Version relation is saved in manifest.json.
Inventory: {len(z['files'])} files, {manifest['total_compressed_bytes']:,} compressed bytes.
Downloaded only the two prespecified small archives; both provided checksums verified.
Full SHA256 values, URLs and file sizes are recorded in manifest.json.

Counts: {counts}. All {len(profiles)} table members were streamed through to EOF.
Member names, columns, missing-value denominators, first rows and PIP counts are in archive_profiles.json.
Delimiter is detected from contents; .tsv members use commas.

Gene-level schema: gene, ACAT_p, top_MarkerID, top_pval. No beta, SE, effect allele, frequency or donor count.
SuSiE schema: SNP, Index, PIP, Credible_Set, gene, celltype. No beta, SE or frequency.
Variant IDs parse as chromosome:one-based-position:allele:allele; REF/ALT meanings still need reference validation.
No oriented observed betas can be constructed from these two archives.
Missing and convergence-failure rows are counted rather than silently discarded.

Observed cells: {', '.join(cells)}.
Explicit proposed cell mapping: {CELL_MAP}; labels validated during selection.
SuSiE credible sets are 95% sets per Zenodo metadata, not certain causal labels.
No donor-level data are requested or inferred. Sample counts and normalization are not recoverable from these schemas.

Official project source inspection finds GENCODE v44 basic annotation paths in
Scanpy/combine_files_add_gene_info.py; this supports release choice but does not establish one unique canonical TSS per gene.
GRCh38, exact reference alleles, beta allele convention and assay normalization require further verification against full summaries and project code before efficacy scoring.
Full 16.6-GB summary archive deferred because inference credentials are missing; download_large_files=false.
The preprint linked by the record is https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v2 .
Its sample counts have not been independently verified here.
'''
    (out/'data_audit.md').write_text(text)
    return manifest
