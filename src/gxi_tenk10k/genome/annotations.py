"""GENCODE GTF transcript TSSs; refuse multiple distinct canonical TSSs."""
import gzip
import re


def load_tss(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    genes = {}
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            c = line.rstrip().split('\t')
            if len(c) != 9 or c[2] != 'transcript':
                continue
            attrs = re.findall(r'(\w+) "([^"]+)"', c[8])
            d = dict(attrs)
            if not any(k == 'tag' and v in ('MANE_Select', 'Ensembl_canonical') for k, v in attrs):
                continue
            gene = d['gene_id'].split('.')[0]
            genes.setdefault(gene, set()).add((c[0], int(c[3] if c[6] == '+' else c[4]), c[6]))
    return {g: next(iter(v)) for g, v in genes.items() if len(v) == 1}
