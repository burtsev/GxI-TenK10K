"""Indexed FASTA retrieval. Caller must provide verified GRCh38 provenance."""
import subprocess
import re
from .variants import oriented_pair


def window(tss1, strand, flank):
    # Minus orientation reverses base indices: shift plus window one base to
    # keep the oriented TSS exactly at flank for either gene strand.
    if strand not in ('+', '-') or flank < 4599:
        raise ValueError('Invalid strand/flank')
    start = tss1 - 1 - flank + (strand == '-')
    if start < 0:
        raise ValueError('Window crosses contig boundary')
    return start, start + 2 * flank


def fasta_sequence(path, contig, start0, end0):
    if not re.fullmatch(r'(chr)?([0-9]+|X|Y|MT|M)', contig):
        raise ValueError('Invalid contig')
    result = subprocess.run(['samtools', 'faidx', str(path), f'{contig}:{start0+1}-{end0}'], capture_output=True, text=True, check=True)
    seq = ''.join(result.stdout.splitlines()[1:]).upper()
    if len(seq) != end0-start0 or set(seq)-set('ACGTN'):
        raise ValueError('Reference retrieval failed or contig boundary')
    return seq


def prepare(path, row, flank):
    start, end = window(int(row['tss1']), row['strand'], flank)
    sequence = fasta_sequence(path, row['chrom'], start, end)
    return oriented_pair(sequence, start, int(row['pos1']), row['ref'], row['alt'], int(row['tss1']), row['strand'])
