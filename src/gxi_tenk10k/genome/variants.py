"""All edits occur in genomic plus coordinates before strand transformation."""

def reverse_complement(seq):
    if set(seq) - set('ACGTN'):
        raise ValueError('Invalid DNA')
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def apply_snv(sequence, start0, pos1, ref, alt):
    if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
        raise ValueError('Expected distinct biallelic SNV')
    index = pos1 - 1 - start0
    if not 0 <= index < len(sequence):
        raise ValueError('Variant outside input')
    if sequence[index] != ref:
        raise ValueError('REF mismatch')
    return sequence[:index] + alt + sequence[index + 1:]


def oriented_pair(sequence, start0, pos1, ref, alt, tss1, strand):
    alt_seq = apply_snv(sequence, start0, pos1, ref, alt)
    tss_index = tss1 - 1 - start0
    variant_index = pos1 - 1 - start0
    if not 0 <= tss_index < len(sequence) or strand not in ('+', '-'):
        raise ValueError('Invalid TSS/strand')
    if strand == '-':
        return reverse_complement(sequence), reverse_complement(alt_seq), len(sequence)-1-tss_index, len(sequence)-1-variant_index
    return sequence, alt_seq, tss_index, variant_index


def covered(index, window, length):
    return (isinstance(window, list) and len(window) == 2 and
            all(isinstance(x, int) for x in window) and 0 <= window[0] <= index < window[1] <= length)
