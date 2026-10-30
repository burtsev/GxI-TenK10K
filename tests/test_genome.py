import pytest
from gxi_tenk10k.genome.variants import apply_snv, oriented_pair, reverse_complement, covered
from gxi_tenk10k.genome.reference import window, fasta_sequence
from gxi_tenk10k.genome.annotations import load_tss
from gxi_tenk10k.data.parse import variant, orient_beta, read_table, schema
import io


def test_plus_and_minus():
    s = 'AACCGGTT'
    plus = oriented_pair(s, 100, 103, 'C', 'A', 105, '+')
    minus = oriented_pair(s, 100, 103, 'C', 'A', 105, '-')
    assert plus == ('AACCGGTT', 'AAACGGTT', 4, 2)
    assert minus == (reverse_complement(s), reverse_complement(plus[1]), 3, 5)
    assert len(plus[0]) == len(plus[1])


@pytest.mark.parametrize('pos,ref,alt', [(100,'A','C'), (109,'A','C'), (103,'A','C'), (103,'C','C'), (103,'C','CT'), (103,'C','N')])
def test_invalid_edits(pos,ref,alt):
    with pytest.raises(ValueError):
        apply_snv('AACCGGTT',100,pos,ref,alt)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_tss_geometry(strand):
    flank, tss = 4599, 10000
    start,end = window(tss, strand, flank)
    seq = 'A'*(end-start)
    ref,alt,index,vindex = oriented_pair(seq,start,tss,'A','C',tss,strand)
    assert index == flank == vindex
    assert len(ref) == 9198
    assert ref[index] != alt[index]


def test_scored_bounds():
    assert covered(3,[3,7],10)
    assert not covered(7,[3,7],10)
    assert not covered(3,None,10)
    assert not covered(3,[3,11],10)


def test_same_orientation_on_minus():
    assert orient_beta(2,'G','A','G') == 2
    assert orient_beta(2,'A','A','G') == -2
    with pytest.raises(ValueError): orient_beta(2,'T','A','G')
    with pytest.raises(ValueError): orient_beta(float('nan'),'A','A','G')


def test_parse_delimiter_and_unknown_schema():
    r,delimiter = read_table(io.BytesIO(b'gene,ACAT_p,top_MarkerID,top_pval\ng,.1,1:1:A:C,.01\n'))
    assert delimiter == ',' and schema(r.fieldnames) == 'gene'
    assert next(r)['top_MarkerID'] == '1:1:A:C'
    with pytest.raises(ValueError): schema(['guess'])
    assert variant('chr1:123:A:G') == ('1',123,'A','G')
    for malformed in ['1:0:A:G','1:10:A:A','nonsense']:
        with pytest.raises(ValueError): variant(malformed)


def test_fasta_boundaries(tmp_path):
    import subprocess
    path = tmp_path/'ref.fa'; path.write_text('>1\nAACCGGTT\n')
    subprocess.run(['samtools','faidx',str(path)],check=True)
    assert fasta_sequence(path,'1',2,6) == 'CCGG'
    with pytest.raises(ValueError): fasta_sequence(path,'1',0,10)


def test_ambiguous_annotation_rejected(tmp_path):
    path=tmp_path/'genes.gtf'
    path.write_text('1\ts\ttranscript\t100\t200\t.\t+\t.\tgene_id "g.1"; tag "Ensembl_canonical";\n1\ts\ttranscript\t101\t200\t.\t+\t.\tgene_id "g.1"; tag "MANE_Select";\n')
    assert load_tss(path) == {}
