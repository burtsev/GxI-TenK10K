"""Read-only validation of the preregistered input and implementation binding."""
import json
from pathlib import Path
from ..util import file_hash

def validate_selection(out: str | Path) -> dict:
    out=Path(out);lock=json.load(open(out/'context_selection.json'))
    for name,expected in lock['file_hashes'].items():
        if file_hash(out/name)!=expected:raise ValueError('Frozen selection input changed: '+name)
    for name,expected in json.load(open(out/'code_manifest.json')).items():
        if file_hash(name)!=expected:raise ValueError('Frozen implementation changed: '+name)
    return lock

def validate_reference(out: str | Path, fasta: str | Path = 'data/raw/gencode44/GRCh38.primary_assembly.genome.fa') -> None:
    source=json.load(open(Path(out)/'data_verification.json'));fasta=Path(fasta)
    if file_hash(fasta)!=source['reference']['sha256']:raise ValueError('Frozen GRCh38 FASTA changed')
    if file_hash(str(fasta)+'.fai')!=source['reference_index_sha256']:raise ValueError('Frozen FASTA coordinate index changed')
