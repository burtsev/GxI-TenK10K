"""Sequential inference, reference reuse via cache, strict effective coverage."""
from ..genome.reference import prepare
from ..genome.variants import covered
from ..util import digest
from ..gi.client import cache_key


def score_pair(client, fasta, row, flank, context):
    ref, alt, tss, index = prepare(fasta, row, flank)
    r = client.predict(ref, context, tss)
    a = client.predict(alt, context, tss)
    if r['window'] != a['window']:
        raise ValueError('REF/ALT effective windows differ')
    if not covered(index, r['window'], len(ref)):
        raise ValueError('Variant outside effective scored window')
    # Legacy model reports a slice, not actual token coverage; primary adapter
    # only permits models whose trained_window_bp is null.
    return dict(row, delta=a['score']-r['score'], ref_score=r['score'], alt_score=a['score'],
                context=context, model=client.model, tss_index=tss, variant_index=index,
                ref_sequence_hash=digest(ref), alt_sequence_hash=digest(alt),
                ref_cache_key=cache_key(ref, context, client.model, tss, client.contract['info']['version'], digest(client.contract), client.base),
                alt_cache_key=cache_key(alt, context, client.model, tss, client.contract['info']['version'], digest(client.contract), client.base),
                scored_start=r['window'][0], scored_end=r['window'][1], status='scored')
