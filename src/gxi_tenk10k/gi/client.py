"""Adapter for the archived live OpenAPI; bounded retries and disk checkpoints."""
import json
import math
import os
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
import httpx
from ..util import digest, save
from ..genome.variants import covered


class Blocked(RuntimeError):
    pass


def cache_key(sequence, context, model, tss, api_version, contract_hash, base_url):
    return digest(dict(sequence_hash=digest(sequence), context=context, model=model,
                       tss_index=tss, api_version=api_version, contract_hash=contract_hash,
                       base_url=base_url, options={'description': context}))


def validate_response(data, sequence, context, model, tss):
    payload, meta = data['data'], data['meta']
    echo = payload['input']
    if payload['task'] != 'expression' or payload['model'] != model:
        raise ValueError('Task/model mismatch')
    if echo['description'] != context or echo['tss_index'] != tss or meta['sequence_length'] != len(sequence):
        raise ValueError('Context/TSS/length mismatch')
    prediction = payload['prediction']
    if prediction['unit'].replace(' ', '') != 'log(TPM+1)':
        raise ValueError('Unverified response unit')
    score = float(prediction['expression_log_tpm'])
    if not math.isfinite(score):
        raise ValueError('Nonfinite prediction')
    counts = meta['task_specific_counts']
    if counts['tss_index'] != tss:
        raise ValueError('Scored TSS mismatch')
    scored = counts.get('scored_window')
    if not covered(tss, scored, len(sequence)):
        raise ValueError('Unknown/invalid scored window')
    return {'score': score, 'window': scored, 'model': model, 'unit': prediction['unit']}


class Client:
    def __init__(self, cache, contract, model, max_requests=1600, timeout=120,
                 retries=5, base_url=None, transport=None, sleep=time.sleep):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.contract = contract
        self.model = model
        self.max_requests = max_requests
        self.calls = 0
        self.retries = retries
        self.sleep = sleep
        self.base = (base_url or os.getenv('GI_BASE_URL', 'https://api.genomicintelligence.ai')).rstrip('/')
        if not self.base.startswith('https://'):
            raise ValueError('HTTPS API required')
        self.key = os.getenv('GI_API_KEY')
        self.http = httpx.Client(base_url=self.base, timeout=timeout, transport=transport)
        self.ledger = self.cache / 'requests.jsonl'
        # Budget includes attempts from interrupted runs, including failed calls.
        self.spent = sum(1 for _ in self.ledger.open()) if self.ledger.exists() else 0

    def predict(self, sequence, context, tss):
        if not 9198 <= len(sequence) <= 500000 or not 4599 <= tss <= len(sequence)-4599:
            raise ValueError('Input geometry violates live contract')
        if set(sequence)-set('ACGTN') or not context.strip():
            raise ValueError('Invalid sequence/context')
        key = cache_key(sequence, context, self.model, tss, self.contract['info']['version'], digest(self.contract), self.base)
        path = self.cache / (key + '.json')
        if path.exists():
            return validate_response(json.loads(path.read_text()), sequence, context, self.model, tss)
        if not self.key:
            raise Blocked('GI_API_KEY absent')
        body = {'sequence': sequence, 'model': self.model, 'tss_index': tss, 'options': {'description': context}}
        for attempt in range(self.retries+1):
            if self.spent >= self.max_requests:
                raise Blocked('Cumulative request budget exhausted')
            # Persist before dispatch so interruption cannot erase expenditure.
            with self.ledger.open('a') as f:
                f.write(json.dumps({'key': key, 'attempt': attempt, 'model': self.model,
                                    'timestamp': datetime.now(timezone.utc).isoformat()})+'\n')
            self.calls += 1
            self.spent += 1
            try:
                r = self.http.post('/v1/tasks/expression/predict', json=body,
                                   headers={'Authorization': 'Bearer '+self.key})
            except httpx.TransportError:
                # Ambiguous billing: no automatic replay of a transport failure.
                raise Blocked('Transport failure; request may have been billed') from None
            with (self.cache/'responses.jsonl').open('a') as f:
                f.write(json.dumps({'key': key, 'status': r.status_code,
                    'rate_headers': {k: v for k, v in r.headers.items() if k.startswith('ratelimit') or k == 'retry-after'}})+'\n')
            if r.status_code in (429, 500, 502, 503, 504) and attempt < self.retries:
                delay = r.headers.get('retry-after', str(min(2**attempt, 30)))
                try:
                    seconds = float(delay)
                except ValueError:
                    seconds = (parsedate_to_datetime(delay)-datetime.now(timezone.utc)).total_seconds()
                if seconds > 60:
                    raise Blocked('Retry-After exceeds bounded wait; resume later')
                self.sleep(max(seconds, 0))
                continue
            if r.status_code != 200:
                raise Blocked(f'GI HTTP {r.status_code}; see redacted response ledger')
            data = r.json()
            result = validate_response(data, sequence, context, self.model, tss)
            save(path, data)
            return result
        raise Blocked('Retries exhausted')
