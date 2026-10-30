"""Auditable stage-gated command line. Synthetic tests never enter live reports."""
import argparse
import json
import os
import sys
import zipfile
from pathlib import Path
import pandas as pd
import yaml
from scipy.stats import false_discovery_control
from .audit import audit
from .data.parse import read_table, orient_beta, variant
from .data.zenodo import members
from .pipeline.select import screen, CELL_MAP
from .pipeline.score import score_pair
from .analysis.metrics import evaluate, split_loci
from .gi.client import Client, Blocked
from .report import report
from .util import save, digest, file_hash

COLUMNS = ['gene', 'locus', 'cell', 'SNP', 'chrom', 'pos1', 'ref', 'alt', 'beta_alt',
           'effect_allele', 'PIP', 'source_file', 'genome_build', 'model', 'status']


def select(config, out):
    frames, genes = [], []
    for filename in config['sources']['initial_archives']:
        with zipfile.ZipFile(Path(config['paths']['raw'])/filename) as a:
            for info in members(a):
                with a.open(info) as b:
                    reader, _ = read_table(b)
                    frame = pd.DataFrame(reader)
                if 'PIP' in frame:
                    frame['source_file'] = filename+':'+info.filename
                    frames.append(frame)
                elif 'ACAT_p' in frame:
                    frame['celltype'] = Path(info.filename).name.split('_all_cis')[0]
                    p = pd.to_numeric(frame.ACAT_p, errors='coerce')
                    ok = p.between(0, 1)
                    frame['q'] = float('nan')
                    frame.loc[ok, 'q'] = false_discovery_control(p[ok])
                    genes.append(frame)
    susie = pd.concat(frames, ignore_index=True)
    significant = pd.concat(genes, ignore_index=True)
    significant = significant[significant.q <= .05][['gene', 'celltype']].drop_duplicates()
    requested = [CELL_MAP[c] for c in config['selection']['requested_cell_types']]
    if not set(requested) <= set(susie.celltype):
        raise ValueError('Requested cells absent from observed release')
    subset = susie[susie.celltype.isin(requested)].merge(significant, on=['gene', 'celltype'], validate='many_to_one')
    screened = screen(subset, config['study']['seed'], config['selection']['target_variant_gene_pairs'], config['selection']['min_pip_preferred'])
    threshold = config['selection']['min_pip_preferred']
    if len(screened) < config['selection']['target_variant_gene_pairs']:
        threshold = config['selection']['min_pip_fallback']
        screened = screen(subset, config['study']['seed'], config['selection']['target_variant_gene_pairs'], threshold)
    screened.to_csv(out/'screened_candidates.tsv', sep='\t', index=False)
    excluded = screened[['gene', 'SNP', 'celltype', 'source_file']].copy()
    excluded['reason'] = 'missing_oriented_beta;unvalidated_reference_TSS;unknown_effective_coverage'
    excluded.to_csv(out/'exclusions.tsv', sep='\t', index=False)
    pd.DataFrame(columns=COLUMNS).to_csv(out/'candidates.tsv', sep='\t', index=False)
    pd.DataFrame(columns=COLUMNS).to_csv(out/'controls.tsv', sep='\t', index=False)
    pip_numeric = pd.to_numeric(subset.PIP, errors='raise')
    high = subset[pip_numeric >= threshold]
    snv_pattern = r'(?:chr)?(?:[0-9]+|X|Y):[1-9][0-9]*:[ACGT]:[ACGT]'
    snvs = high[high.SNP.str.fullmatch(snv_pattern)]
    nonpal = snvs[~snvs.SNP.str.contains(r':(?:A:T|T:A|C:G|G:C)$', regex=True)]
    flow = pd.DataFrame([
        ('source_susie', len(susie)), ('requested_cells', int(susie.celltype.isin(requested).sum())),
        ('significant_gene_fdr_05', len(subset)), ('pip_threshold', len(high)),
        ('biallelic_snv', len(snvs)), ('nonpalindromic_snv', len(nonpal)),
        ('sampled_distinct_genes', len(screened)), ('fully_eligible', 0)], columns=['stage', 'rows'])
    flow.to_csv(out/'selection_flow.tsv', sep='\t', index=False)
    status = {'source_susie_rows': len(susie), 'requested_cell_rows': int(susie.celltype.isin(requested).sum()),
              'gene_fdr_05_rows': len(subset), 'pip_threshold': threshold,
              'screened_pairs': len(screened), 'screened_genes': screened.gene.nunique(),
              'screened_by_cell': screened.celltype.value_counts().to_dict(),
              'eligible_candidates': 0, 'eligible_controls': 0, 'reference_checked': 0,
              'excluded_screened_pairs': len(excluded), 'config_hash': digest(config),
              'note': 'Screening only; all screened rows fail downstream eligibility. Source filtering counts are aggregate, not eligible efficacy denominators.'}
    save(out/'selection.json', status)
    return status


def readiness(args, config):
    reasons = []
    if not os.getenv(config['api']['key_env']):
        reasons.append('GI_API_KEY absent')
    for field in ['observations', 'fasta', 'reference_manifest', 'quota_preflight']:
        if not getattr(args, field, None):
            reasons.append(field+' not supplied')
    return reasons


def infer(args, config, out, limit, stage='pilot'):
    reasons = readiness(args, config)
    if reasons:
        return {'status': 'blocked', 'reasons': reasons, 'n_scored': 0, 'requests': 0}
    quota = json.loads(Path(args.quota_preflight).read_text())
    # Explicit recorded cost/allowance, not a claim inferred from burst headers.
    if not quota.get('cost_approved') or int(quota.get('remaining_requests', 0)) < args.max_requests:
        raise Blocked('Quota/cost preflight insufficient')
    provenance = json.loads(Path(args.reference_manifest).read_text())
    if provenance.get('assembly') != 'GRCh38' or provenance.get('annotation_release') != 'GENCODE v44' or provenance.get('sha256') != file_hash(args.fasta):
        raise ValueError('Unverified reference provenance')
    observations = pd.read_csv(args.observations, sep='\t')
    required = set(COLUMNS)-{'model', 'status', 'beta_alt'} | {'tss1', 'strand', 'beta', 'effect_allele', 'maf', 'distance', 'gc', 'category', 'effect_allele_source'}
    if not required <= set(observations):
        raise ValueError('Incomplete verified observations: '+','.join(sorted(required-set(observations))))
    if observations[list(required)].isna().any().any() or not observations.genome_build.eq('GRCh38').all():
        raise ValueError('Missing or mismatched observation metadata')
    observations['beta_alt'] = [orient_beta(r.beta, r.effect_allele, r.ref, r.alt) for r in observations.itertuples()]
    observations['split'] = split_loci(observations, config['study']['seed'])
    contract = json.loads((out/'gi_openapi.json').read_text())
    catalog = json.loads((out/'gi_models.json').read_text())
    model = config['models']['primary']
    spec = next(m['bio_spec'] for m in catalog['models'] if m['id'] == model)
    if spec['trained_window_bp'] is not None:
        raise Blocked('Primary requires verifiable actual token coverage')
    flank = spec['recommended_flank_bp']
    predictions, exclusions = [], []
    observations = observations.drop_duplicates(['SNP', 'gene', 'cell'])
    if not observations.category.isin(['candidate', 'control']).all():
        raise ValueError('Unknown cohort category')
    if stage == 'smoke':
        nominated = {'1:36478697:C:A', '16:56971389:C:T'}
        observations = observations[observations.SNP.isin(nominated)]
        if not nominated <= set(observations.SNP):
            raise Blocked('Both protocol smoke variants require verified observations')
    else:
        observations = observations[(observations.category == 'control') |
                                    (observations.PIP >= config['selection']['min_pip_fallback'])]
    if stage == 'pilot':
        pair_ids = observations[['SNP', 'gene', 'category']].drop_duplicates()
        counts = pair_ids.category.value_counts()
        if counts.get('candidate', 0) < 100 or counts.get('control', 0) < 100:
            raise Blocked('Pilot requires 100 distinct candidates and 100 distinct controls')
        ids = pd.concat([g.head(100) for _, g in pair_ids.groupby('category')])
    else:
        ids = observations[['SNP', 'gene']].drop_duplicates().head(limit)
    chosen = observations.merge(ids, on=list(ids.columns), validate='many_to_one')
    if chosen[['SNP', 'gene']].drop_duplicates().shape[0] < limit:
        raise Blocked('Insufficient distinct variant-gene pairs for stage')
    if stage in ('trial', 'pilot') and chosen.cell.nunique() < 4:
        raise Blocked('Stage requires at least four observed cells')
    if chosen.groupby(['SNP', 'gene']).cell.nunique().max() > 4:
        raise ValueError('More than four contexts per pair')
    for row in chosen.to_dict('records'):
        chrom, pos, ref, alt = variant(row['SNP'])
        if (chrom, pos, ref, alt) != (str(row['chrom']).removeprefix('chr'), int(row['pos1']), row['ref'], row['alt']):
            raise ValueError('Variant ID and coordinate/allele columns disagree')
        if len(ref) != 1 or len(alt) != 1 or {ref, alt} in ({'A', 'T'}, {'C', 'G'}):
            raise ValueError('Primary requires nonpalindromic biallelic SNVs')
        if not 0 <= row['maf'] <= .5 or not 0 <= row['gc'] <= 1 or not 0 <= row['PIP'] <= 1:
            raise ValueError('Invalid matching covariates/PIP')
        if stage != 'smoke' and row['category'] == 'candidate' and row['PIP'] < config['selection']['min_pip_fallback']:
            raise ValueError('Candidate below prespecified fallback PIP')
    if args.dry_run:
        return {'status': 'dry_run', 'rows': len(chosen), 'worst_case_requests': len(chosen)*2}
    client = Client(config['paths']['cache'], contract, model, args.max_requests, config['api']['timeout_seconds'], config['api']['max_retries'])
    for row in chosen.to_dict('records'):
        context = config['models']['context_template'].format(cell_type=row['cell'])
        try:
            prediction = score_pair(client, args.fasta, row, flank, context)
            predictions.append(prediction)
        except (Blocked, ValueError, KeyError) as e:
            exclusions.append(dict(SNP=row['SNP'], gene=row['gene'], cell=row['cell'], reason=str(e)))
        # Every pair is persisted immediately; inference itself reuses hashed responses.
        pd.DataFrame(predictions).to_parquet(out/'predictions.parquet', index=False)
        pd.DataFrame(exclusions).to_csv(out/'inference_exclusions.tsv', sep='\t', index=False)
    chosen.to_parquet(out/'observations.parquet', index=False)
    return {'status': 'complete' if len(predictions) == len(chosen) else 'incomplete',
            'n_scored': len(predictions), 'n_attempted': len(chosen), 'n_excluded': len(exclusions), 'requests': client.calls}


def evaluation(config, out):
    if (out/'live_run_manifest.json').exists():
        from .analysis.live import summarize_live
        return summarize_live(out)
    file = out/'predictions.parquet'
    frame = pd.read_parquet(file) if file.exists() else pd.DataFrame(columns=['beta_alt', 'delta', 'locus', 'gene', 'cell'])
    if frame.empty:
        frame = pd.DataFrame({k: pd.Series(dtype='float64' if k in ['beta_alt', 'delta'] else 'object') for k in ['beta_alt', 'delta', 'locus', 'gene', 'cell']})
    a = config['analysis']
    metrics = evaluate(frame, a['near_zero_delta_epsilon'], a['bootstrap_replicates'], a['permutation_replicates'], config['study']['seed'])
    if not frame.empty and 'category' in frame:
        metrics['by_category'] = {category: evaluate(g, a['near_zero_delta_epsilon'], a['bootstrap_replicates'], a['permutation_replicates'], config['study']['seed']) for category, g in frame.groupby('category')}
        metrics['by_cell'] = {cell: evaluate(g, a['near_zero_delta_epsilon'], a['bootstrap_replicates'], a['permutation_replicates'], config['study']['seed']) for cell, g in frame.groupby('cell')}
        metrics['primary'] = metrics['by_category'].get('candidate')
    save(out/'metrics.json', metrics)
    return metrics


def main(argv=None):
    args_in = sys.argv[1:] if argv is None else argv
    if args_in and args_in[0] == "followup":
        from .followup.cli import main as followup_main
        return followup_main(args_in[1:])
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='command', required=True)
    for command in ['audit', 'select', 'smoke', 'infer', 'evaluate', 'report', 'run']:
        p = subs.add_parser(command, help=command+' benchmark stage')
        p.add_argument('--config', default='configs/pilot.yaml')
        p.add_argument('--out', type=Path, default=Path('results/pilot-20261029'))
        p.add_argument('--dry-run', action='store_true')
        p.add_argument('--max-requests', type=int, default=1600)
        p.add_argument('--observations', help='Verified normalized TSV; see docs/IMPLEMENTATION.md')
        p.add_argument('--fasta')
        p.add_argument('--reference-manifest')
        p.add_argument('--quota-preflight')
        p.add_argument('--verified-live', action='store_true', help='Use audited full-summary workflow and documented 40960-bp geometry amendment; hard cap 760 ledger attempts')
    args = parser.parse_args(argv)
    config = yaml.safe_load(Path(args.config).read_text())
    if not 0 <= args.max_requests <= config['api']['max_requests']:
        parser.error('max-requests must be within configured budget')
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    if args.command == 'run' and args.verified_live:
        from .pipeline.live import run_live
        args.max_requests = min(args.max_requests, 760)
        args.prepare_only = False
        run_live(args)
        evaluation(config, out)
        report(out)
        return 0 if json.loads((out/'stages.json').read_text())['pilot']['status'] == 'complete' else 2
    binding = digest({'config': config, 'observations': file_hash(args.observations) if args.observations else None,
                      'reference_manifest': file_hash(args.reference_manifest) if args.reference_manifest else None})
    binding_file = out/'execution_binding.json'
    if args.command in ('run', 'smoke', 'infer'):
        if binding_file.exists() and json.loads(binding_file.read_text())['hash'] != binding:
            parser.error('Output directory belongs to different inputs; use a new --out and reuse the shared cache')
        save(binding_file, {'hash': binding})
    stages = json.loads((out/'stages.json').read_text()) if (out/'stages.json').exists() else {}
    command = args.command
    try:
        if command in ['audit', 'run']:
            audit(config, out); stages['audit'] = {'status': 'complete'}
        if command in ['select', 'run']:
            select(config, out); stages['selection'] = {'status': 'provisional_only'}
        if command in ['smoke', 'infer', 'run']:
            names = ['smoke', 'trial', 'pilot'] if command == 'run' else [command]
            for name in names:
                limit = {'smoke': 2, 'trial': 20, 'pilot': 200, 'infer': 200}[name]
                if name in ('trial', 'pilot', 'infer') and stages.get('smoke', {}).get('status') != 'complete':
                    result = {'status': 'blocked', 'reasons': ['live smoke gate not passed']+readiness(args, config), 'n_scored': 0, 'requests': 0}
                elif name in ('pilot', 'infer') and stages.get('trial', {}).get('status') != 'complete':
                    result = {'status': 'blocked', 'reasons': ['trial gate not passed'], 'n_scored': 0, 'requests': 0}
                elif readiness(args, config):
                    result = {'status': 'blocked', 'reasons': readiness(args, config), 'n_scored': 0, 'requests': 0}
                else:
                    stage_out = out/name
                    stage_out.mkdir(exist_ok=True)
                    for filename in ['gi_openapi.json', 'gi_models.json']:
                        (stage_out/filename).write_bytes((out/filename).read_bytes())
                    result = infer(args, config, stage_out, limit, 'pilot' if name == 'infer' else name)
                    if (stage_out/'predictions.parquet').exists():
                        (out/'predictions.parquet').write_bytes((stage_out/'predictions.parquet').read_bytes())
                stages[name] = result
                save(out/'stages.json', stages)
        if command in ['evaluate', 'run']:
            if not (out/'predictions.parquet').exists():
                import pyarrow as pa
                import pyarrow.parquet as pq
                fields = [(k, pa.float64() if k in ('delta', 'beta_alt', 'ref_score', 'alt_score') else pa.string()) for k in COLUMNS + ['delta', 'ref_score', 'alt_score']]
                schema = pa.schema(fields, metadata={b'schema_version': b'1', b'coordinates': b'pos1=one-based;windows=zero-based-half-open', b'effect': b'genomic-ALT-minus-REF;log(TPM+1)', b'genome': b'GRCh38'})
                table = pa.Table.from_pylist([], schema=schema)
                pq.write_table(table, out/'predictions.parquet')
                pq.write_table(table, out/'observations.parquet')
                (out/'requests.jsonl').touch()
            evaluation(config, out)
        save(out/'stages.json', stages)
        if command in ['report', 'run']:
            report(out)
        print(json.dumps({'command': command, 'output': str(out), 'stages': stages}, indent=2))
        if command in ('smoke', 'infer', 'run') and any(s.get('status') in ('blocked', 'incomplete') for s in stages.values()):
            return 2
        return 0
    except (ValueError, Blocked, OSError, KeyError) as e:
        save(out/'stage_error.json', {'stage': command, 'error': str(e)})
        print(f'{command} failed: {e}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
