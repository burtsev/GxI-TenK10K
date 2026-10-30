"""Export a fixed small, nonsecret evidence set for deliberate version control."""
import argparse
import json
import shutil
from pathlib import Path
from gxi_tenk10k.util import file_hash, save

FILES = ['report.md', 'data_audit.md', 'manifest.json', 'archive_profiles.json',
         'environment.json', 'gi_models.json', 'discovery.json', 'selection.json',
         'selection_flow.tsv', 'screened_candidates.tsv', 'candidates.tsv', 'controls.tsv',
         'exclusions.tsv', 'metrics.json', 'stages.json', 'test_results.txt',
         'smoke_reference_audit.json', 'official_source_evidence.json', 'execution_binding.json', 'command_history.json']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path('results/pilot-20261029'))
    parser.add_argument('--dest',type=Path,default=Path('docs/runs/pilot-20261029'))
    args=parser.parse_args();args.dest.mkdir(parents=True,exist_ok=True)
    for name in FILES:
        source=args.source/name
        if source.stat().st_size>100_000:raise ValueError('Evidence exceeds small-file cap: '+name)
        shutil.copyfile(source,args.dest/name)
    contract=json.loads((args.source/'gi_openapi.json').read_text())
    names=['ExpressionPredictRequest','ExpressionOptions','ExpressionData','ExpressionInput','ExpressionCounts']
    save(args.dest/'gi_contract_excerpt.json',{'original_sha256':file_hash(args.source/'gi_openapi.json'),
        'info':contract['info'], 'schemas':{k:contract['components']['schemas'][k] for k in names},
        'predict_description':contract['paths']['/v1/tasks/expression/predict']['post']['description'],
        'security_schemes':contract['components']['securitySchemes']})
    save(args.dest/'export_manifest.json',{'files':{name:file_hash(args.dest/name) for name in FILES+['gi_contract_excerpt.json']},
        'note':'Raw archives, DNA/reference snapshots, full OpenAPI/docs and empty Parquet outputs remain in ignored local directories. See report for limitations.'})
    print(args.dest)

if __name__=='__main__':main()
