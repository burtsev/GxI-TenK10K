"""Export measured sequence-free evidence and figure files for review/commit."""
import json,shutil
from pathlib import Path
from gxi_tenk10k.util import save,file_hash


def main():
    source=Path('results/live-20261029');dest=Path('docs/runs/live-20261029');dest.mkdir(parents=True,exist_ok=True)
    files=['report.md','measured_predictions.tsv','candidates.tsv','controls.tsv','exclusions.tsv','metrics.json',
           'trial_metrics.json','stages.json','full_summary_audit.json','reference_manifest.json','live_run_manifest.json',
           'environment.json','geometry_amendment.json','execution_amendments.json','comparison_stage.json','comparison_exclusions.tsv',
           'trial_cell_contrasts.tsv','cell_contrast_metrics.json','smoke_geometry_40960.json','smoke_geometry_20480.json',
           'gi_models.json','test_results.txt','command_history.json','manifest.json']
    import pandas as pd
    comparisons=pd.read_parquet(source/'trial_comparison_predictions.parquet')
    comparisons.to_csv(source/'trial_comparison_predictions.tsv',sep='\t',index=False)
    files.append('trial_comparison_predictions.tsv')
    for name in files:
        p=source/name
        if p.stat().st_size>500_000:raise ValueError('Evidence exceeds size cap: '+name)
        shutil.copyfile(p,dest/name)
    history=source/'initial_unmatched_pass'
    histdest=dest/'initial_unmatched_pass';histdest.mkdir(exist_ok=True)
    for name in ['metrics.json','stages.json','trial_metrics.json','exclusions.tsv']:
        shutil.copyfile(history/name,histdest/name)
    historical=pd.read_parquet(history/'predictions.parquet')
    # Keep required provenance last: trailing empty TSV fields are valid data,
    # but placing a nonmissing column last avoids whitespace-tool false alarms.
    historical[[c for c in historical.columns if c!='source_file']+['source_file']].to_csv(histdest/'measured_predictions.tsv',sep='\t',index=False)
    (dest/'figures').mkdir(exist_ok=True)
    for name in ['benchmark.png','diagnostics.png']:
        shutil.copyfile(source/'figures'/name,dest/'figures'/name)
    save(dest/'export_manifest.json',{'files':{str(p.relative_to(dest)):file_hash(p) for p in sorted(dest.rglob('*')) if p.is_file() and p.name!='export_manifest.json'},
        'note':'Sequence-free observed/predicted effects and aggregate evidence only. Raw sources, FASTA, archive subsets and API cache remain ignored.'})
    print(dest)

if __name__=='__main__':main()
