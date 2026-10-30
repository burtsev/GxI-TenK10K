"""Reporting-only I/O workaround for the frozen QA manifest serializer.

The frozen QA rewrites an expanding response-hash dictionary on every iteration.
Suppress intermediate writes and retain the identical final sorted-key contents.
All source bindings, response hashes, checks and scientific metrics are preserved.
"""
import json
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
from gxi_tenk10k.followup import qa as qa_module
from gxi_tenk10k.followup.report import report
from gxi_tenk10k.followup.binding import validate_selection
from gxi_tenk10k.util import file_hash,save

def main(out=Path('results/followup-20261030')):
    out=Path(out);validate_selection(out);metrics_hash=file_hash(out/'metrics.json');original_save=qa_module.save;state={}
    def final_only(path,value):
        if Path(path).name=='response_hashes.json':
            if 'expected' not in state:state['expected']=pd.read_parquet(out/'historical_registry_for_future_studies.parquet',columns=['cache_key']).cache_key.nunique()
            if len(value)<state['expected']:return
        return original_save(path,value)
    qa_module.save=final_only
    try:result=report(out)
    finally:qa_module.save=original_save
    assert file_hash(out/'metrics.json')==metrics_hash
    hashes=json.load(open(out/'response_hashes.json'));assert len(hashes)==state['expected']
    note={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'scope':'Reporting-only output serialization optimization after completion of all scientific analyses','frozen_source_files_changed':False,'scientific_metrics_sha256_before_and_after':metrics_hash,'response_hash_count':len(hashes),'final_response_hash_manifest_sha256':file_hash(out/'response_hashes.json'),'workaround_source_sha256':file_hash(__file__),'method':'Suppress intermediate full-manifest writes; execute every frozen QA hash calculation/check and write final dictionary once.'}
    save(out/'reporting_io_amendment.json',note);save('docs/runs/followup-20261030/reporting_io_amendment.json',note)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
