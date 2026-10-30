import json
import subprocess
import sys
from pathlib import Path
import pandas as pd
import pytest
import yaml
from gxi_tenk10k.cli import main


@pytest.mark.parametrize('command',['audit','select','smoke','infer','evaluate','report','run'])
def test_cli_help(command):
    r=subprocess.run([str(Path(sys.executable).parent/'gxi-tenk10k'),command,'--help'],capture_output=True,text=True)
    assert r.returncode==0 and '--config' in r.stdout


def test_missing_credentials_stages(tmp_path,monkeypatch):
    monkeypatch.delenv('GI_API_KEY',raising=False)
    assert main(['smoke','--out',str(tmp_path)])==2
    stages=json.loads((tmp_path/'stages.json').read_text())
    assert stages['smoke']['status']=='blocked' and stages['smoke']['requests']==0
    assert main(['infer','--out',str(tmp_path),'--dry-run'])==2
    stages=json.loads((tmp_path/'stages.json').read_text())
    assert stages['infer']['status']=='blocked'
    assert main(['evaluate','--out',str(tmp_path)])==0
    m=json.loads((tmp_path/'metrics.json').read_text())
    assert m['n_scored']==0 and m['sign_concordance'] is None


def test_empty_stage_cannot_succeed(tmp_path,monkeypatch):
    from argparse import Namespace
    from gxi_tenk10k.cli import infer
    from gxi_tenk10k.gi.client import Blocked
    from gxi_tenk10k.util import file_hash
    monkeypatch.setenv('GI_API_KEY','synthetic-offline')
    fasta=tmp_path/'ref.fa';fasta.write_text('>1\nACGT\n')
    provenance=tmp_path/'reference.json';provenance.write_text(json.dumps({'assembly':'GRCh38','annotation_release':'GENCODE v44','sha256':file_hash(fasta)}))
    quota=tmp_path/'quota.json';quota.write_text(json.dumps({'cost_approved':True,'remaining_requests':1600}))
    observations=tmp_path/'observations.tsv'
    cols=['gene','locus','cell','SNP','chrom','pos1','ref','alt','effect_allele','PIP','source_file','genome_build','tss1','strand','beta','maf','distance','gc','category','effect_allele_source']
    pd.DataFrame(columns=cols).to_csv(observations,sep='\t',index=False)
    (tmp_path/'gi_openapi.json').write_text(json.dumps({'info':{'version':'synthetic'}}))
    (tmp_path/'gi_models.json').write_text(json.dumps({'models':[{'id':'g0-expression-8192','bio_spec':{'trained_window_bp':None,'recommended_flank_bp':40960}}]}))
    args=Namespace(observations=str(observations),fasta=str(fasta),reference_manifest=str(provenance),quota_preflight=str(quota),max_requests=1600,dry_run=False)
    config=yaml.safe_load(open('configs/pilot.yaml'));config['paths']['cache']=str(tmp_path/'cache')
    with pytest.raises(Blocked,match='100 distinct'):
        infer(args,config,tmp_path,200)
