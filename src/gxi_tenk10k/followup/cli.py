"""CLI for the independent follow-up; historical CLI stays immutable."""
import argparse,asyncio,json
from pathlib import Path
from ..util import save
from .recovery import recover
from .data import prescreen,extract
from .freeze import freeze
from .audit import audit
from .runner import score,export,freeze_analysis
from .analysis import dev_select,evaluate

def main(argv):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['audit','freeze','score-dev','select-context','score-test','evaluate','report','resume','diagnostics','geometry'])
    p.add_argument('--out',type=Path,default=Path('results/followup-20261030'));p.add_argument('--config',default='configs/followup.yaml');a=p.parse_args(argv);out=a.out
    if a.stage=='audit':result=asyncio.run(audit(out))
    elif a.stage=='freeze':
        if not (out/'recovery/recovery_audit.json').exists():recover(out/'recovery')
        if not (out/'candidate_reserve.parquet').exists():prescreen(out)
        extract(out);result=freeze(out);freeze_analysis(out)
        from .prompts import registry
        registry(out)
    elif a.stage=='score-dev':result={'rows':len(asyncio.run(score(out,'development')))}
    elif a.stage=='select-context':result=dev_select(out)
    elif a.stage=='score-test':result={'rows':len(asyncio.run(score(out,'test')))}
    elif a.stage=='diagnostics':result={'rows':len(asyncio.run(score(out,'development',diagnostics=True)))}
    elif a.stage=='geometry':
        from .secondary import geometry
        result=asyncio.run(geometry(out))
    elif a.stage=='evaluate':result=evaluate(out)
    elif a.stage=='report':
        from .report import report
        result=report(out)
    else:
        if not (out/'context_selection.json').exists():asyncio.run(score(out,'development'));dev_select(out)
        asyncio.run(score(out,'test'));evaluate(out)
        from .report import report
        result=report(out)
    print(json.dumps(result,default=str,indent=2));return 0
