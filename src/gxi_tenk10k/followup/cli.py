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
    import yaml
    config=yaml.safe_load(Path(a.config).read_text());canonical=yaml.safe_load(Path('configs/followup.yaml').read_text())
    if config!=canonical:raise ValueError('This resolved study requires configs/followup.yaml; changes need an explicit new preregistered study')
    if a.stage=='audit':result=asyncio.run(audit(out))
    elif a.stage=='freeze':
        if not (out/'recovery/recovery_audit.json').exists():recover(out/'recovery')
        if not (out/'candidate_reserve.parquet').exists():prescreen(out)
        extract(out);result=freeze(out);freeze_analysis(out)
        from .verification import verify
        verify(out)
        from .prompts import registry
        registry(out)
        from .secondary import freeze_geometry
        freeze_geometry(out)
    elif a.stage=='score-dev':result={'rows':len(asyncio.run(score(out,'development')))}
    elif a.stage=='select-context':result=dev_select(out)
    elif a.stage=='score-test':result={'rows':len(asyncio.run(score(out,'test')))}
    elif a.stage=='diagnostics':result={'rows':len(asyncio.run(score(out,'development',diagnostics=True)))}
    elif a.stage=='geometry':
        from .secondary import geometry
        result=asyncio.run(geometry(out))
    elif a.stage=='evaluate':
        result=evaluate(out)
        from .analysis import development_inference
        from .secondary import secondary
        development_inference(out);secondary(out)
    elif a.stage=='report':
        from .report import report
        result=report(out)
    else:
        if not (out/'context_selection.json').exists():
            asyncio.run(score(out,'development'));asyncio.run(score(out,'development',diagnostics=True));dev_select(out)
        asyncio.run(score(out,'test'))
        from .secondary import geometry,secondary
        from .analysis import development_inference
        asyncio.run(geometry(out));evaluate(out);development_inference(out);secondary(out)
        from .report import report
        result=report(out)
    print(json.dumps(result,default=str,indent=2));return 0
