"""Run the verified live workflow; export GI_API_KEY before invocation."""
import argparse
from pathlib import Path
from gxi_tenk10k.pipeline.live import run_live


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('results/live-20261029'))
    parser.add_argument('--max-requests',type=int,default=760)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--config',default='configs/pilot.yaml')
    run_live(parser.parse_args())

if __name__=='__main__':main()
