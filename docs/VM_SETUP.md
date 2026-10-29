# Linux VM setup and sizing

**Recommended selection from supplied provider list:** 16 CPUs / 64 GiB RAM; 500-GB SSD; Ubuntu 24.04 LTS; no GPU. For 100-variant pilot only, 8 CPUs / 32 GiB RAM with 100–200 GB SSD may suffice. Full 170.6-GB compressed record requires considerably larger disk (2 TB or more) and possibly 32–64 CPUs / 128–256 GiB RAM for concurrent parquet processing. RAM and disk are separate: do not assume provider includes 500 GB disk. GI inference is remote, so VM vCPUs do not speed up remote API inference.

## System setup
```bash
sudo apt update
sudo apt install -y git curl wget unzip jq build-essential python3 python3-venv python3-pip samtools tabix zstd pigz tmux
df -h
free -h
nproc
git clone https://github.com/burtsev/GxI-TenK10K.git
cd GxI-TenK10K
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
# After Codex creates pyproject.toml:
pip install -e '.[dev]'
# Supply GI_API_KEY securely via shell secrets manager / environment, never commit it.
export GI_BASE_URL=https://api.genomicintelligence.ai
mkdir -p data/raw data/cache data/processed results
```

Use `tmux` for long live sessions and `df -h` before downloads. If using attached/block disk, mount it to `data/` and `results/`; allow several times the ZIP member size for decompress + parquet + intermediates. If memory consumption spikes, use DuckDB and pyarrow iterators with bounded batches and a strict temp directory.

The following are **expected commands once Codex implements the CLI**, not commands that already work:
```bash
gxi-tenk10k audit --config configs/pilot.yaml
gxi-tenk10k select --config configs/pilot.yaml
gxi-tenk10k smoke --config configs/pilot.yaml
gxi-tenk10k infer --config configs/pilot.yaml --dry-run
gxi-tenk10k run --config configs/pilot.yaml --max-requests 1600
gxi-tenk10k report --config configs/pilot.yaml
```

Credentials: `GI_API_KEY` is required for online inference. The API's exact authentication header, endpoint and request schema must be confirmed from live OpenAPI. No credentials in source code, printed logs or tracked YAML. Set `DATA_DIR` and `RESULTS_DIR` if needed; capture immutable dataset/model versions and hashes. Respect official API quotas, privacy/licensing and request costs. The tested endpoint must be a real current one, not an assumption from these instructions.
