# pysmips

**Cached [SMIPS](https://www.tern.org.au/) daily soil-moisture cubes for
Australia — download once per pixel-day, never twice.** The Soil
Moisture Integration and Prediction System (CSIRO/TERN) publishes a
national ~1 km daily soil-water product as one Cloud-Optimised GeoTIFF
per product per day on the TERN datastore. Every pixel-day this machine
ever downloads lands in one sparse, chunk-indexed store, so repeat
requests, overlapping AOIs, extended date ranges and new products reuse
everything already fetched. Part of the
[Borevitz Lab](https://biology.anu.edu.au/research/research-groups/borevitz-group-plant-genomics-climate-adaption) ecosystem.

## How it works

```
{data_root}/smips_store/
├── index.db       # SQLite: populated (product, day, chunk) cells + upstream holes
└── smips.zarr/
    ├── totalbucket   # one sparse (time, y, x) array per product
    └── smindex ...
```

- SMIPS v1.0 sits on **one fixed grid** (4110 × 3474 px, ~0.01°,
  EPSG:4326) that has not changed since 2005-01-01, so the lattice is
  hardcoded (`pysmips.grid`). Every COG the store opens is checked
  against it; a mismatch raises instead of writing.
- Any bbox × date range maps deterministically to a set of
  (day, 128 × 128-px chunk) cells. `Store.get_ds(bbox, start, end)`
  diffs them against the ledger and downloads **only the missing
  cells**: one remote COG open per day, one integer-aligned windowed
  read per chunk — no resampling, ever. Days are fetched concurrently
  (`workers`, default 8).
- Four chunks nest in each 512 × 512 COG block, so a small AOI does not
  drag a 500 km block into the store, while adjacent chunks of one day
  share the block read.
- Reads return the **exact AOI window** (not whole chunks — a SMIPS chunk
  is ~130 km across) as an `xarray.Dataset` on `(time, lat, lon)`.
- A day the datastore does not hold (HTTP 404) is recorded as *absent*
  once it is older than `SMIPS.publish_lag_days` (30) and reads as NaN;
  a younger 404 is "not published yet" and is retried next time. Any
  other failure propagates rather than poisoning the store.
- Pixel reads require a TERN API key (listings are public) — set
  `tern_api_key` in `~/.config/Troi.json`, `TROI_TERN_KEY`, or pass
  `api_key=` per call. Keys are free from <https://account.tern.org.au/>.

## Products

| key | upstream | units | from |
|---|---|---|---|
| `totalbucket` | `totalbucket/` | mm (total profile soil water) | 2005-01-01 |
| `smindex` | `SMindex/` | 0–1 saturation fraction | 2005-01-01 |
| `bucket1` | `bucket1/` | mm | 2015-01-01 |
| `bucket2` | `bucket2/` | mm | 2015-01-01 |
| `deepD` | `deepD/` | mm | 2015-01-01 |
| `runoff` | `runoff/` | mm | 2015-01-01 |

Days before a product's first publication are skipped without a
request and read as NaN.

## Usage

The core API is **troi-agnostic** — a bbox and dates:

```python
from datetime import date
from pysmips.store import Store

store = Store()
bbox = [147.30, -35.52, 147.62, -35.10]   # [W, S, E, N]

ds = store.get_ds(bbox, date(2020, 1, 1), date(2020, 12, 31))
ds['totalbucket']                          # (time, lat, lon) DataArray, mm

ds = store.get_ds(bbox, date(2020, 1, 1), date(2020, 12, 31),
                  products=('totalbucket', 'smindex', 'bucket1'))

store.fill(bbox, date(2020, 1, 1), date(2020, 12, 31))   # → 0: already local
store.absent_days('totalbucket', date(2020, 1, 1), date(2020, 12, 31))  # upstream holes
```

Pass `log=print` to `fill`/`get_ds` for one progress line per 64-day
time chunk on long fills.

Pipelines that speak the shared `troi.Troi` use the adapters:

```python
ds = store.get_ds_troi(troi)
```

`download_smips(troi)` remains as a thin wrapper.

## Performance

Live measurements against the TERN datastore, `workers=8`:

| Scenario | Downloaded | Time |
|---|---|---|
| Cold fill — Kyeamba (32 × 42 px), 7 days | 14 cells | 2.1 s |
| Cold fill — Murrumbidgee envelope (581 × 251 px, 18 chunks/day), one year | 6 588 cells | 44 s (≈ 8 days/s) |
| Same request again | nothing | **0.03 s** |
| Sub-AOI inside a filled envelope, any dates covered | nothing | **0.01 s** |
| Read the envelope year (366 × 251 × 581) | — | 1.3 s |
| Read a 20-year point series | — | ≈ 1 s |

Store footprint: ≈ 340 MB per envelope-year of `totalbucket` (≈ 7 GB
for the full 2005–2025 record of one product over an 82 000 km²
catchment). A 20-year envelope fill is ≈ 15 min per product. Absolute
times vary with network and TERN load; the zeros are the point — they
are ledger lookups, no network involved.

## Install

### pip

```bash
pip install git+https://github.com/thestochasticman/pysmips.git
```

Dependencies (the `troi-core` core from PyPI, plus rasterio / xarray /
zarr ≥ 3) are declared in `pyproject.toml` and installed automatically.

### From source

```bash
git clone https://github.com/thestochasticman/pysmips.git
cd pysmips
pip install -e .
```

Package design (shared across the lab's packages — no inheritance,
composition only):

- **`Troi`** (from `troi`) — identity: what region, what dates.
- **`SMIPS`** (`pysmips.smips`) — config: endpoint, product catalog,
  publication lag.
- **`Paths`** (`pysmips.paths`) — derived locations of the store for a
  given `Config`.
- **`grid`** — the fixed SMIPS lattice, chunk and time-axis math (pure,
  offline-testable).
- **`Store`** (`pysmips.store`) — ties them together.

## Test

```bash
# offline (pure math + synthetic store):
python pysmips/grid.py     # True
python pysmips/smips.py    # True
python pysmips/paths.py    # True
python pysmips/store.py    # True

# live (small real reads from TERN — needs tern_api_key):
python pysmips/download_smips.py  # True
```
