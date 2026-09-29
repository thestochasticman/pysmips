"""Fetch the SMIPS cube for a troi -- via the machine-wide store.

Thin compatibility wrapper: the heavy lifting (ledger diffing,
concurrent windowed COG reads, missing-day handling) lives in
:class:`pysmips.store.Store`. Kept as a module so a familiar
``download_smips(troi)`` entry point exists.
"""
import xarray as xr
from troi import Troi
from pysmips.smips import SMIPS, defaultsmips
from pysmips.store import DEFAULT_PRODUCTS

def download_smips(troi: Troi, products=DEFAULT_PRODUCTS, api_key: str = None,
                   smips: SMIPS = defaultsmips, log=None) -> xr.Dataset:
    """Return the SMIPS ``(time, lat, lon)`` cube for ``troi``.

    Fetches only the (day, chunk) cells no previous request has
    populated -- repeat, overlapping and extended queries re-download
    nothing.

    Args:
        troi: The :class:`troi.Troi` (bbox + inclusive date range).
        products: Store keys (default ``('totalbucket',)``); see
            :class:`pysmips.smips.SMIPS` for the catalog.
        api_key: TERN API key; falls back to ``config.tern_api_key``.
        smips: Endpoint/catalog configuration; defaults to the bundled one.
        log: Optional callable for one progress line per time chunk.

    Returns:
        xarray.Dataset with one variable per product on the native ~1 km grid.
    """
    from pysmips.store import Store
    store = Store(config=troi.config, smips=smips)
    return store.get_ds_troi(troi, products=products, api_key=api_key, log=log)


def test_live_fetch_and_dedup():
    """Live: cold fetch covers a week of totalbucket; repeat, overlapping
    and shifted requests fetch nothing; extending fetches only the
    extension; smindex is a 0-1 fraction; the grid check passes."""
    import time
    import numpy as np
    import tempfile
    from datetime import date
    from troi import Config, config as global_config
    from pysmips.store import Store

    tmpdir = tempfile.mkdtemp(prefix='pysmips_live_test_')
    cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir, tern_api_key=global_config.tern_api_key)
    store = Store(config=cfg)
    bbox = [147.30, -35.52, 147.62, -35.10]          # Kyeamba Creek, ~32 x 42 px
    start, end = date(2020, 6, 1), date(2020, 6, 7)

    t = time.time()
    fetched = store.fill(bbox, start, end)
    cold = time.time() - t
    if fetched < 7:                                   # >= one chunk per day
        return False
    ds = store.get_ds(bbox, start, end)
    tb = ds['totalbucket'].values
    if tb.shape[0] != 7 or not np.isfinite(tb).any() or not (0 < np.nanmean(tb) < 2000):
        return False
    if store.fill(bbox, start, end) != 0:                              # identical repeat
        return False
    if store.fill([bbox[0] + 0.05, bbox[1], bbox[2] + 0.05, bbox[3]], start, end) != 0:
        return False                                                   # shifted, same chunks
    ext = store.fill(bbox, start, date(2020, 6, 9))                    # extension only
    if not (0 < ext <= 2 * len(grid_chunks(bbox))):
        return False
    smi = store.get_ds(bbox, start, end, products=('smindex',))['smindex'].values
    if not (np.nanmin(smi) >= 0 and np.nanmax(smi) <= 1):
        return False
    print(f'  cold fill {fetched} cells / 7 days in {cold:.1f}s; store at {tmpdir}')
    return True

def grid_chunks(bbox):
    from pysmips import grid
    return grid.chunks_in_window(grid.window_for_bbox(bbox))


def test():
    from troi import config
    if not config.tern_api_key:
        print('SKIPPED: set tern_api_key in ~/.config/Troi.json '
              '(or TROI_TERN_KEY) to run the live suite')
        return None
    return test_live_fetch_and_dedup()


if __name__ == '__main__':
    result = test()
    if result is not None:
        print(result)
