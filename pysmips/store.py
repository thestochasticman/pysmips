"""One machine-wide SMIPS store that fills itself on demand.

Every SMIPS pixel-day this machine ever downloads lands in a single
sparse Zarr store, one ``(time, y, x)`` array per product on the SMIPS
national ~1 km grid (:mod:`pysmips.grid`):

    {config.tmp_dir}/smips_store/
    ├── index.db       # SQLite ledger: fetched (product, day, chunk) cells
    └── smips.zarr/
        ├── totalbucket   # sparse (time, y, x); only written chunks exist
        └── smindex ...

``Store.get_ds(bbox, start, end)`` diffs the requested (day, chunk)
cells per product against the ledger, fetches only the missing ones --
each day is one remote COG open plus one integer-aligned windowed read
per missing chunk, never a resample -- then reads the exact AOI window.
Days are fetched concurrently; Zarr writes happen once per (time-chunk,
spatial-chunk) block on the calling thread.

A day the datastore does not have (HTTP 404 -- in practice only the two
or three most recent, not yet published) is simply not fetched: nothing
is recorded, it reads as NaN, and the next fill asks again. Any other
failure (401 from a stale key, timeout, 5xx) propagates -- recording it
as done would poison the store with silent NaN.

Every COG opened is checked against the hardcoded grid; a mismatch
raises rather than writes. Pixel reads require a TERN API key
(``config.tern_api_key`` or the ``api_key`` argument).
"""
import sqlite3
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timedelta, timezone
from os import makedirs

import numpy as np
import pandas as pd
import xarray as xr
import zarr
from attrs import frozen, field

from troi import Config, config as default_config
from pysmips import grid
from pysmips.paths import Paths
from pysmips.smips import SMIPS, defaultsmips

_SCHEMA = """
CREATE TABLE IF NOT EXISTS cells (
    product    TEXT NOT NULL,
    day        TEXT NOT NULL,
    cy         INTEGER NOT NULL,
    cx         INTEGER NOT NULL,
    written_at TEXT NOT NULL,
    PRIMARY KEY (product, day, cy, cx)
) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS cells_by_product_day ON cells(product, day);
"""

DEFAULT_PRODUCTS = ('totalbucket',)
_MISSING = object()       # sentinel: the day's COG is not on the datastore (404)


class GridMismatch(RuntimeError):
    """A source COG is not on the hardcoded SMIPS grid."""


@frozen
class Store:
    """The machine-wide SMIPS store: one grid, one ledger, zero re-downloads.

    Composed from :class:`troi.Config` (where the store lives, and the
    TERN API key) and :class:`pysmips.smips.SMIPS` (endpoint + product
    catalog). No inheritance.

    Example:
        ```python
        from datetime import date
        from pysmips.store import Store

        store = Store()
        ds = store.get_ds(bbox, date(2020, 1, 1), date(2020, 12, 31))
        ds = store.get_ds(bbox, start, end, products=('totalbucket', 'smindex'))
        ds['totalbucket']          # (time, lat, lon) DataArray, mm
        ```
    """

    config: Config = default_config
    smips: SMIPS = defaultsmips
    workers: int = 8                 # concurrent per-day COG reads
    paths: Paths = field(init=False)

    paths.default(lambda s: Paths(s.config))

    def __attrs_post_init__(s):
        makedirs(s.paths.root, exist_ok=True)

    def _db(s) -> sqlite3.Connection:
        db = sqlite3.connect(s.paths.index_db)
        db.execute('PRAGMA journal_mode=WAL')
        db.executescript(_SCHEMA)
        return db

    def _api_key(s, api_key: str = None) -> str:
        api_key = api_key or s.config.tern_api_key
        if not api_key:
            raise ValueError(
                'Set tern_api_key in ~/.config/Troi.json or pass api_key parameter'
            )
        return api_key

    def _array(s, product: str, mode: str = 'a') -> zarr.Array:
        root = zarr.open_group(s.paths.store, mode=mode)
        try:
            return root[product]
        except KeyError:
            return root.create_array(
                product,
                shape=(grid.NDAYS, grid.HEIGHT, grid.WIDTH),
                chunks=(grid.TCHUNK, grid.CHUNK, grid.CHUNK),
                dtype='float32',
                fill_value=np.nan,
            )

    # -- fill -------------------------------------------------------------

    def fill(s, bbox: list[float], start: date, end: date,
             products=DEFAULT_PRODUCTS, api_key: str = None, log=None) -> int:
        """Ensure every (day, chunk) cell of every requested product
        covering ``bbox`` x ``[start, end]`` is populated.

        Troi-agnostic. Returns the number of cells actually downloaded --
        0 means the request was already fully covered and no network was
        touched. Days before a product's first publication are skipped
        without a request; days the datastore does not hold are skipped
        after one. ``log`` is an optional callable receiving one progress
        line per time chunk.
        """
        for p in products:
            s.smips.check(p)
        window = grid.window_for_bbox(bbox)
        chunks = grid.chunks_in_window(window)
        fetched = 0
        db = s._db()
        try:
            for product in products:
                p_start = max(start, s.smips.first_day(product))
                if end < p_start:
                    continue
                i0, i1 = grid.day_index(p_start), grid.day_index(end)
                done = set(db.execute(
                    'SELECT day, cy, cx FROM cells WHERE product = ? AND day BETWEEN ? AND ?',
                    (product, str(p_start), str(end)),
                ).fetchall())
                missing = {}                     # day -> [chunks]
                for i in range(i0, i1 + 1):
                    d = grid.date_of(i)
                    need = [c for c in chunks if (str(d), *c) not in done]
                    if need:
                        missing[d] = need
                if not missing:
                    continue
                key = s._api_key(api_key)      # only now: no key needed for a full hit
                arr = s._array(product)
                for tc in grid.time_chunks(i0, i1):
                    t0, t1 = grid.tchunk_range(tc)
                    days = [d for d in missing if t0 <= grid.day_index(d) < t1]
                    if not days:
                        continue
                    n = s._fill_block(db, arr, product, tc, days, missing, key)
                    fetched += n
                    if log:
                        log(f'{product} {grid.date_of(t0)}..{grid.date_of(t1 - 1)}: '
                            f'{n} cells over {len(days)} days')
            return fetched
        finally:
            db.close()

    def _fill_block(s, db, arr, product: str, tc: int, days: list[date],
                    missing: dict, api_key: str) -> int:
        """Fetch every missing day of one time chunk concurrently, then
        write each touched spatial chunk once and record the ledger."""
        t0, t1 = grid.tchunk_range(tc)
        with ThreadPoolExecutor(max_workers=s.workers) as ex:
            got = dict(zip(days, ex.map(
                lambda d: s._fetch_day(product, d, missing[d], api_key), days)))

        now = datetime.now(timezone.utc).isoformat()
        rows, fetched = [], 0
        touched = {c for d in days for c in missing[d] if got[d] is not _MISSING}
        for cy, cx in sorted(touched):
            r0, r1, c0, c1 = grid.chunk_window(cy, cx)
            block = arr[t0:t1, r0:r1, c0:c1]
            for d in days:
                data = got[d]
                if data is _MISSING or (cy, cx) not in data:
                    continue
                block[grid.day_index(d) - t0] = data[(cy, cx)]
                rows.append((product, str(d), cy, cx, now))
                fetched += 1
            arr[t0:t1, r0:r1, c0:c1] = block
        with db:
            db.executemany(
                'INSERT OR REPLACE INTO cells (product, day, cy, cx, written_at) '
                'VALUES (?, ?, ?, ?, ?)', rows)
        return fetched

    def _fetch_day(s, product: str, day: date, chunks: list, api_key: str):
        """One remote COG open + one windowed read per requested chunk.
        Returns ``{(cy, cx): float32 array}``, or ``_MISSING`` on a 404 --
        the only failure that means "no data" rather than "try again"."""
        import rasterio
        from rasterio.windows import Window
        url = s.smips.url(product, day)
        env = dict(GDAL_HTTP_HEADERS=f'x-api-key: {api_key}',
                   GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
                   GDAL_HTTP_MAX_RETRY='4', GDAL_HTTP_RETRY_DELAY='2')
        try:
            with rasterio.Env(**env), rasterio.open(url) as src:
                s._check_grid(src, url)
                out = {}
                for cy, cx in chunks:
                    r0, r1, c0, c1 = grid.chunk_window(cy, cx)
                    data = src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0)).astype('float32')
                    if src.nodata is not None:
                        data = np.where(data == src.nodata, np.nan, data)
                    out[(cy, cx)] = data
                return out
        except rasterio.errors.RasterioIOError as e:
            if '404' in str(e):
                return _MISSING
            raise

    @staticmethod
    def _check_grid(src, url: str) -> None:
        t = tuple(src.transform)[:6]
        if (src.width, src.height) != (grid.WIDTH, grid.HEIGHT) or \
                any(abs(a - b) > 1e-9 for a, b in zip(t, grid.TRANSFORM)):
            raise GridMismatch(
                f'{url} is {src.width}x{src.height} with transform {t}; the store '
                f'grid is {grid.WIDTH}x{grid.HEIGHT} {grid.TRANSFORM}. Refusing to '
                f'write -- SMIPS has been regridded upstream.')

    # -- read -------------------------------------------------------------

    def get_ds(s, bbox: list[float], start: date, end: date,
               products=DEFAULT_PRODUCTS, api_key: str = None, log=None) -> xr.Dataset:
        """Return the SMIPS cube for ``bbox`` x ``[start, end]``, downloading
        only what's missing first.

        Troi-agnostic -- the data layer of the package. Pipelines that
        speak :class:`troi.Troi` use :meth:`get_ds_troi`.

        Args:
            bbox: ``[west, south, east, north]`` in EPSG:4326.
            start, end: Inclusive dates.
            products: Store keys from :class:`pysmips.smips.SMIPS`.
            api_key: TERN API key; falls back to ``config.tern_api_key``.

        Returns:
            xarray.Dataset with dims ``(time, lat, lon)`` -- the exact AOI
            window on the native grid -- and one variable per product.
            Days the datastore does not hold (or before a product's first
            day) are NaN.
        """
        s.fill(bbox, start, end, products=products, api_key=api_key, log=log)
        window = grid.window_for_bbox(bbox)
        row0, row1, col0, col1 = window
        i0, i1 = grid.day_index(start), grid.day_index(end)
        lat, lon = grid.coords_for_window(window)
        time = pd.date_range(start, end, freq='D')
        data_vars = {}
        root = zarr.open_group(s.paths.store, mode='a')
        for product in products:
            if product in root:
                block = root[product][i0:i1 + 1, row0:row1, col0:col1]
            else:                                   # nothing ever fetched (e.g. all pre-2015)
                block = np.full((i1 - i0 + 1, row1 - row0, col1 - col0), np.nan, 'float32')
            data_vars[product] = (('time', 'lat', 'lon'), block,
                                  {'units': s.smips.units(product)})
        return xr.Dataset(
            data_vars, coords={'time': time, 'lat': lat, 'lon': lon},
            attrs={'crs': 'EPSG:4326', 'source': 'SMIPS v1.0 (CSIRO/TERN)',
                   'url': s.smips.base_url},
        )

    # -- Troi adapters (the reproducibility layer speaks Troi) ----------

    def fill_troi(s, troi, products=DEFAULT_PRODUCTS, api_key: str = None, log=None) -> int:
        """:meth:`fill` for a :class:`troi.Troi`."""
        return s.fill(troi.bbox, troi.start, troi.end, products=products,
                      api_key=api_key, log=log)

    def get_ds_troi(s, troi, products=DEFAULT_PRODUCTS, api_key: str = None,
                    log=None) -> xr.Dataset:
        """:meth:`get_ds` for a :class:`troi.Troi`."""
        return s.get_ds(troi.bbox, troi.start, troi.end, products=products,
                        api_key=api_key, log=log)


# -- offline tests (synthetic cells, no network) ----------------------------

_TEST_BBOX = [147.30, -35.52, 147.62, -35.10]   # Kyeamba Creek


def _tmp_store() -> Store:
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysmips_store_test_')
    return Store(config=Config(out_dir=tmpdir, tmp_dir=tmpdir))


def _prime(store: Store, product: str, bbox, start: date, end: date,
           value: float = 1.0):
    """Populate bbox's chunks for every day in [start, end] directly,
    bypassing the network."""
    arr = store._array(product)
    window = grid.window_for_bbox(bbox)
    days = [start + timedelta(days=k) for k in range((end - start).days + 1)]
    rows = []
    for cy, cx in grid.chunks_in_window(window):
        r0, r1, c0, c1 = grid.chunk_window(cy, cx)
        for d in days:
            arr[grid.day_index(d), r0:r1, c0:c1] = value        # flat field per day
            rows.append((product, str(d), cy, cx, 'synthetic'))
    db = store._db()
    with db:
        db.executemany(
            'INSERT OR REPLACE INTO cells (product, day, cy, cx, written_at) '
            'VALUES (?, ?, ?, ?, ?)', rows)
    db.close()


def test_synthetic_write_read_roundtrip():
    store = _tmp_store()
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 5), 123.0)
    ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 5))
    return (
        ds['totalbucket'].shape[0] == 5
        and float(ds['totalbucket'][0, 0, 0]) == 123.0
        and ds.lat[0] > ds.lat[-1]
        and ds['totalbucket'].shape[1] < grid.CHUNK       # exact AOI, not whole chunks
        and ds['totalbucket'].attrs['units'] == 'mm'
    )


def test_fill_skips_populated_cells():
    store = _tmp_store()
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 31))
    return store.fill(_TEST_BBOX, date(2010, 1, 5), date(2010, 1, 20)) == 0


def test_missing_day_is_skipped_not_recorded():
    """A 404 day (simulated) writes nothing to the ledger, reads as NaN,
    and is asked for again on the next fill; fetched days are recorded."""
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysmips_store_test_')
    store = Store(config=Config(out_dir=tmpdir, tmp_dir=tmpdir, tern_api_key='synthetic'))
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2), 5.0)
    calls = []

    def fake_fetch(self, product, day, chunks, api_key):
        calls.append(day)
        return _MISSING if day == date(2010, 1, 3) else {
            c: np.full((grid.CHUNK, grid.CHUNK), 7.0, 'float32') for c in chunks}

    real = Store._fetch_day
    setattr(Store, '_fetch_day', fake_fetch)
    try:
        n1 = store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 4))
        n2 = store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 4))
        ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 4))   # re-asks Jan 3 again
    finally:
        setattr(Store, '_fetch_day', real)
    nchunks = len(grid.chunks_in_window(grid.window_for_bbox(_TEST_BBOX)))
    return (
        n1 == nchunks                                    # only Jan 4 was fetched
        and n2 == 0                                      # ...and Jan 3 re-asked but not counted
        and calls == [date(2010, 1, 3), date(2010, 1, 4), date(2010, 1, 3), date(2010, 1, 3)]
        and np.isnan(ds['totalbucket'][2]).all()
        and float(ds['totalbucket'][3, 0, 0]) == 7.0
        and float(ds['totalbucket'][1, 0, 0]) == 5.0
    )


def test_products_are_independent():
    """A populated totalbucket must not satisfy a smindex request."""
    store = _tmp_store()
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2))
    db = store._db()
    row = db.execute("SELECT 1 FROM cells WHERE product = 'smindex'").fetchone()
    db.close()
    return row is None


def test_pre_publication_days_need_no_network():
    """bucket1 begins 2015: a 2010 request is satisfied (all NaN) offline."""
    store = _tmp_store()
    n = store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 10), products=('bucket1',))
    ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 10), products=('bucket1',))
    return n == 0 and np.isnan(ds['bucket1']).all()


def test_unknown_product_raises():
    store = _tmp_store()
    try:
        store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2), products=('vibes',))
    except ValueError:
        return True
    return False


def test_missing_key_raises_before_network():
    store = _tmp_store()  # config with no tern_api_key
    try:
        store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2))
    except ValueError as e:
        return 'tern_api_key' in str(e)
    return False


def test_grid_check_rejects_regrid():
    class Fake:
        width, height = grid.WIDTH, grid.HEIGHT
        transform = (grid.XRES * 2, 0.0, grid.X0, 0.0, -grid.YRES, grid.Y_TOP)
    try:
        Store._check_grid(Fake(), 'fake://')
    except GridMismatch:
        return True
    return False


def test():
    return all([
        test_synthetic_write_read_roundtrip(),
        test_fill_skips_populated_cells(),
        test_missing_day_is_skipped_not_recorded(),
        test_products_are_independent(),
        test_pre_publication_days_need_no_network(),
        test_unknown_product_raises(),
        test_missing_key_raises_before_network(),
        test_grid_check_rejects_regrid(),
    ])


if __name__ == '__main__':
    print(test())
