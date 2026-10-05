"""One machine-wide SMIPS store that fills itself on demand.

Every SMIPS pixel-day this machine ever downloads lands in a single
sparse Zarr store, one ``(time, y, x)`` array per product on the SMIPS
national ~1 km grid (:mod:`pysmips.grid`):

    {config.tmp_dir}/smips_store/
    ├── smips.zarr/
    │   ├── totalbucket   # sparse (time, y, x); only written chunks exist
    │   └── smindex ...
    ├── ledger/           # one JSON marker per written (product, time-chunk, chunk)
    │   └── totalbucket/00085/012_034.json   # {"days": "0110...1"}  (64 chars)
    ├── absent/           # upstream 404s, dated, for the gaps() audit only
    │   └── totalbucket/2026/2026-10-04.json
    └── claims/           # cross-node mutex directories while a block is being written

The ledger is files, not a database (see ``troi/docs/ledger.md``): on
Gadi's Lustre, file locks are node-local, so SQLite is unsafe when
several PBS jobs share a store, and no server can run there. A marker
is committed by atomic rename, so a reader sees it whole or not at all;
a block is written under a :class:`troi.ledger.Claim` so two jobs never
read-modify-write the same Zarr block at once.

``Store.get_ds(bbox, start, end)`` diffs the requested (day, chunk)
cells per product against the ledger, fetches only the missing ones --
each day is one remote COG open plus one integer-aligned windowed read
per missing chunk, never a resample -- then reads the exact AOI window.
Days are fetched concurrently; Zarr writes happen once per (time-chunk,
spatial-chunk) block on the calling thread, each followed by its marker.

A day the datastore does not have (HTTP 404 -- in practice the two or
three most recent, not yet published) is not fetched: nothing enters the
ledger, it reads as NaN, and the next fill asks again. The 404 is
recorded under ``absent/`` with the time it was seen, so ``gaps()`` can
say *why* a day is empty; that record never suppresses a fetch. Any
other failure (401 from a stale key, timeout, 5xx) propagates -- after
the days that did succeed in the same block have been written -- because
recording it as done would poison the store with silent NaN.

Every COG opened is checked against the hardcoded grid; a mismatch
raises rather than writes. Pixel reads require a TERN API key
(``config.tern_api_key`` or the ``api_key`` argument). Nothing is ever
resampled here: ``get_ds`` returns native pixels with ``crs``,
``transform``, ``nodata`` and ``native_res_m`` attrs so a consumer can
regrid reproducibly.
"""
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta, timezone
from os import makedirs

import numpy as np
import pandas as pd
import xarray as xr
import zarr
from attrs import frozen, field

from troi import Config, config as default_config
from troi.ledger import Markers, Claims, Claim, ensure_array, Gap, GapReport
from troi.absent import clamp_end, record_absent, absent_age_days, gdal_no_404_cache
from pysmips import grid
from pysmips.paths import Paths
from pysmips.smips import SMIPS, defaultsmips

DEFAULT_PRODUCTS = ('totalbucket',)
_MISSING = object()       # sentinel: the day's COG is not on the datastore (404)
NATIVE_RES_M = 1000       # nominal; a SMIPS pixel is ~1.1 km N-S and 0.9-1.1 km E-W


class GridMismatch(RuntimeError):
    """A source COG is not on the hardcoded SMIPS grid."""


def _block_parts(product: str, tc: int, cy: int, cx: int) -> tuple:
    return (product, f'{tc:05d}', f'{cy:03d}_{cx:03d}')


def _absent_parts(product: str, day: date) -> tuple:
    return (product, str(day.year), str(day))


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
        store.gaps(bbox, start, end).summary()
        ```
    """

    config: Config = default_config
    smips: SMIPS = defaultsmips
    workers: int = 8                 # concurrent per-day COG reads
    lease_s: float = 600.0           # a block claim older than this is presumed abandoned
    paths: Paths = field(init=False)

    paths.default(lambda s: Paths(s.config))

    def __attrs_post_init__(s):
        makedirs(s.paths.root, exist_ok=True)

    # -- ledger -------------------------------------------------------------

    @property
    def _ledger(s) -> Markers:
        return Markers(s.paths.ledger)

    @property
    def _absent(s) -> Markers:
        return Markers(s.paths.absent)

    def _mask(s, product: str, tc: int, cy: int, cx: int) -> str:
        """0/1 string, one char per day of time chunk ``tc`` ('0' = not fetched)."""
        t0, t1 = grid.tchunk_range(tc)
        m = s._ledger.read(_block_parts(product, tc, cy, cx))
        return m['days'] if m else '0' * (t1 - t0)

    def _mark(s, product: str, tc: int, cy: int, cx: int, offsets) -> None:
        """OR the given day offsets (within the time chunk) into the block's
        marker. Call only while holding the block's claim."""
        bits = list(s._mask(product, tc, cy, cx))
        for k in offsets:
            bits[k] = '1'
        s._ledger.write(_block_parts(product, tc, cy, cx), {'days': ''.join(bits)})

    def _missing(s, product: str, tc: int, chunks, lo: int, hi: int) -> dict:
        """``{day: [chunks]}`` not yet in the ledger for day indices
        ``[lo, hi]`` of time chunk ``tc``."""
        t0, _ = grid.tchunk_range(tc)
        present = s._ledger.list(_block_parts(product, tc, 0, 0)[:-1])
        missing = {}
        for cy, cx in chunks:
            leaf = _block_parts(product, tc, cy, cx)[-1]
            mask = s._mask(product, tc, cy, cx) if leaf in present else None
            for i in range(lo, hi + 1):
                if mask is None or mask[i - t0] == '0':
                    missing.setdefault(grid.date_of(i), []).append((cy, cx))
        return missing

    def _api_key(s, api_key: str = None) -> str:
        api_key = api_key or s.config.tern_api_key
        if not api_key:
            raise ValueError(
                'Set tern_api_key in ~/.config/Troi.json or pass api_key parameter'
            )
        return api_key

    def _array(s, product: str, mode: str = 'a') -> zarr.Array:
        root = zarr.open_group(s.paths.store, mode=mode)
        if mode == 'r':
            return root[product]
        return ensure_array(
            s.paths.root, root, product,
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

        Troi-agnostic. Returns the number of cells actually downloaded.
        ``end`` is clamped to today. Days before a product's first
        publication are skipped without a request; days the datastore
        does not hold are skipped after one (and asked again next time),
        so 0 means every cell that exists upstream was already local.
        ``log`` is an optional callable receiving one progress line per
        time chunk.

        Safe to run from many processes and nodes at once: each
        (time-chunk, spatial-chunk) block is written under a
        :class:`troi.ledger.Claim`.
        """
        for p in products:
            s.smips.check(p)
        end = clamp_end(end)
        if end < start:
            return 0
        window = grid.window_for_bbox(bbox)
        chunks = grid.chunks_in_window(window)
        fetched = 0
        for product in products:
            p_start = max(start, s.smips.first_day(product))
            if end < p_start:
                continue
            i0, i1 = grid.day_index(p_start), grid.day_index(end)
            for tc in grid.time_chunks(i0, i1):
                t0, t1 = grid.tchunk_range(tc)
                lo, hi = max(i0, t0), min(i1, t1 - 1)
                if not s._missing(product, tc, chunks, lo, hi):
                    continue                     # full hit: no claim, no key, no network
                key = s._api_key(api_key)
                arr = s._array(product)
                keys = [(product, tc, cy, cx) for cy, cx in chunks]
                with Claims(s.paths.root, keys, lease_s=s.lease_s) as claims:
                    missing = s._missing(product, tc, chunks, lo, hi)   # re-diff under the claim
                    if not missing:
                        continue
                    n = s._fill_block(arr, product, tc, missing, key, claims)
                fetched += n
                if log:
                    log(f'{product} {grid.date_of(lo)}..{grid.date_of(hi)}: '
                        f'{n} cells over {len(missing)} days')
        return fetched

    def _fill_block(s, arr, product: str, tc: int, missing: dict, api_key: str,
                    claims: Claims) -> int:
        """Fetch every missing day of one time chunk concurrently, then
        write each touched spatial chunk once and mark it in the ledger.

        A 404 day is recorded under ``absent/`` and skipped. Any other
        per-day failure is held back until the days that did succeed are
        written and marked, then re-raised."""
        t0, t1 = grid.tchunk_range(tc)
        got, errors = {}, []
        with claims.keepalive(), ThreadPoolExecutor(max_workers=s.workers) as ex:
            futures = {ex.submit(s._fetch_day_result, product, d, missing[d], api_key): d
                       for d in sorted(missing)}
            for f in as_completed(futures):
                d, r = futures[f], f.result()
                if r is _MISSING:
                    record_absent(s._absent, _absent_parts(product, d), unit=str(d))
                elif isinstance(r, Exception):
                    errors.append((d, r))
                else:
                    got[d] = r

        fetched = 0
        touched = {c for d in got for c in missing[d]}
        for cy, cx in sorted(touched):
            r0, r1, c0, c1 = grid.chunk_window(cy, cx)
            block = arr[t0:t1, r0:r1, c0:c1]
            offsets = []
            for d, data in got.items():
                if (cy, cx) not in data:
                    continue
                k = grid.day_index(d) - t0
                block[k] = data[(cy, cx)]
                offsets.append(k)
            arr[t0:t1, r0:r1, c0:c1] = block          # zarr 3 commits each chunk by rename
            s._mark(product, tc, cy, cx, offsets)      # ...then the ledger says so
            fetched += len(offsets)
        if errors:
            d, e = errors[0]
            raise RuntimeError(
                f'{len(errors)} day(s) of {product} failed in time chunk {tc} '
                f'(first: {d}); {fetched} cells from other days were written') from e
        return fetched

    def _fetch_day_result(s, product, day, chunks, api_key):
        """:meth:`_fetch_day`, returning the exception instead of raising so
        one bad day cannot abort its siblings in the pool."""
        try:
            return s._fetch_day(product, day, chunks, api_key)
        except Exception as e:                        # noqa: BLE001 -- re-raised by caller
            return e

    def _fetch_day(s, product: str, day: date, chunks: list, api_key: str):
        """One remote COG open + one windowed read per requested chunk.
        Returns ``{(cy, cx): float32 array}``, or ``_MISSING`` on a 404 --
        the only failure that means "no data" rather than "try again"."""
        import rasterio
        from rasterio.windows import Window
        url = s.smips.url(product, day)
        env = dict(GDAL_HTTP_HEADERS=f'x-api-key: {api_key}',
                   GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
                   GDAL_HTTP_MAX_RETRY='4', GDAL_HTTP_RETRY_DELAY='2',
                   **gdal_no_404_cache(s.smips.base_url))   # a 404 must not stick for the process
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
            # A standalone 404 ("HTTP response code: 404"), not the digits of
            # a date inside the URL (…_20240404.tif) of some other failure.
            if re.search(r'(?<!\d)404(?!\d)', str(e)):
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

    # -- audit --------------------------------------------------------------

    def gaps(s, bbox: list[float], start: date, end: date,
             products=DEFAULT_PRODUCTS, today: date = None) -> GapReport:
        """Which (product, day, chunk) cells of the request are not in the
        store, and why. Touches no network.

        Statuses: ``before_product_start`` and ``after_today`` are facts
        about the request; ``absent_upstream`` means the last fill was told
        404 (``age_days`` since); ``claimed_in_progress`` means another
        process holds the block right now; ``never_fetched`` is the one
        that should be zero after a fill.
        """
        for p in products:
            s.smips.check(p)
        today = today or datetime.now(timezone.utc).date()
        chunks = grid.chunks_in_window(grid.window_for_bbox(bbox))
        ndays = (end - start).days + 1
        gaps = []
        for product in products:
            first = s.smips.first_day(product)
            absent = s._absent
            for tc in grid.time_chunks(grid.day_index(max(start, grid.EPOCH)),
                                       grid.day_index(min(end, grid.HORIZON))):
                t0, t1 = grid.tchunk_range(tc)
                lo = max(grid.day_index(max(start, grid.EPOCH)), t0)
                hi = min(grid.day_index(min(end, grid.HORIZON)), t1 - 1)
                masks = {c: s._mask(product, tc, *c) for c in chunks}
                claimed = {c for c in chunks
                           if os.path.isdir(Claim(s.paths.root, (product, tc, *c)).dir)}
                for i in range(lo, hi + 1):
                    d = grid.date_of(i)
                    absent_marker = None
                    for c in chunks:
                        if masks[c][i - t0] == '1':
                            continue
                        unit = (product, str(d), *c)
                        if d < first:
                            gaps.append(Gap(unit, 'before_product_start'))
                        elif d > today:
                            gaps.append(Gap(unit, 'after_today'))
                        else:
                            if absent_marker is None:
                                absent_marker = absent.read(_absent_parts(product, d)) or False
                            if absent_marker:
                                gaps.append(Gap(unit, 'absent_upstream',
                                                age_days=absent_age_days(absent_marker, today)))
                            elif c in claimed:
                                gaps.append(Gap(unit, 'claimed_in_progress'))
                            else:
                                gaps.append(Gap(unit, 'never_fetched'))
        return GapReport(expected=len(products) * ndays * len(chunks), gaps=tuple(gaps))

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
            window on the native grid, never resampled -- and one variable
            per product. Days the datastore does not hold (or before a
            product's first day, or after today) are NaN. Attrs carry the
            georeferencing contract: ``crs``, ``transform`` (six affine
            numbers of this window), ``nodata`` and ``native_res_m``.
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
        transform = (grid.XRES, 0.0, grid.X0 + col0 * grid.XRES,
                     0.0, -grid.YRES, grid.Y_TOP - row0 * grid.YRES)
        return xr.Dataset(
            data_vars, coords={'time': time, 'lat': lat, 'lon': lon},
            attrs={'crs': 'EPSG:4326', 'transform': list(transform), 'nodata': None,
                   'native_res_m': NATIVE_RES_M,
                   'source': 'SMIPS v1.0 (CSIRO/TERN)', 'url': s.smips.base_url},
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

    def gaps_troi(s, troi, products=DEFAULT_PRODUCTS) -> GapReport:
        """:meth:`gaps` for a :class:`troi.Troi`."""
        return s.gaps(troi.bbox, troi.start, troi.end, products=products)


# -- offline tests (synthetic cells, no network) ----------------------------

_TEST_BBOX = [147.30, -35.52, 147.62, -35.10]   # Kyeamba Creek


def _tmp_store(**config_kw) -> Store:
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysmips_store_test_')
    return Store(config=Config(out_dir=tmpdir, tmp_dir=tmpdir, **config_kw))


def _days(start: date, end: date) -> list[date]:
    return [grid.date_of(i) for i in range(grid.day_index(start), grid.day_index(end) + 1)]


def _prime(store: Store, product: str, bbox, start: date, end: date,
           value: float = 1.0, mark: bool = True):
    """Populate bbox's chunks for every day in [start, end] directly,
    bypassing the network. ``mark=False`` leaves the ledger untouched
    (simulates a crash between the Zarr write and the marker)."""
    arr = store._array(product)
    window = grid.window_for_bbox(bbox)
    days = _days(start, end)
    for cy, cx in grid.chunks_in_window(window):
        r0, r1, c0, c1 = grid.chunk_window(cy, cx)
        for d in days:
            arr[grid.day_index(d), r0:r1, c0:c1] = value        # flat field per day
        if mark:
            for tc in grid.time_chunks(grid.day_index(start), grid.day_index(end)):
                t0, t1 = grid.tchunk_range(tc)
                offs = [grid.day_index(d) - t0 for d in days if t0 <= grid.day_index(d) < t1]
                store._mark(product, tc, cy, cx, offs)


def _fake_fetch(value: float = 7.0, missing_days=(), calls=None, fail_days=()):
    """A stand-in for ``Store._fetch_day`` with no network."""
    def fake(self, product, day, chunks, api_key):
        if calls is not None:
            calls.append(day)
        if day in missing_days:
            return _MISSING
        if day in fail_days:
            raise OSError(f'synthetic failure for {day}')
        return {c: np.full((grid.CHUNK, grid.CHUNK), value, 'float32') for c in chunks}
    return fake


class _patched:
    """``with _patched(Store, '_fetch_day', fake):`` -- restore afterwards."""
    def __init__(self, cls, name, value):
        self.cls, self.name, self.value = cls, name, value

    def __enter__(self):
        self.real = getattr(self.cls, self.name)
        setattr(self.cls, self.name, self.value)

    def __exit__(self, *exc):
        setattr(self.cls, self.name, self.real)


def _leftovers(store: Store) -> list[str]:
    """Temp files and claim directories that should not survive a fill."""
    out = []
    for base, dirs, files in os.walk(store.paths.root):
        out += [f for f in files if f.endswith('.tmp')]
        if base.endswith('/claims'):
            out += dirs
    return out


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


def test_get_ds_attrs_follow_the_contract():
    """transform/crs attrs describe exactly the returned window."""
    store = _tmp_store()
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1))
    ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1))
    a, b, c, d, e, f = ds.attrs['transform']
    return (
        ds.attrs['crs'] == 'EPSG:4326' and ds.attrs['nodata'] is None
        and abs(c + 0.5 * a - float(ds.lon[0])) < 1e-9        # first pixel centre
        and abs(f + 0.5 * e - float(ds.lat[0])) < 1e-9
        and abs(a - grid.XRES) < 1e-12 and abs(e + grid.YRES) < 1e-12
        and ds.attrs['native_res_m'] == NATIVE_RES_M
    )


def test_fill_skips_populated_cells():
    store = _tmp_store()
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 31))
    return store.fill(_TEST_BBOX, date(2010, 1, 5), date(2010, 1, 20)) == 0


def test_missing_day_is_recorded_absent_and_still_reasked():
    """A 404 day (simulated) writes nothing to the ledger, reads as NaN,
    leaves a dated absent marker, and is asked for again on the next
    fill; fetched days are recorded."""
    store = _tmp_store(tern_api_key='synthetic')
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2), 5.0)
    calls = []
    with _patched(Store, '_fetch_day', _fake_fetch(7.0, missing_days=(date(2010, 1, 3),), calls=calls)):
        n1 = store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 4))
        n2 = store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 4))
        ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 4))   # re-asks Jan 3 again
    nchunks = len(grid.chunks_in_window(grid.window_for_bbox(_TEST_BBOX)))
    absent = store._absent.read(_absent_parts('totalbucket', date(2010, 1, 3)))
    return (
        n1 == nchunks                                    # only Jan 4 was fetched
        and n2 == 0                                      # ...and Jan 3 re-asked but not counted
        and sorted(calls[:2]) == [date(2010, 1, 3), date(2010, 1, 4)]   # concurrent: any order
        and calls[2:] == [date(2010, 1, 3), date(2010, 1, 3)]
        and np.isnan(ds['totalbucket'][2]).all()
        and float(ds['totalbucket'][3, 0, 0]) == 7.0
        and float(ds['totalbucket'][1, 0, 0]) == 5.0
        and absent['status'] == 404 and absent['unit'] == '2010-01-03'
        and _leftovers(store) == []
    )


def test_one_failed_day_does_not_discard_its_siblings():
    """A non-404 failure raises, but only after the other days of the block
    are written and marked; the failed day stays never_fetched."""
    store = _tmp_store(tern_api_key='synthetic')
    with _patched(Store, '_fetch_day', _fake_fetch(3.0, fail_days=(date(2010, 1, 2),))):
        try:
            store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 3))
            return False
        except RuntimeError as e:
            raised = isinstance(e.__cause__, OSError)
    report = store.gaps(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 3))
    nchunks = len(grid.chunks_in_window(grid.window_for_bbox(_TEST_BBOX)))
    return (raised and report.counts['never_fetched'] == nchunks
            and len(report.gaps) == nchunks
            and all(g.unit[1] == '2010-01-02' for g in report.gaps)
            and _leftovers(store) == [])


def test_crash_before_marker_is_refetched():
    """Zarr data without a marker is not trusted: the day is asked again
    and then marked."""
    store = _tmp_store(tern_api_key='synthetic')
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1), 9.0, mark=False)
    before = store.gaps(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1))
    calls = []
    with _patched(Store, '_fetch_day', _fake_fetch(2.0, calls=calls)):
        store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1))
    after = store.gaps(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1))
    ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 1))
    return (before.counts['never_fetched'] > 0 and calls == [date(2010, 1, 1)]
            and after.complete and float(ds['totalbucket'][0, 0, 0]) == 2.0)


def test_end_is_clamped_to_today():
    store = _tmp_store(tern_api_key='synthetic')
    today = datetime.now(timezone.utc).date()
    calls = []
    with _patched(Store, '_fetch_day', _fake_fetch(1.0, calls=calls)):
        store.fill(_TEST_BBOX, today, today + timedelta(days=30))
    return calls == [today]


def _worker_fill(tmp_dir, bbox, start, end, value):
    store = Store(config=Config(out_dir=tmp_dir, tmp_dir=tmp_dir, tern_api_key='synthetic'))
    with _patched(Store, '_fetch_day', _fake_fetch(value)):
        store.fill(bbox, start, end)


def test_two_processes_fill_overlapping_ranges():
    """Two processes fill overlapping bboxes and dates of the same store at
    once. Afterwards every cell is marked exactly once, the data is
    intact, and no claim or temp file is left behind."""
    import multiprocessing as mp
    store = _tmp_store(tern_api_key='synthetic')
    tmp = store.config.tmp_dir
    west = [147.30, -35.52, 147.50, -35.10]
    east = [147.40, -35.52, 147.62, -35.10]
    ctx = mp.get_context('fork')
    ps = [ctx.Process(target=_worker_fill, args=(tmp, west, date(2010, 1, 1), date(2010, 1, 10), 1.0)),
          ctx.Process(target=_worker_fill, args=(tmp, east, date(2010, 1, 5), date(2010, 1, 15), 2.0))]
    for p in ps:
        p.start()
    for p in ps:
        p.join(timeout=300)
    report = store.gaps(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 15))
    ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 15))
    vals = ds['totalbucket'].values
    return (all(p.exitcode == 0 for p in ps) and report.complete
            and not np.isnan(vals).any() and set(np.unique(vals)) <= {1.0, 2.0}
            and _leftovers(store) == [])


def test_gaps_classifies_every_status():
    store = _tmp_store()
    today = date(2026, 10, 5)
    _prime(store, 'bucket1', _TEST_BBOX, date(2015, 1, 3), date(2015, 1, 3))
    record_absent(store._absent, _absent_parts('bucket1', date(2015, 1, 2)), unit='2015-01-02')
    held = Claim(store.paths.root, ('bucket1', grid.time_chunks(grid.day_index(date(2015, 1, 4)),
                                                                 grid.day_index(date(2015, 1, 4)))[0],
                                     *grid.chunks_in_window(grid.window_for_bbox(_TEST_BBOX))[0])).acquire()
    try:
        r = store.gaps(_TEST_BBOX, date(2014, 12, 31), date(2015, 1, 4), products=('bucket1',), today=today)
    finally:
        held.release()
    nchunks = len(grid.chunks_in_window(grid.window_for_bbox(_TEST_BBOX)))
    c = r.counts
    r_future = store.gaps(_TEST_BBOX, today, today.replace(day=6), products=('totalbucket',), today=today)
    ages = {g.age_days for g in r.gaps if g.status == 'absent_upstream'}
    return (
        r.expected == 5 * nchunks
        and c['before_product_start'] == nchunks            # 2014-12-31
        and c['absent_upstream'] == nchunks                 # 2015-01-02
        and c['claimed_in_progress'] == 2                   # the held chunk, on 01-01 and 01-04
        and c['never_fetched'] == 2 * nchunks - 2           # the other chunks of those two days
        and all(a >= 0 for a in ages)
        and r_future.counts['after_today'] == nchunks and r_future.counts['never_fetched'] == nchunks
        and not r.complete
    )


def test_products_are_independent():
    """A populated totalbucket must not satisfy a smindex request."""
    store = _tmp_store()
    _prime(store, 'totalbucket', _TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2))
    r = store.gaps(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 2), products=('smindex',))
    return len(r.gaps) == r.expected and r.counts['never_fetched'] == r.expected


def test_pre_publication_days_need_no_network():
    """bucket1 begins 2015: a 2010 request is satisfied (all NaN) offline."""
    store = _tmp_store()
    n = store.fill(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 10), products=('bucket1',))
    ds = store.get_ds(_TEST_BBOX, date(2010, 1, 1), date(2010, 1, 10), products=('bucket1',))
    return n == 0 and bool(np.isnan(ds['bucket1']).all())


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
        test_get_ds_attrs_follow_the_contract(),
        test_fill_skips_populated_cells(),
        test_missing_day_is_recorded_absent_and_still_reasked(),
        test_one_failed_day_does_not_discard_its_siblings(),
        test_crash_before_marker_is_refetched(),
        test_end_is_clamped_to_today(),
        test_two_processes_fill_overlapping_ranges(),
        test_gaps_classifies_every_status(),
        test_products_are_independent(),
        test_pre_publication_days_need_no_network(),
        test_unknown_product_raises(),
        test_missing_key_raises_before_network(),
        test_grid_check_rejects_regrid(),
    ])


if __name__ == '__main__':
    print(test())
