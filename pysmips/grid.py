"""The fixed SMIPS grid every stored pixel lives on, plus the time axis.

SMIPS v1.0 is published as one national COG per product per day on a
single EPSG:4326 lattice (4110 x 3474 pixels, ~0.01 degree, ~1 km) that
has not changed since 2005-01-01 -- every product and every year carries
the identical affine transform, so it is hardcoded here (like the
Copernicus DEM in pycopdem) rather than read from the source (like SLGA
in pyslga). The store still verifies each COG it opens against these
constants and refuses to write if TERN ever regrids.

Spatial chunks are 128 x 128 pixels: four nest exactly in each 512 x 512
COG block, so a chunk is one integer-aligned windowed read with no
resampling, and a small AOI does not drag a 500 km block into the store.
Time chunks are 64 days, so a 20-year point series reads ~115 chunks.

The time axis is fixed: day 0 is 2005-01-01 (the product's first day)
and the array runs to HORIZON. Zarr chunks that were never written cost
nothing on disk, so the long horizon is free.

All functions here are pure -- no I/O, no store access.
"""
from datetime import date, timedelta

X0, Y_TOP = 112.904998779, -9.005000113999998        # top-left corner
XRES, YRES = 0.009997566018978103, 0.009997121616580312
WIDTH, HEIGHT = 4110, 3474
COG_BLOCK = 512                                      # source COG tiling
CHUNK = 128                                          # store spatial chunk
TCHUNK = 64                                          # store time chunk (days)
EPOCH = date(2005, 1, 1)                             # day index 0
HORIZON = date(2049, 12, 31)
NDAYS = (HORIZON - EPOCH).days + 1

TRANSFORM = (XRES, 0.0, X0, 0.0, -YRES, Y_TOP)       # rasterio order (a b c d e f)


def day_index(d: date) -> int:
    """Index of ``d`` on the time axis (0 = 2005-01-01)."""
    i = (d - EPOCH).days
    if not 0 <= i < NDAYS:
        raise ValueError(f'{d} is outside the store time axis {EPOCH}..{HORIZON}')
    return i


def date_of(i: int) -> date:
    """Inverse of :func:`day_index`."""
    return EPOCH + timedelta(days=int(i))


def window_for_bbox(bbox: list[float]) -> tuple[int, int, int, int]:
    """Exact pixel window ``(row0, row1, col0, col1)`` covering ``bbox``,
    snapped outward to pixel edges and clipped to the raster.

    Unlike pyslga/pycopdem this is *not* chunk-aligned: SMIPS chunks are
    ~130 km across, and a read must return the AOI, not its chunks.
    Chunk alignment is the fill's concern (:func:`chunks_in_window`).
    """
    west, south, east, north = bbox
    col0 = max(0, int((west - X0) / XRES // 1))
    col1 = min(WIDTH, -int(-((east - X0) / XRES) // 1))
    row0 = max(0, int((Y_TOP - north) / YRES // 1))
    row1 = min(HEIGHT, -int(-((Y_TOP - south) / YRES) // 1))
    if col1 <= col0 or row1 <= row0:
        raise ValueError(f'bbox {bbox} does not intersect the SMIPS grid')
    return (row0, row1, col0, col1)


def chunks_in_window(window: tuple[int, int, int, int]) -> list[tuple[int, int]]:
    """All spatial chunk ids ``(cy, cx)`` intersecting a pixel window."""
    row0, row1, col0, col1 = window
    return [
        (cy, cx)
        for cy in range(row0 // CHUNK, -(-row1 // CHUNK))
        for cx in range(col0 // CHUNK, -(-col1 // CHUNK))
    ]


def chunk_window(cy: int, cx: int) -> tuple[int, int, int, int]:
    """Pixel window of one spatial chunk, clipped to the raster (edge
    chunks are partial)."""
    return (
        cy * CHUNK, min(HEIGHT, (cy + 1) * CHUNK),
        cx * CHUNK, min(WIDTH, (cx + 1) * CHUNK),
    )


def time_chunks(i0: int, i1: int) -> list[int]:
    """Time-chunk ids covering day indices ``[i0, i1]`` inclusive."""
    return list(range(i0 // TCHUNK, i1 // TCHUNK + 1))


def tchunk_range(tc: int) -> tuple[int, int]:
    """Day-index range ``[t0, t1)`` of a time chunk, clipped to the axis."""
    return (tc * TCHUNK, min(NDAYS, (tc + 1) * TCHUNK))


def coords_for_window(window: tuple[int, int, int, int]):
    """Pixel-centre coordinate arrays ``(lat, lon)`` for a window
    (lat descending)."""
    import numpy as np
    row0, row1, col0, col1 = window
    lon = X0 + (np.arange(col0, col1) + 0.5) * XRES
    lat = Y_TOP - (np.arange(row0, row1) + 0.5) * YRES
    return lat, lon


_BBOX = [147.30, -35.52, 147.62, -35.10]   # Kyeamba Creek


def test_window_contains_bbox():
    row0, row1, col0, col1 = window_for_bbox(_BBOX)
    west, south, east, north = _BBOX
    return (
        X0 + col0 * XRES <= west and X0 + col1 * XRES >= east
        and Y_TOP - row0 * YRES >= north and Y_TOP - row1 * YRES <= south
        and (col1 - col0) < 40 and (row1 - row0) < 50     # exact, not chunked
    )


def test_chunks_nest_in_cog_blocks():
    return COG_BLOCK % CHUNK == 0 and all(
        (chunk_window(cy, cx)[0] % CHUNK == 0 and chunk_window(cy, cx)[2] % CHUNK == 0)
        for cy, cx in chunks_in_window(window_for_bbox(_BBOX))
    )


def test_edge_chunk_is_clipped():
    r0, r1, c0, c1 = chunk_window(HEIGHT // CHUNK, WIDTH // CHUNK)
    return r1 == HEIGHT and c1 == WIDTH and r1 > r0 and c1 > c0


def test_time_axis_roundtrip():
    d = date(2010, 6, 15)
    return (
        day_index(EPOCH) == 0 and date_of(day_index(d)) == d
        and time_chunks(0, TCHUNK - 1) == [0]
        and time_chunks(TCHUNK - 1, TCHUNK) == [0, 1]
        and tchunk_range(0) == (0, TCHUNK)
    )


def test_out_of_axis_raises():
    try:
        day_index(date(2004, 12, 31))
    except ValueError:
        return True
    return False


def test_disjoint_bbox_raises():
    try:
        window_for_bbox([10.0, 50.0, 11.0, 51.0])
    except ValueError:
        return True
    return False


def test_overlapping_bboxes_share_chunks():
    a = window_for_bbox(_BBOX)
    b = window_for_bbox([_BBOX[0] + 0.05, _BBOX[1], _BBOX[2] + 0.05, _BBOX[3]])
    return len(set(chunks_in_window(a)) & set(chunks_in_window(b))) > 0


def test():
    return all([
        test_window_contains_bbox(),
        test_chunks_nest_in_cog_blocks(),
        test_edge_chunk_is_clipped(),
        test_time_axis_roundtrip(),
        test_out_of_axis_raises(),
        test_disjoint_bbox_raises(),
        test_overlapping_bboxes_share_chunks(),
    ])


if __name__ == '__main__':
    print(test())
