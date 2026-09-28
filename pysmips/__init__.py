# Light-weight exports only: pysmips.store (and the download_smips wrapper)
# pull in rasterio/zarr/xarray, so those stay behind explicit submodule imports.
from pysmips.paths import Paths
from pysmips.smips import SMIPS, defaultsmips

__all__ = [
    'Paths',
    'SMIPS',
    'defaultsmips',
]
