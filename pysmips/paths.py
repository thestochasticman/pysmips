"""Derived on-disk locations of the machine-wide SMIPS store.

The store is keyed by :class:`troi.Config` (one store per data root,
shared by every request on this machine). Rule of thumb across the
lab's packages: user-settable inputs -> Config, derived locations ->
Paths. No inheritance -- composition only.
"""
from attrs import frozen, field
from troi import Config, config as default_config


@frozen
class Paths:
    """Where the pysmips store lives for a given Config.

    Attributes:
        config: The :class:`troi.Config` supplying the data root (and the
            TERN API key).
        root: Store directory (``{config.tmp_dir}/smips_store``). Cross-node
            claims live under ``{root}/claims`` (see :mod:`troi.ledger`).
        store: The sparse Zarr store -- one ``(time, y, x)`` array per product.
        ledger: Marker tree of populated blocks:
            ``ledger/{product}/{tc:05d}/{cy:03d}_{cx:03d}.json`` holding a
            0/1 string, one character per day of the time chunk.
        absent: Marker tree of upstream 404s, for the ``gaps()`` audit only:
            ``absent/{product}/{year}/{YYYY-MM-DD}.json``.

    Example:
        ```python
        from pysmips.paths import Paths

        Paths().store  # '~/Downloads/Troi-Tmp/smips_store/smips.zarr'
        ```
    """

    config: Config = default_config

    root: str = field(init=False)
    store: str = field(init=False)
    ledger: str = field(init=False)
    absent: str = field(init=False)

    root.default(lambda s: f'{s.config.tmp_dir}/smips_store')
    store.default(lambda s: f'{s.root}/smips.zarr')
    ledger.default(lambda s: f'{s.root}/ledger')
    absent.default(lambda s: f'{s.root}/absent')


def test_paths_derive_from_config():
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysmips_paths_test_')
    cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    paths = Paths(cfg)
    return (
        paths.root == f'{tmpdir}/smips_store'
        and paths.store == f'{tmpdir}/smips_store/smips.zarr'
        and paths.ledger == f'{tmpdir}/smips_store/ledger'
        and paths.absent == f'{tmpdir}/smips_store/absent'
    )


def test():
    return test_paths_derive_from_config()


if __name__ == '__main__':
    print(test())
