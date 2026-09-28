from datetime import date
from attrs import frozen


@frozen
class SMIPS:
    """Endpoint and product catalog for SMIPS v1.0 on the TERN datastore.

    The Soil Moisture Integration and Prediction System (CSIRO/TERN)
    publishes one national ~1 km COG per product per day. Pixel reads
    require a TERN API key (``x-api-key`` header); directory listings
    are public. Filenames are deterministic, so nothing is resolved from
    listings.

    Products (store key -> directory, filename stem, units, first day):

    ==============  ===========  ======================  ========  ==========
    key             directory    stem                    units     from
    ==============  ===========  ======================  ========  ==========
    ``totalbucket`` totalbucket  smips_totalbucket_mm    mm        2005-01-01
    ``smindex``     SMindex      smips_smi_perc          fraction  2005-01-01
    ``bucket1``     bucket1      smips_bucket1_mm        mm        2015-01-01
    ``bucket2``     bucket2      smips_bucket2_mm        mm        2015-01-01
    ``deepD``       deepD        smips_deepD_mm          mm        2015-01-01
    ``runoff``      runoff       smips_runoff_mm         mm        2015-01-01
    ==============  ===========  ======================  ========  ==========

    ``smindex`` is named ``perc`` upstream but holds a 0-1 saturation
    fraction (checked against the delivered files).
    """

    base_url: str = 'https://data.tern.org.au/model-derived/smips/v1_0'

    # Days younger than this may simply not be published yet: a 404 for
    # them is not recorded as absent, so the next fill retries. Older
    # 404s are genuine holes and are recorded so they are never re-asked.
    publish_lag_days: int = 30

    products = {
        'totalbucket': ('totalbucket', 'smips_totalbucket_mm', 'mm', date(2005, 1, 1)),
        'smindex':     ('SMindex',     'smips_smi_perc',       'fraction', date(2005, 1, 1)),
        'bucket1':     ('bucket1',     'smips_bucket1_mm',     'mm', date(2015, 1, 1)),
        'bucket2':     ('bucket2',     'smips_bucket2_mm',     'mm', date(2015, 1, 1)),
        'deepD':       ('deepD',       'smips_deepD_mm',       'mm', date(2015, 1, 1)),
        'runoff':      ('runoff',      'smips_runoff_mm',      'mm', date(2015, 1, 1)),
    }

    def check(self, product: str) -> str:
        if product not in self.products:
            raise ValueError(
                f'Unknown SMIPS product {product!r}. Known: {sorted(self.products)}'
            )
        return product

    def units(self, product: str) -> str:
        return self.products[self.check(product)][2]

    def first_day(self, product: str) -> date:
        return self.products[self.check(product)][3]

    def url(self, product: str, day: date) -> str:
        """COG URL for one product-day."""
        directory, stem, _, _ = self.products[self.check(product)]
        return f'{self.base_url}/{directory}/{day.year}/{stem}_{day:%Y%m%d}.tif'


defaultsmips = SMIPS()


def test_url_pattern():
    s = SMIPS()
    return (
        s.url('totalbucket', date(2005, 1, 1))
        == 'https://data.tern.org.au/model-derived/smips/v1_0/totalbucket/2005/smips_totalbucket_mm_20050101.tif'
        and s.url('smindex', date(2015, 7, 1)).endswith('/SMindex/2015/smips_smi_perc_20150701.tif')
    )


def test_unknown_product_raises():
    try:
        SMIPS().url('vibes', date(2020, 1, 1))
    except ValueError:
        return True
    return False


def test():
    return all([test_url_pattern(), test_unknown_product_raises()])


if __name__ == '__main__':
    print(test())
