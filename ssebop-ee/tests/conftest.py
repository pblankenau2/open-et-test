import logging

import ee
import pytest


@pytest.fixture(scope="session", autouse=True)
def test_init():
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logging.debug('Test Setup')

    ee.Initialize()

    # Make a simple EE request
    logging.debug(ee.Number(1).getInfo())
