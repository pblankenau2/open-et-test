import ee

ee.Initialize()


class NDVI_ET(object):
    """Earth Engine DisALEXI"""

    def __init__(self, input_image,
                 m=1.25, b=0.0,
                 elevation=ee.Image('USGS/NED')):
        """Initialize an image for computing NDVI ET

        """
        self.toa_image = input_image

        # Use the SPACECRAFT_ID property identify each Landsat type
        self.spacecraft_id = ee.String(self.toa_image.get('SPACECRAFT_ID'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'BQA'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1', 'BQA'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']})
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'bqa']

        self.image = ee.Image(self.toa_image) \
            .select(input_bands.get(self.spacecraft_id), output_bands)

    def compute_etf(self):
        """"""
        return self._get_ndvi(self.image) \
            .multiply(self.m).add(self.b) \
            .rename(['etf'])

    def _get_ndvi(self):
        """Compute NDVI

        Parameters
        ----------
        self.input_image : ee.Image

        Returns
        -------
        ndvi : ee.Image

        """
        return self.input_image.normalizedDifference(['nir', 'red']) \
            .rename(['ndvi'])


class Landsat():
    def __init__()
        # Use the SPACECRAFT_ID property identify each Landsat type
        self.spacecraft_id = ee.String(self.toa_image.get('SPACECRAFT_ID'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'BQA'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1', 'BQA'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA'],
        })
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'bqa']

        self.image = ee.Image(self.toa_image) \
            .select(input_bands.get(self.spacecraft_id), output_bands)


class Sentinel():
    def __init__()
        # Use the SPACECRAFT_NAME property identify each Landsat type
        self.spacecraft_name = ee.String(self.toa_image.get('SPACECRAFT_NAME'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'Sentinel-2A': ['B2', 'B3', 'B4', 'B8', 'B7'],
            'Sentinel-2B': ['B2', 'B3', 'B4', 'B8', 'B7'],
        })
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2']

        self.image = ee.Image(self.toa_image) \
            .select(input_bands.get(self.spacecraft_id), output_bands)
