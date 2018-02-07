import ee

ee.Initialize()


class NDVI_ET(object):
    """GEE based model for computing ETf as a linear function of NDVI"""

    def __init__(self, input_image, m=1.25, b=0.0,
                 elevation=ee.Image('USGS/NED')):
        """Create an NDVI_ET object

        Parameters
        ----------
        input_image :
            Must have "ndvi" and "qa" bands
        m : float
            Slope.
        b : float
            Offset.
        elevation : ee.Image
            Elevation image (not currently being used)

        Notes
        -----
        ETf = m * NDVI + b

        """
        # Unpack the inputs
        self.input_image = input_image
        self._m = m
        self._b = b
        self._elevation = elevation

        # Unpack the input bands?
        # self.ndvi = ee.Image(input_image).select(['ndvi'])

    def etf(self):
        """"""
        return ee.Image(self.input_image).select(['ndvi']) \
            .multiply(self._m).add(self._b) \
            .rename(['etf'])

    @classmethod
    def LandsatTOA(cls, toa_image):
        """Constructs an NDVI_ET object from a Landsat TOA image

        Parameters
        ----------
        toa_image : ee.Image

        Returns
        -------
        A NDVI_ET

        """
        # Use the SPACECRAFT_ID property identify each Landsat type
        spacecraft_id = ee.String(ee.Image(toa_image).get('SPACECRAFT_ID'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10'],
        })
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst']
        rename_image = ee.Image(toa_image) \
            .select(input_bands.get(spacecraft_id), output_bands)

        # Build the input image
        input_image = ee.Image([
            cls._ndvi(rename_image),
            ee.Image(toa_image).select(['BQA'], ['qa'])])

        # Add properties and instantiate class
        return cls(ee.Image(input_image).setMulti({
            'system:time_start': ee.Image(toa_image).get('system:time_start')
        }))

    @classmethod
    def Sentinel2TOA(cls, toa_image):
        """Constructs an NDVI_ET object from a Sentinel TOA image

        Parameters
        ----------
        toa_image : ee.Image

        Returns
        -------
        NDVI_ET

        """

        # Don't distinguish between Sentinel-2 A and B
        # Rename bands to generic names
        # Scale bands to 0-1 (from 0-10000)
        input_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
        output_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
        rename_image = ee.Image(toa_image) \
            .select(input_bands, output_bands) \
            .divide(10000.0)

        # Build the input image
        input_image = ee.Image([
            cls._ndvi(rename_image),
            ee.Image(toa_image).select(['QA60'], ['qa'])])

        # Add properties and instantiate class
        return cls(ee.Image(input_image).setMulti({
            'system:time_start': ee.Image(toa_image).get('system:time_start')
        }))

    @staticmethod
    def _ndvi(toa_image):
        """Compute NDVI from a common TOA image

        Parameters
        ----------
        toa_image : ee.Image

        Returns
        -------
        ee.Image

        """
        return ee.Image(toa_image).normalizedDifference(['nir', 'red']) \
            .rename(['ndvi'])
