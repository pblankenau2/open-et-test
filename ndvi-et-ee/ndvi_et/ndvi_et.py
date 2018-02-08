import ee

ee.Initialize()


class NDVI_ET():
    """GEE based model for computing ETf as a linear function of NDVI"""

    def __init__(self, image, m=1.25, b=0.0, elevation=ee.Image('USGS/NED')):
        """Initialize an NDVI_ET object.

        Parameters
        ----------
        image :
            Must have bands: "ndvi"
        m : float
            Slope (the default is 1.25).
        b : float
            Offset (the default is 0.0).
        elevation : ee.Image
            Elevation image (the default is ee.Image('USGS/NED').
            This parameter is not used in the calculation and is only be used
            to test passing EE objects to the function.

        Notes
        -----
        ETf = m * NDVI + b

        """
        input_image = ee.Image(image)
        self.ndvi = input_image.select(['ndvi'])
        self._m = m
        self._b = b
        self._elevation = elevation

    def etf(self):
        """Compute ETf"""
        return ee.Image(self.ndvi) \
            .multiply(self._m).add(self._b) \
            .rename(['etf'])

    def _ndvi(self):
        """Compute NDVI

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image with 'nir' and 'red bands.

        Returns
        -------
        ee.Image

        """
        return ee.Image(toa_image).normalizedDifference(['nir', 'red']) \
            .rename(['ndvi'])

    @classmethod
    def fromLandsatTOA(cls, toa_image, **kwargs):
        """Constructs an NDVI_ET object from a Landsat TOA image

        Parameters
        ----------
        toa_image : ee.Image

        Returns
        -------
        NDVI_ET

        """
        # Use the SPACECRAFT_ID property identify each Landsat type
        spacecraft_id = ee.String(ee.Image(toa_image).get('SPACECRAFT_ID'))

        # Rename bands to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'BQA'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1', 'BQA'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA'],
        })
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'BQA']
        prep_image = ee.Image(toa_image) \
            .select(input_bands.get(spacecraft_id), output_bands)

        # Build the input image
        # Eventually send the BQA band or a cloud mask through also
        input_image = ee.Image([
            cls._ndvi(prep_image)
        ])

        # Add properties and instantiate class
        input_image = ee.Image(input_image.setMulti({
            'system:time_start': ee.Image(toa_image).get('system:time_start')
        }))

        # Instantiate the class
        return cls(input_image, **kwargs)

    @classmethod
    def fromSentinel2TOA(cls, toa_image, **kwargs):
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
        input_bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'QA60']
        output_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'QA60']
        prep_image = ee.Image(toa_image) \
            .select(input_bands, output_bands) \
            .divide(10000.0)

        # Build the input image
        # Eventually send the BQA band or a cloud mask through also
        input_image = ee.Image([
            cls._ndvi(prep_image)
        ])

        # Add properties and instantiate class
        input_image = ee.Image(input_image.setMulti({
            'system:time_start': ee.Image(toa_image).get('system:time_start')
        }))

        # Instantiate the class
        return cls(input_image, **kwargs)

    @staticmethod
    def _ndvi(toa_image):
        """Compute NDVI

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image with 'nir' and 'red bands.

        Returns
        -------
        ee.Image

        """
        return ee.Image(toa_image).normalizedDifference(['nir', 'red']) \
            .rename(['ndvi'])