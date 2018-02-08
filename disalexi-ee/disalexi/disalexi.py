import ee

ee.Initialize()


class DisALEXI():
    """Earth Engine DisALEXI"""

    def __init__(self, image,
                 elevation=ee.Image('USGS/NED'),
                 landcover=None, lc_type=None,
                 iterations=20
                 ):
        """Initialize an image for computing DisALEXI.

        FilterDate looks at the time_starts, so if the Alexi image
            has a start time of 0 UTC, to get the Alexi image for the
            image date, you may need to move the image date back a day.

        Parameters
        ----------
        image : ee.Image
            DisALEXI input image.
        elevation: ee.Image, optional
            Elevation [m] (the default is ee.Image('USGS/NED')).
        lc_type : {'NLCD' or 'GlobeLand30'}, optional
            Land cover type key word.
        landcover : ee.Image, optional
            Land cover image.
        iterations : int, optional
            Number of iterations of main calculation (the default is 20).

        """
        # input_image = ee.Image(image)
        self.iterations = iterations

        # CGM - Applying cloud mask directly to input image
        #   instead of to a_pt in main TSEB function
        self.cfmask = ee.Image(image).select('cfmask')
        self.mask = self.cfmask.eq(0)
        input_image = ee.Image(image).updateMask(self.mask)

        # Unpack the input bands from the image
        self.albedo = input_image.select('albedo')
        self.lai = input_image.select('lai')
        # self.lai = self.lai.where(lai.mask(), 0.01)
        self.lst = input_image.select('lst')
        self.ndvi = input_image.select('ndvi')

        # Copy system properties
        # self.index = input_image.get('system:index')
        self.time_start = input_image.get('system:time_start')

        # Set server side date/time properties using the 'system:time_start'
        self.datetime = ee.Date(self.time_start)
        self.date = ee.Date(self.datetime.format('yyyy-MM-dd'))
        self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).double()
        self.hour = ee.Number(self.datetime.getFraction('day')).multiply(24)
        self.hour_int = self.hour.floor()
        # Time used in IDL is hours and fractional minutes (no seconds)
        self.time = ee.Date(self.datetime).get('hour').add(
            ee.Date(self.datetime).get('minute').divide(60))

        # Elevation [m]
        if elevation is None:
            self.elevation = ee.Image('USGS/SRTMGL1_003').rename(['elevation'])
        else:
            self.elevation = elevation.rename(['elevation'])

        # Set default land cover image and type
        # For now default to CONUS and use default if image and type were not set
        # GlobeLand30 values need to be set to the lowest even multiple of 10,
        #   since that is currently what is in the landcover.xlsx file.
        # http://www.globallandcover.com/GLC30Download/index.aspx
        if landcover is None and lc_type is None:
            # Using NLCD as default land cover and type
            self.landcover = ee.Image('USGS/NLCD/NLCD2011').select(['landcover'])
            self.lc_type = 'NLCD'
            # Using GlobeLand30 land cover and type
            # self.landcover = ee.Image(
            #     ee.ImageCollection('users/cgmorton/GlobeLand30_2010').mosaic()) \
            #         .divide(10).floor().multiply(10) \
            #         .rename(['landcover'])
            # self.lc_type = 'GLOBELAND30'
        elif landcover is None:
            # What should happen if on only the land cover image is set?
            raise ValueError('landcover must be set if lc_type is set')
        elif lc_type is None:
            # What should happen if on only the land cover type is set?
            # The images could be looked up by type from a default LC image dict
            raise ValueError('lc_type must be set if landcover is set')
        else:
            self.landcover = landcover
            self.lc_type = lc_type

        # ALEXI ET - CONUS
        self.et_coll = ee.ImageCollection(
            'projects/climate-engine/alexi/conus/daily/et')
        self.et_transform = [0.04, 0, -125.0, 0, -0.04, 49.80]
        self.et_crs = 'EPSG:4326'
        # ALEXI ET - Global (not ingested)
        # self.et_coll = ee.ImageCollection('projects/climate-engine/alexi/global/daily/et')
        # self.et_transform = [0.05, 0, -180.0, 0, -0.05, 90]
        # self.et_crs = 'EPSG:4326'

        # Hard coding using CFSR for wind speed
        self.windspeed_coll = ee.ImageCollection('NOAA/CFSV2/FOR6H') \
            .select([
                'u-component_of_wind_height_above_ground',
                'v-component_of_wind_height_above_ground'])

        # Hard coding using MERRA2 solar insolation [W m-2]
        self.rs_hourly_coll = ee.ImageCollection(
            'projects/climate-engine/merra2/rs_hourly')
        self.rs_daily_coll = ee.ImageCollection(
            'projects/climate-engine/merra2/rs_daily')

    def ta(self):
        """"""

        self._set_alexi_et_vars()
        self._set_elevation_vars()
        self._set_landcover_vars()
        self._set_solar_vars()
        self._set_time_vars()
        self._set_weather_vars()

        # Set Tair to have the same geotransform as ALEXI ET
        t_air = self.alexi_et.multiply(0).add(self.lst) \
            .reduceResolution(reducer=ee.Reducer.mean(), maxPixels=30000) \
            .reproject(crs=self.et_crs, crsTransform=self.et_transform)

        return ee.Image(t_air.set('system:time_start', self.time_start)) \
            .rename(['t_air'])

    def et(self, t_air):
        """Compute Landsat scale DisALEXI ET

        Parameters
        ----------
        t_air: ee.Image
            ALEXI ET scale air temperature (K)

        Returns
        -------
        image : ee.Image
            DisALEXI ET image

        """
        self._set_elevation_vars()
        self._set_alexi_et_vars()
        self._set_landcover_vars()
        self._set_solar_vars()
        self._set_time_vars()
        self._set_weather_vars()

        # DEADBEEF Compute ET as 1.25 * NDVI
        et = self.ndvi.multiply(1.25)

        return ee.Image(et.set('system:time_start', self.time_start)) \
            .rename(['et'])

    def _set_alexi_et_vars(self):
        """Extract ALEXI ET image for the target image time"""
        self.alexi_et = ee.Image(ee.ImageCollection(self.et_coll) \
            .filterDate(self.date, self.date.advance(1, 'day')).first())
        self.alexi_et = self.alexi_et.rename(['alexi_et'])

    def _set_elevation_vars(self):
        """Compute elevation derived variables"""
        self.pressure = self.elevation \
            .expression(
                '101.3 * (((293.0 - 0.0065 * z) / 293.0) ** 5.26)',
                {'z': self.elevation}) \
            .rename(['pressure'])

    def _set_landcover_vars(self):
        """Compute Land Cover / LAI derived variables

        Eventually add code to fall back on default values
            aleafv: 0.9, aleafn: 0.9, aleafl: 0.9
            adeadv: 0.2, adeadn: 0.2, adeadl: 0.2
            hc_min: 0.1, hc_max: 0.5, xl: 0.5, clump: 0.99

        Parameters
        ----------
        lai : ee.Image
            Leaf area index
        landcover : ee.Image
            Landcover
        lc_type : string
            Landcover type (choices are "NLCD" or "GlobeLand30")

        """

        # DEADBEEF - Hardcoding values for testing
        self.aleafv = self.landcover.multiply(0).add(0.9)
        self.aleafn = self.landcover.multiply(0).add(0.9)
        self.aleafl = self.landcover.multiply(0).add(0.9)
        self.adeadv = self.landcover.multiply(0).add(0.2)
        self.adeadn = self.landcover.multiply(0).add(0.2)
        self.adeadl = self.landcover.multiply(0).add(0.2)
        self.leaf_width = self.landcover.multiply(0).add(0.5)
        self.clump = self.landcover.multiply(0).add(0.99)
        self.hc = self.lai.multiply(0).add(0.4)

    def _set_solar_vars(self, interpolate_flag=True):
        """Extract MERRA2 solar images for the target image time"""

        # Interpolate rs hourly image at image time
        # Hourly Rs is time average so time starts are 30 minutes early
        # Move image time 30 minutes earlier to simplify filtering/interpolation
        # This interpolation scheme will only work for hourly data
        if interpolate_flag:
            interp_dt = self.datetime.advance(-0.5, 'hour')
            # time_a = interp_time
            # time_b = interp_time
            rs_a_img = ee.Image(self.rs_hourly_coll \
                .filterDate(interp_dt.advance(-1, 'hour'), interp_dt).first())
            rs_b_img = ee.Image(self.rs_hourly_coll \
                .filterDate(interp_dt, interp_dt.advance(1, 'hour')).first())
            t_a = ee.Number(rs_a_img.get('system:time_start'))
            t_b = ee.Number(rs_b_img.get('system:time_start'))
            self.rs1 = rs_b_img.subtract(rs_a_img) \
                .multiply(interp_dt.millis().subtract(t_a).divide(t_b.subtract(t_a))) \
                .add(rs_a_img) \
                .rename(['rs'])
        else:
            self.rs1 = ee.Image(
                ee.ImageCollection(self.rs_hourly_coll \
                    .filterDate(self.date, self.date.advance(1, 'day')) \
                    .filter(ee.Filter.calendarRange(self.hour_int, self.hour_int, 'hour'))
                    .first())) \
                .rename(['rs'])
        self.rs24 = ee.Image(
            ee.ImageCollection(self.rs_daily_coll \
                .filterDate(self.date, self.date.advance(1, 'day')) \
                .first())) \
            .rename(['rs'])

    def _set_time_vars(self):
        """Compute time and position related variables

        CGM - The zs returned by this function is not used in the original
            Python code.
        The hour in the original call was the integer hour, even though it
            seems like it should be the float hour.

        """
        # DEADBEEF - Hardcoding values for testing
        self.t_end = ee.Image.constant(11.02448526443880)
        self.t_end = ee.Image.constant(26.01087501850882)
        self.sol_zenith = ee.Image.constant(0.45641128977509)

    def _set_weather_vars(self):
        """Compute weather derived variables (only wind from CFSv2 for now)

        Assume input image is a CFSv2 collection
        Assume wind speed image has two bands and compute magnitude
        It would probably make more sense to compute a single windspeed image
            in the input collection.
        For simplicity, computing the mean of all images in the UTC day

        Do we need daily, 6hr, or interpolated instantaneous data?

        CFSv2 units q: kg/kg, p: Pa, ta: K, wind: m/s
        """
        windspeed_img = ee.Image(self.windspeed_coll \
            .filterDate(self.date, self.date.advance(1, 'day')).mean())
        self.windspeed = windspeed_img \
            .expression('sqrt(b(0) ** 2 + b(1) ** 2)') \
            .rename(['windspeed'])

    @classmethod
    def fromLandsatTOA(cls, toa_image, **kwargs):
        """Constructs an NDVI_ET object from a Landsat TOA image

        Parameters
        ----------
        toa_image : ee.Image
            A raw Landsat Collection 1 TOA image.

        Returns
        -------
        DisALEXI

        """

        # Use the SPACECRAFT_ID property identify each Landsat type
        spacecraft_id = ee.String(ee.Image(toa_image).get('SPACECRAFT_ID'))

        # Rename bands to generic names
        # Rename thermal band "k" coefficients to generic names
        input_bands = ee.Dictionary({
            'LANDSAT_5': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6', 'BQA'],
            'LANDSAT_7': ['B1', 'B2', 'B3', 'B4', 'B5', 'B7', 'B6_VCID_1', 'BQA'],
            'LANDSAT_8': ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'BQA']})
        output_bands = [
            'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'lst', 'BQA']
        k1 = ee.Dictionary({
            'LANDSAT_5': 'K1_CONSTANT_BAND_6',
            'LANDSAT_7': 'K1_CONSTANT_BAND_6_VCID_1',
            'LANDSAT_8': 'K1_CONSTANT_BAND_10'})
        k2 = ee.Dictionary({
            'LANDSAT_5': 'K2_CONSTANT_BAND_6',
            'LANDSAT_7': 'K2_CONSTANT_BAND_6_VCID_1',
            'LANDSAT_8': 'K2_CONSTANT_BAND_10'})
        prep_image = ee.Image(toa_image) \
            .select(input_bands.get(spacecraft_id), output_bands) \
            .set('k1_constant', ee.Number(ee.Image(toa_image).get(k1.get(spacecraft_id)))) \
            .set('k2_constant', ee.Number(ee.Image(toa_image).get(k2.get(spacecraft_id))))

        # Build the input image
        input_image = ee.Image([
            cls._albedo(prep_image),
            cls._bqa_cfmask(prep_image),
            cls._lai(prep_image),
            cls._lst(prep_image),
            cls._ndvi(prep_image)])

        # Add properties and instantiate class
        input_image = ee.Image(input_image.setMulti({
            # 'SCENE_ID': ee.Image(toa_image).get('system:index'),
            'system:index': ee.Image(toa_image).get('system:index'),
            'system:time_start': ee.Image(toa_image).get('system:time_start')
        }))

        # Instantiate the class
        return cls(input_image, **kwargs)


    @staticmethod
    def _albedo(toa_image):
        """Compute total shortwave broadband albedo following [Liang2001]

        Parameters
        ----------
        toa_image : ee.Image

        Returns
        -------
        ee.Image

        Notes
        -----
        The Python code had the following line and comment:
            "bands = [1, 3, 4, 5, 7]  # dont use blue"
        IDL code and [Liang2001] indicate that the green band is not used.
        Coefficients were derived for Landsat 7 ETM+, but were found to be
            "suitable" to Landsat 4/5 TM also.

        References
        ----------
        .. [Liang2001] Shunlin Liang (2001),
            Narrowband to broadband conversions of land surface albedo -
            I Algorithms, Remote Sensing of Environment,
            Volume 76, Issue2, Pages 213-238,
            http://doi.org/10.1016/S0034-4257(00)00205-4
        """
        bands = ['blue', 'red', 'nir', 'swir1', 'swir2']
        coef = [0.356, 0.130, 0.373, 0.085, 0.072]
        return ee.Image(toa_image).select(bands) \
            .multiply(coef).reduce(ee.Reducer.sum()).subtract(0.0018) \
            .rename(['albedo'])

    @staticmethod
    def _bqa_cfmask(toa_image):
        """Extract CFmask from Landsat Collection 1 BQA band

        https://landsat.usgs.gov/collectionqualityband

        Confidence values
        00 = "Not Determined" = Algorithm did not determine the status of this condition
        01 = "No" = Algorithm has low to no confidence that this condition exists
            (0-33 percent confidence)
        10 = "Maybe" = Algorithm has medium confidence that this condition exists
            (34-66 percent confidence)
        11 = "Yes" = Algorithm has high confidence that this condition exists
            (67-100 percent confidence

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image.

        Returns
        -------
        ee.Image

        """
        bqa_image = ee.Image(toa_image).select(['BQA'])

        def getQABits(bqa_image, start, end, newName):
            """
            From Tyler's function
            https://ee-api.appspot.com/#97ab9a8f694b28128a5a5ca2e2df7841
            """
            pattern = 0
            for i in range(start, end + 1):
                pattern += int(2 ** i)
            return bqa_image.select([0], [newName]) \
                .bitwise_and(pattern).right_shift(start)

        # Extract the various masks from the QA band
        fill_mask = getQABits(bqa_image, 0, 0, 'designated_fill')
        # drop_mask = getQABits(bqa_image, 1, 1, 'dropped_pixel')
        # Landsat 8 only
        # terrain_mask = getQABits(bqa_image, 1, 1, 'terrain_occlusion')
        # saturation_mask = getQABits(
        #     bqa_image, 2, 3, 'saturation_confidence').gte(2)
        # cloud_mask = getQABits(bqa_image, 4, 4, 'cloud')
        cloud_mask = getQABits(bqa_image, 5, 6, 'cloud_confidence').gte(2)
        shadow_mask = getQABits(bqa_image, 7, 8, 'shadow_confidence').gte(3)
        snow_mask = getQABits(bqa_image, 9, 10, 'snow_confidence').gte(3)
        # Landsat 8 only
        # cirrus_mask = getQABits(bqa_image, 11, 12, 'cirrus_confidence').gte(3)

        # Convert masks to old style Fmask values
        # 0 - Clear land
        # 1 - Clear water
        # 2 - Cloud shadow
        # 3 - Snow
        # 4 - Cloud
        return fill_mask \
            .add(shadow_mask.multiply(2)) \
            .add(snow_mask.multiply(3)) \
            .add(cloud_mask.multiply(4)) \
            .rename(['cfmask'])

    # def _bqa_cloud_mask(toa_image):
    #     """Apply collection 1 CFMask cloud mask to a daily Landsat TOA image

    #     https://landsat.usgs.gov/collectionqualityband
    #     https://code.earthengine.google.com/356a3580096cca315785d0859459abbd

    #     Confidence values
    #     00 = "Not Determined" = Algorithm did not determine the status of this condition
    #     01 = "No" = Algorithm has low to no confidence that this condition exists (0-33 percent confidence)
    #     10 = "Maybe" = Algorithm has medium confidence that this condition exists (34-66 percent confidence)
    #     11 = "Yes" = Algorithm has high confidence that this condition exists (67-100 percent confidence

    #     """
    #     qa_img = toa_image.select(['BQA'])
    #
    #     # Extracting cloud masks from BQA using rightShift() and  bitwiseAnd().
    #     # Cloud (med & high confidence), snow, shadow, cirrus (Landsat 8), and fill.
    #     # Low confidence clouds tend to be the FMask buffer.
    #     # Use high confidence for snow, shadow and cirrus.
    #     cloud_mask = qa_img.rightShift(4).bitwiseAnd(1).neq(0) \
    #         .And(qa_img.rightShift(5).bitwiseAnd(3).gte(2)) \
    #         .Or(qa_img.rightShift(7).bitwiseAnd(3).gte(3)) \
    #         .Or(qa_img.rightShift(9).bitwiseAnd(3).gte(3)) \
    #         .Or(qa_img.rightShift(11).bitwiseAnd(3).gte(3)) \
    #         .Or(qa_img.bitwiseAnd(1).neq(0))
    #     return toa_image.updateMask(cloud_mask.Not())

    @staticmethod
    def _lai(toa_image):
        """Compute LAI using METRIC NDVI / LAI empirical equation

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image.

        Returns
        -------
        ee.Image

        """
        return DisALEXI._ndvi(toa_image) \
            .pow(3).multiply(7.0).clamp(0, 6) \
            .rename(['lai'])

    @staticmethod
    def _emissivity(toa_image):
        """Compute METRIC narrowband emissivity

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image.

        Returns
        -------
        ee.Image

        """
        lai = DisALEXI._lai(toa_image)
        ndvi = DisALEXI._ndvi(toa_image)

        # Initial values are for NDVI > 0 and LAI <= 3
        return lai.divide(300).add(0.97) \
            .where(ndvi.lte(0), 0.99) \
            .where(ndvi.gt(0).And(lai.gt(3)), 0.98)

    @staticmethod
    def _lst(toa_image):
        """Compute emissivity corrected land surface temperature (LST)
        from brightness temperature.

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image.

        Returns
        -------
        ee.Image

        Notes
        -----
        Note, the coefficients were derived from a small number of scenes in
        southern Idaho [Allen2007] and may not be appropriate for other areas.

        References
        ----------
        .. [ALlen2007a] R. Allen, M. Tasumi, R. Trezza (2007),
            Satellite-Based Energy Balance for Mapping Evapotranspiration with
            Internalized Calibration (METRIC) Model,
            Journal of Irrigation and Drainage Engineering, Vol 133(4),
            http://dx.doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)

        """
        # Get properties from image
        k1 = ee.Number(ee.Image(toa_image).get('k1_constant'))
        k2 = ee.Number(ee.Image(toa_image).get('k2_constant'))

        ts_brightness = ee.Image(toa_image).select(['lst'])
        emissivity = DisALEXI._emissivity(toa_image)

        # First back out radiance from brightness temperature
        # Then recalculate emissivity corrected Ts
        thermal_rad_toa = ts_brightness.expression(
            'k1 / (exp(k2 / ts_brightness) - 1)',
            {'ts_brightness': ts_brightness, 'k1': k1, 'k2': k2})

        # tnb = 0.866   # narrow band transmissivity of air
        # rp = 0.91     # path radiance
        # rsky = 1.32   # narrow band clear sky downward thermal radiation
        rc = thermal_rad_toa.expression(
            '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
            {
                'thermal_rad_toa': thermal_rad_toa,
                'emiss': emissivity,
                'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32})
        lst = rc.expression(
            'k2 / log(emiss * k1 / rc + 1)',
            {'emiss': emissivity, 'rc': rc, 'k1': k1, 'k2': k2})

        return lst.rename(['lst'])

    def _ndvi(toa_image):
        """Compute NDVI

        Parameters
        ----------
        toa_image : ee.Image
            Renamed TOA image with 'red' and 'nir'.

        Returns
        -------
        ee.Image

        """
        return ee.Image(toa_image).normalizedDifference(['nir', 'red']) \
            .rename(['ndvi'])
