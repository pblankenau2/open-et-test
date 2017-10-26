import ee

ee.Initialize()


class DisALEXI(object):
    """Earth Engine DisALEXI"""

    def __init__(self, image,
                 elevation=ee.Image('USGS/NED'),
                 landcover=None, lc_type=None,
                 iterations=20):
        """Initialize an image for computing DisALEXI

        FilterDate looks at the time_starts, so if the Alexi image
            has a start time of 0 UTC, to get the Alexi image for the
            image date, you may need to move the image date back a day.

        Parameters
        ----------
        image : ee.Image
            Prepped image
        elevation: ee.Image
            Elevation [m]
            (the default is ee.Image('USGS/NED'))
        lc_type : str
            Land cover type (choices are 'NLCD' or 'GlobeLand30')
        landcover : ee.Image
            Land cover
        iterations : int
            Number of iterations of main calculation
            (the default is 20)
        """
        # self.image = ee.Image(image)
        self.iterations = iterations

        # Set server side date/time properties using the 'system:time_start'
        self.datetime = ee.Date(ee.Image(image).get('system:time_start'))
        self.date = ee.Date(self.datetime.format('yyyy-MM-dd'))
        self.doy = ee.Number(self.datetime.getRelative('day', 'year')).add(1).double()
        self.hour = ee.Number(self.datetime.getFraction('day')).multiply(24)
        self.hour_int = self.hour.floor()
        # Time used in IDL is hours and fractional minutes (no seconds)
        self.time = ee.Date(self.datetime).get('hour').add(
            ee.Date(self.datetime).get('minute').divide(60))

        # CGM - Applying cloud mask directly to input image
        #   instead of to a_pt in main TSEB function
        self.cfmask = image.select('cfmask')
        self.mask = self.cfmask.eq(0)
        input_image = ee.Image(image).updateMask(self.mask)

        # Get input bands from the image
        self.albedo = input_image.select('albedo')
        self.lai = input_image.select('lai')
        # self.lai = self.lai.where(lai.mask(), 0.01)
        self.lst = input_image.select('lst')
        self.ndvi = input_image.select('ndvi')

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

    def compute_ta(self):
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

        return t_air

    def compute_et(self, t_air):
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

        return ee.Image(et).rename(['et'])

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
