import logging
import sys

import ee

ee.Initialize()

system_properties = ['system:index', 'system:time_start']


class SSEBop():
    def __init__(self,
                 image,
                 k_factor=1.25,
                 dt_source='ASSET',
                 elev_source='ASSET',
                 tcorr_source='SCENE',
                 tmax_source='DAYMET',
                 tdiff_buffer_value=10
                 ):
        """Initialize an image for computing SSEBop

        Args:
            image (ee.Image): A prepped SSEBop input image.
                Image must have the following bands: "ndvi" and "lst".
            k_factor (float): Scale factor for ETf (or ETo) values.
                Used to convert ETf from "ETrF" (0-1) to EToF (0-1.2) to match ETo
                (Or to convert ETo to ETr)
            dt_source (str): 'ASSET'
            elev_source (str): 'ASSET', 'GTOPO', 'NED', 'SRTM'
            tcorr_source (str, float): 'SCENE', 'MONTHLY', or a constant value
            tmax_source (str): 'DAYMET', 'GRIDMET', 'TOPOWX_MEDIAN',
                or a constant value
            tdiff_buffer_value (int): Cloud mask buffer using Tdiff (in K)
        """
        self.image = image

        # Build SCENE_ID from the system:index (possibly merged)
        scene_id = ee.String(image.get('system:index'))
        scene_id = ee.List(scene_id.split('_')).slice(-3)
        self.scene_id = ee.String(scene_id.get(0)).cat('_') \
            .cat(ee.String(scene_id.get(1))).cat('_') \
            .cat(ee.String(scene_id.get(2)))

        # Build WRS2_TILE from the scene_id
        self.wrs2_tile = ee.String('p').cat(self.scene_id.slice(5, 8)) \
            .cat('r').cat(self.scene_id.slice(8, 11))

        # Set server side date/time properties using the 'system:time_start'
        self.date = ee.Date(ee.Image(image).get('system:time_start'))
        self.year = ee.Number(self.date.get('year'))
        self.month = ee.Number(self.date.get('month'))
        self.start_date = ee.Date(self._date_to_time_0utc(self.date))
        self.end_date = self.start_date.advance(1, 'day')
        self.doy = ee.Number(self.date.getRelative('day', 'year')).add(1).int()

        # self.lst = ee.Image(self.image).select('lst')
        # self.ndvi = ee.Image(self.image).select('ndvi')

        # Input parameters
        self.dt_source = dt_source
        self.elev_source = elev_source
        self.tcorr_source = tcorr_source
        self.tmax_source = tmax_source
        self.k_factor = k_factor
        self.tdiff_buffer_value = tdiff_buffer_value

    def compute_etf(self):
        """Compute SSEBop ETf for a single image

        Apply Tdiff cloud mask buffer (mask values of 0 are set to nodata)

        """

        # Get input images and ancillary data needed to compute SSEBop ETf
        lst = ee.Image(self.image).select('lst')
        # lst = self.lst
        tcorr, tcorr_index = self._get_tcorr()
        tmax = ee.Image(self._get_tmax())
        dt = ee.Image(self._get_dt())

        # Compute SSEBop ETf
        etf = tmax \
            .expression(
                '(tmax * tcorr + dt - lst) / dt',
                {'tmax': tmax, 'dt': dt, 'lst': lst, 'tcorr': tcorr}) \
            .clamp(0, 2) \
            .multiply(self.k_factor) \
            .updateMask(tmax.subtract(lst).lte(self.tdiff_buffer_value)) \
            .rename(['etf']) \
            .copyProperties(self.image, system_properties) \
            .setMulti({'TCORR': tcorr, 'TCORR_INDEX': tcorr_index})
        return ee.Image(etf)

    def _get_dt(self):
        """"""
        if self._is_number(self.dt_source):
            dt_coll = ee.ImageCollection([
                ee.Image.constant(self.dt_source).set('DOY', self.doy)])
        elif self.dt_source.upper() == 'ASSET':
            dt_coll = ee.ImageCollection('projects/usgs-ssebop/daymet_dt_median') \
                .filter(ee.Filter.calendarRange(self.doy, self.doy, 'day_of_year'))
        else:
            logging.error('\nInvalid dT: {}\n'.format(self.dt_source))
            sys.exit()
        return ee.Image(dt_coll.first())

    def _get_elev(self):
        """"""
        if self._is_number(self.elev_source):
            elev_image = ee.Image.constant(ee.Number(self.elev_source))
        elif self.elev_source.upper() == 'ASSET':
            elev_image = ee.Image('projects/usgs-ssebop/srtm_1km')
        elif self.elev_source.upper() == 'GTOPO':
            elev_image = ee.Image('USGS/GTOPO30')
        elif self.elev_source.upper() == 'NED':
            elev_image = ee.Image('USGS/NED')
        elif self.elev_source.upper() == 'SRTM':
            elev_image = ee.Image('CGIAR/SRTM90_V4')
        else:
            logging.error('\nUnsupported elev_source: {}\n'.format(
                self.elev_source))
            sys.exit()
        return elev_image.select([0], ['elev'])

    def _get_tcorr(self):
        """Get Tcorr from pre-computed assets for each Tmax source

        This function is setup to get Tcorr based on the following priority:
          1) user specificed Tcorr value
          2) pre-compted scene specific Tcorr if available
          3) pre-computed WRS2 tile (path/row) monthly Tcorr if avaiable
          4) global default Tcorr of 0.978

        The Tcorr collections are a function the Tmax dataset that is used.
        Each Tcorr collection has an "INDEX" which specificed the priority
        (see table below) with lower values being preferred.
        The Tcorr prioritization is done my merging the collections,
        sorting based on the INDEX, and taking the first object.

        Tcorr INDEX values indicate
          0 - Scene specific Tcorr
          1 - Mean monthly Tcorr per WRS2 tile
          2 - Default Tcorr
          3 - User defined Tcorr

        Returns:
            ee.Number: Tcorr value
            ee.Number: Tcorr INDEX value
        """

        if self._is_number(self.tcorr_source):
            tcorr = ee.Number(self.tcorr_source)
            tcorr_index = ee.Number(3)
        elif (self.tcorr_source.upper() == 'SCENE' and
              self.tmax_source.upper() == 'DAYMET'):
            default_coll = ee.FeatureCollection([
                ee.Feature(None, {'INDEX': 2, 'TCORR': 0.978})])
            month_coll = ee.FeatureCollection(
                    'projects/usgs-ssebop/tcorr_daymet_monthly') \
                .filterMetadata('WRS2_TILE', 'equals', self.wrs2_tile) \
                .filterMetadata('MONTH', 'equals', self.month)
            scene_coll = ee.FeatureCollection(
                    'projects/usgs-ssebop/tcorr_daymet') \
                .filterMetadata('SCENE_ID', 'equals', self.scene_id)
            tcorr_coll = ee.FeatureCollection(
                default_coll.merge(month_coll).merge(scene_coll)).sort('INDEX')
            tcorr_ftr = ee.Feature(tcorr_coll.first())
            tcorr = ee.Number(tcorr_ftr.get('TCORR'))
            tcorr_index = ee.Number(tcorr_ftr.get('INDEX'))
        elif (self.tcorr_source.upper() == 'MONTH' and
              self.tmax_source.upper() == 'DAYMET'):
            default_coll = ee.FeatureCollection([
                ee.Feature(None, {'INDEX': 2, 'TCORR': 0.978})])
            month_coll = ee.FeatureCollection(
                'projects/usgs-ssebop/tcorr_daymet_monthly') \
                .filterMetadata('WRS2_TILE', 'equals', self.wrs2_tile) \
                .filterMetadata('MONTH', 'equals', self.month)
            tcorr_coll = ee.FeatureCollection(
                default_coll.merge(month_coll)).sort('INDEX')
            tcorr_ftr = ee.Feature(tcorr_coll.first())
            tcorr = ee.Number(tcorr_ftr.get('TCORR'))
            tcorr_index = ee.Number(tcorr_ftr.get('INDEX'))
        else:
            logging.error(
                '\nInvalid tcorr_source/tmax_source: {} / {}\n'.format(
                    self.tcorr_source, self.tmax_source))
            sys.exit()
        return tcorr, tcorr_index

    def _get_tmax(self):
        if self._is_number(self.tmax_source):
            tmax_coll = ee.ImageCollection([
                ee.Image.constant(self.tmax_source).rename(['tmax'])])
        elif self.tmax_source.upper() == 'DAYMET':
            # DAYMET does not include Dec 31st on leap years
            # Adding one extra date to end date to avoid errors
            tmax_coll = ee.ImageCollection('NASA/ORNL/DAYMET_V3') \
                .filterDate(self.start_date, self.end_date.advance(1, 'day')) \
                .select(['tmax']) \
                .map(self._c_to_k)
        else:
            logging.error('\nUnsupported tmax_source: {}\n'.format(
                self.tmax_source))
            sys.exit()

        return ee.Image(tmax_coll.first())

    # Eventually move to common or utils
    def _c_to_k(image):
        """Convert temperature from C to K"""
        return image.add(273.15) \
            .copyProperties(image, system_properties)

    def _date_to_time_0utc(date):
        """Get the 0 UTC time_start for a date

        Extra operations are needed since update() does not set milliseconds to 0.

        Args:
            date (ee.Date):

        Returns:
            ee.Number
        """
        return date.update(hour=0, minute=0, second=0).millis()\
            .divide(1000).floor().multiply(1000)

    def _is_number(x):
        try:
            float(x)
            return True
        except:
            return False
