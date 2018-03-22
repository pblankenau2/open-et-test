import logging
import sys

import ee

ee.Initialize()

system_properties = ['system:index', 'system:time_start']


def interp_et_coll(et_reference_coll, et_fraction_coll, interp_days=64,
                   interp_type='linear'):
    """Generate daily ETa collection from ETo and ETf collections

    Non-mappable function.

    Args:
        et_reference_coll (ee.ImageCollection): Daily Reference ET images
        et_fraction_coll (ee.ImageCollection): Landsat ETf images
        interp_days (int):
        interp_type (str): Interpolation type.
            Choices are: "linear"
            Defaults to "linear"

    Returns:
        ee.ImageCollection() of daily ET images
    """
    # Add TIME_0UTC as a separate image band for quality mosaic
    interp_etf_coll = et_fraction_coll.map(add_time_bands)

    # Filters for joining the neighboring Landsat images in time
    # Need to add one extra day since ETo time_start may be offset from ETf
    # DEADBEEF - Does it need to be added to both filters?
    # DEADBEEF - Which one should be < / > and which should be <= / >=
    # We should probably use TIME_0UTC here instead of system:time_start
    prev_filter = ee.Filter.And(
        ee.Filter.maxDifference(
            difference=(interp_days + 1) * 24 * 60 * 60 * 1000,
            leftField='system:time_start', rightField='system:time_start'),
        ee.Filter.greaterThan(
            leftField='system:time_start', rightField='system:time_start'))
    next_filter = ee.Filter.And(
        ee.Filter.maxDifference(
            difference=(interp_days + 1) * 24 * 60 * 60 * 1000.0,
            leftField='system:time_start', rightField='system:time_start'),
        ee.Filter.lessThanOrEquals(
            leftField='system:time_start', rightField='system:time_start'))

    # Join the neighboring Landsat images in time
    et_reference_coll = ee.ImageCollection(
        ee.Join.saveAll('prev', 'system:time_start', True).apply(
            et_reference_coll, interp_etf_coll, prev_filter))
    et_reference_coll = ee.ImageCollection(
        ee.Join.saveAll('next', 'system:time_start', False).apply(
            et_reference_coll, interp_etf_coll, next_filter))

    # Interpolate
    if interp_type.lower() == 'linear':
        return ee.ImageCollection(et_reference_coll.map(linear_et_func))
    else:
        logging.error('\nERROR: Invalid interpolation type: {}'.format(
            interp_type))
        sys.exit()
        # return ee.ImageCollection([])


def linear_et_func(image):
    """Linearly interpolate daily ETf and multiply by daily reference ET

    Function should be mapped over the daily daily ETo/ETr collection
        (i.e. ETo/ETr collection)

    Image Bands:
        et_reference:

    Image Properties:
        prev: list of images with bands: etf and time
        next: list of images with bands: etf and time
    """
    et_reference_image = ee.Image(image).select('et_reference')

    time_0utc = date_to_time_0utc(ee.Date(image.get('system:time_start')))
    time_image = ee.Image.constant(time_0utc).double().rename(['time'])

    # For mosaic, joined images were sorted with closest image in time last
    prev_qm_image = ee.ImageCollection.fromImages(
        ee.List(et_reference_image.get('prev'))).mosaic()
    next_qm_image = ee.ImageCollection.fromImages(
        ee.List(et_reference_image.get('next'))).mosaic()

    prev_etf_image = ee.Image(prev_qm_image.select('etf')).double()
    next_etf_image = ee.Image(next_qm_image.select('etf')).double()
    prev_time_image = ee.Image(prev_qm_image.select('time')).double()
    next_time_image = ee.Image(next_qm_image.select('time')).double()

    # Fill masked values with values from the opposite image
    # Something like this is needed to ensure there are always two images
    #   to interpolate between
    # For large data gaps, this will cause a flat line instead of a ramp
    prev_time_mosaic = ee.Image(ee.ImageCollection.fromImages([
        next_time_image, prev_time_image]).mosaic())
    next_time_mosaic = ee.Image(ee.ImageCollection.fromImages([
        prev_time_image, next_time_image]).mosaic())
    prev_etf_mosaic = ee.Image(ee.ImageCollection.fromImages([
        next_etf_image, prev_etf_image]).mosaic())
    next_etf_mosaic = ee.Image(ee.ImageCollection.fromImages([
        prev_etf_image, next_etf_image]).mosaic())

    # Calculate time ratio of Landsat image between other cloud free
    #   Landsat images images
    time_ratio_image = time_image.subtract(prev_time_mosaic) \
        .divide(next_time_mosaic.subtract(prev_time_mosaic))

    # Interpolate NDVI values to the current image time
    interp_etf_image = next_etf_mosaic.subtract(prev_etf_mosaic) \
        .multiply(time_ratio_image).add(prev_etf_mosaic)

    return interp_etf_image.multiply(et_reference_image).select([0], ['et']) \
        .copyProperties(image, system_properties)


def aggregate_daily(image_coll, start_date, end_date, agg_type='mean'):
    """Aggregate images by day

    This function should be used to mosaic Landsat images from same path
        but different rows
    Aggregation is currently hardcoded to 'mean'
    system:time_start of returned images will be 0 UTC (not the image time)

    Args:
        image_coll (ee.ImageCollection): Input image collection.
        start_date (date/number/string): Start date.
            Needs to be an EE readable date (ISO Date string or milliseconds)
        end_date (date/number/string): End date.
            Needs to be an EE readable date (ISO Date string or milliseconds)
        agg_type (str): Aggregation type.
            Currently only 'mean' is supported for daily aggregations

    Returns:
        ee.ImageCollection()
    """
    # Build a collection of date "features" to join to
    date_list = ee.List.sequence(
        ee.Date(start_date).millis(),
        ee.Date(end_date).millis(),
        # ee.Date(end_date).advance(1, 'day').millis(),
        24 * 3600 * 1000)
    def set_date(time):
        return ee.Feature(None, {
            'system:index': ee.Date(time).format('yyyy-MM-dd'),
            'system:time_start': ee.Number(time).int64(),
            'DATE': ee.Date(time).format('yyyy-MM-dd')})

    # Add a date property to the image collection
    def set_image_date(img):
        return ee.Image(img.set(
            'DATE', ee.Date(img.get('system:time_start')).format('yyyy-MM-dd')))

    join_coll = ee.FeatureCollection(
        ee.Join.saveAll('join').apply(
            ee.FeatureCollection(date_list.map(set_date)),
            ee.ImageCollection(image_coll.map(set_image_date)),
            ee.Filter.equals(leftField='DATE', rightField='DATE')))

    def aggregate_func(ftr):
        # The composite image time will be 0 UTC (not Landsat time)
        # if agg_type.lower() == 'mean':
        return ee.Image(
            ee.ImageCollection.fromImages(ftr.get('join')).mean()
                .copyProperties(ftr, system_properties + ['DATE']))

    return ee.ImageCollection(join_coll.map(aggregate_func))


def add_time_bands(image):
    """Add TIME_0UTC as a separate image band for quality mosaic

    Mask time band with image mask
    Intentionally using TIME_0UTC (instead of system:time_start)
        so that joins and interpolation happen evenly per day
    """
    time_0utc = date_to_time_0utc(ee.Date(image.get('system:time_start')))
    return image.addBands([
        image.select([0]).double().multiply(0).add(time_0utc).rename(['time'])])


def date_to_time_0utc(date):
    """Get the 0 UTC time_start for a date

    Extra operations are needed since update() does not set milliseconds to 0.

    Args:
        date (ee.Date):

    Returns:
        ee.Number
    """
    return date.update(hour=0, minute=0, second=0).millis() \
        .divide(1000).floor().multiply(1000)
