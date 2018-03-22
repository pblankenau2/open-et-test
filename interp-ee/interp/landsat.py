
import datetime
import logging
import sys

import ee

import wrs2

ee.Initialize()


def get_landsat_coll(wrs2_tile_list, start_date, end_date=None, cloud_cover=100,
                     landsat5_flag=True, landsat7_flag=True,
                     landsat8_flag=True, landsat_type='TOA'):
    """Initialize Landsat collection for a date range and WRS2 tile

    This function can NOT be mapped over an Earth Engine collection

    Args:
        wrs2_tile_list (list): Landsat WRS2 tiles to include in collection
        start_date (str): ISO date format (YYYY-MM-DD)
        end_date (str): ISO Format date (YYYY-MM-DD) (inclusive/optional).
            If not set or None, end date will be one day after start date.
            Default is None.
        cloud_cover (float): Maximum ACCA cloud cover percentage.
            Filters using the cloud cover percentage value in the metadata.
            Default is 100 (return all scenes).
        landsat5_flag (bool): if True, include Landsat 4/5 images.
            Default is True.
        landsat7_flag (bool): if True, include Landsat 7 images.
            Default is True.
        landsat8_flag (bool): if True, include Landsat 8 images.
            Default is True.
        landsat_type (str): Landsat image type.
            Choices are 'TOA', 'RAD' (and eventually SR).
            Default is 'TOA'.

    Returns:
        ee.ImageCollection() of Landsat images
    """

    if start_date == end_date or not end_date:
        # Set to 1 day past start date
        start_dt = datetime.datetime.strptime(start_date, '%Y-%m-%d')
        end_date = (start_dt + datetime.timedelta(days=1)).strftime('%Y-%m-%d')

    wrs2_tile_geom = ee.Geometry.MultiPoint([
        wrs2.centroids[wrs2_tile]
        for wrs2_tile in wrs2_tile_list
        if wrs2_tile in wrs2.centroids.keys()], 'EPSG:4326')

    if landsat_type == 'TOA':
        l8_coll_str = 'LANDSAT/LC08/C01/T1_RT_TOA'
        l7_coll_str = 'LANDSAT/LE07/C01/T1_RT_TOA'
        l5_coll_str = 'LANDSAT/LT05/C01/T1_TOA'
    elif landsat_type == 'RAD':
        l8_coll_str = 'LANDSAT/LC08/C01/T1_RT'
        l7_coll_str = 'LANDSAT/LE07/C01/T1_RT'
        l5_coll_str = 'LANDSAT/LT05/C01/T1'

    # If more than one Landsat type, merge collections
    # Otherwise return a single collection
    # This is super messy, but the goal is to only build and
    #   merge collections if absolutely necessary.
    # This will probably only be helpful for generating Tcorr
    if sum([landsat8_flag, landsat7_flag, landsat5_flag]) > 1:
        landsat_coll = ee.ImageCollection([])
        landsat_coll = ee.ImageCollection(landsat_coll.merge(
            ee.ImageCollection(l8_coll_str) \
                .filterDate(start_date, end_date) \
                .filterDate('2013-03-24', datetime.datetime.today().date().isoformat()) \
                .filterBounds(wrs2_tile_geom) \
                .filterMetadata('CLOUD_COVER_LAND', 'less_than', cloud_cover) \
                .filterMetadata('DATA_TYPE', 'equals', 'L1TP')))
        landsat_coll = ee.ImageCollection(landsat_coll.merge(
            ee.ImageCollection(l7_coll_str) \
                .filterDate(start_date, end_date) \
                .filterBounds(wrs2_tile_geom) \
                .filterMetadata('CLOUD_COVER_LAND', 'less_than', cloud_cover) \
                .filterMetadata('DATA_TYPE', 'equals', 'L1TP')))
        landsat_coll = ee.ImageCollection(landsat_coll.merge(
            ee.ImageCollection(l5_coll_str) \
                .filterDate(start_date, end_date) \
                .filterBounds(wrs2_tile_geom) \
                .filterMetadata('CLOUD_COVER_LAND', 'less_than', cloud_cover) \
                .filterMetadata('DATA_TYPE', 'equals', 'L1TP')))
    elif landsat8_flag:
        landsat_coll = ee.ImageCollection(l8_coll_str) \
            .filterDate(start_date, end_date) \
            .filterBounds(wrs2_tile_geom) \
            .filterMetadata('CLOUD_COVER_LAND', 'less_than', cloud_cover) \
            .filterMetadata('DATA_TYPE', 'equals', 'L1TP')
    elif landsat7_flag:
        landsat_coll = ee.ImageCollection(l7_coll_str) \
            .filterDate(start_date, end_date) \
            .filterBounds(wrs2_tile_geom) \
            .filterMetadata('CLOUD_COVER_LAND', 'less_than', cloud_cover) \
            .filterMetadata('DATA_TYPE', 'equals', 'L1TP')
    elif landsat5_flag:
        landsat_coll = ee.ImageCollection(l5_coll_str) \
            .filterDate(start_date, end_date) \
            .filterBounds(wrs2_tile_geom) \
            .filterMetadata('CLOUD_COVER_LAND', 'less_than', cloud_cover) \
            .filterMetadata('DATA_TYPE', 'equals', 'L1TP')

    return landsat_coll


def get_landsat_image(wrs2_tile_list, start_date, end_date=None,
                      cloud_cover=100, landsat_type='TOA'):
    """Initialize Landsat image for a date range and WRS2 tile

    This function can NOT be mapped over an Earth Engine collection

    Args:
        wrs2_tile_list (list): Landsat WRS2 tiles to include in collection
        start_date (str): ISO date format (YYYY-MM-DD)
        end_date (str): ISO Format date (YYYY-MM-DD) (inclusive/optional).
            If not set or None, end date will be one day after start date.
            Default is None.
        cloud_cover (float): Maximum ACCA cloud cover percentage.
            Filters using the cloud cover percentage value in the metadata.
            Default is 100 (return all scenes).
        landsat_type (str): Landsat image type.
            Choices are 'TOA', 'RAD' (and eventually SR).
            Default is 'TOA'.

    Returns:
        Single Landsat ee.Image()
    """
    # It would make more sense to have a dedicated function for looking up a
    # Landsat image by ID instead of filtering the collection to a single image
    image_l5_flag = False
    image_l7_flag = False
    image_l8_flag = False
    if start_date >= '2013-03-24':
        image_l8_flag = True
    if start_date >= '1999-01-01':
        image_l7_flag = True
    if end_date and end_date <= '2011-12-31':
        image_l5_flag = True
    elif not end_date and start_date <= '2011-12-31':
        image_l5_flag = True

    # Build the Landsat collection
    landsat_coll = ee.ImageCollection(get_landsat_coll(
        wrs2_tile_list=wrs2_tile_list, start_date=start_date,
        end_date=end_date, cloud_cover=cloud_cover,
        landsat5_flag=image_l5_flag, landsat7_flag=image_l7_flag,
        landsat8_flag=image_l8_flag, landsat_type=landsat_type))

    return ee.Image(landsat_coll.first())


def landsat_bqa_cloud_mask_func(img):
    """Extract collection 1 CFMask cloud mask

    Output image is structured to be applied directly with updateMask()
      i.e. 0 is cloud, 1 is cloud free

    https://landsat.usgs.gov/collectionqualityband
    https://code.earthengine.google.com/356a3580096cca315785d0859459abbd

    Confidence values
    00 = "Not Determined" = Algorithm did not determine the status of this condition
    01 = "No" = Algorithm has low to no confidence that this condition exists (0-33 percent confidence)
    10 = "Maybe" = Algorithm has medium confidence that this condition exists (34-66 percent confidence)
    11 = "Yes" = Algorithm has high confidence that this condition exists (67-100 percent confidence
    """
    qa_img = ee.Image(img.select(['BQA']))

    # Extracting cloud masks from BQA using rightShift() and  bitwiseAnd()
    # Cloud (med & high confidence), then snow, then shadow, then fill
    # Low confidence clouds tend to be the FMask buffer
    cloud_mask = qa_img.rightShift(4).bitwiseAnd(1).neq(0) \
        .And(qa_img.rightShift(5).bitwiseAnd(3).gte(2)) \
        .Or(qa_img.rightShift(7).bitwiseAnd(3).gte(3)) \
        .Or(qa_img.rightShift(9).bitwiseAnd(3).gte(3)) \
        .Or(qa_img.rightShift(11).bitwiseAnd(3).gte(3))
    return cloud_mask.Not()


def parse_landsat_id(system_index):
    """Return the components of a Landsat collection 1 product ID

    LT05_PPPRRR_YYYYMMDD
    Format matches EE collection 1 system:index
    """
    sensor = system_index[0:4]
    path = int(system_index[5:8])
    row = int(system_index[8:11])
    year = int(system_index[12:16])
    month = int(system_index[16:18])
    day = int(system_index[18:20])
    return sensor, path, row, year, month, day
