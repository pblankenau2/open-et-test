#--------------------------------
# Name:         interpolate_image_export.py
# Purpose:      Export interpolated ETa/ETf/ETo/count images
#--------------------------------

import argparse
from builtins import input
import datetime
from functools import reduce
import glob
import json
import logging
import math
from operator import mul
import os
import pprint
import re
import subprocess
import sys
import time

import ee
from osgeo import ogr, osr

# Import interpolation functions
# This is an awful way of getting the parent folder into the path
# We really should package this up as a module with a setup.py
open_et_test_path = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(__file__)))))
print(open_et_test_path)
sys.path.insert(0, os.path.join(open_et_test_path, 'interp-ee', 'interp'))
import inputs
import interpolate
import landsat
import utils

# Import ET Models
sys.path.insert(0, os.path.join(open_et_test_path, 'eeflux-ee'))
import eeflux


def main(ini_path=None, overwrite_flag=False):
    """Export annual ET/ETrF/ETr/count images

    Parameters
    ----------
    ini_path : str
        Input file path.
    overwrite_flag : bool, optional
        If True, overwrite existing files (the default is False)

    Returns
    -------
    None

    """
    logging.info('\nExport annual ET/ETrF/ETr/count image')

    # Read config file
    ini = inputs.read(ini_path)
    inputs.parse_section(ini, section='INPUTS')
    inputs.parse_section(ini, section='INTERPOLATE')
    inputs.parse_section(ini, section='EXPORT')
    inputs.parse_section(ini, section=ini['INPUTS']['et_model'])

    # Use study area spatial reference if not set explicitly in INI
    if ini['EXPORT']['output_osr'] is None:
        # Get output coordinate system from study area shapefile
        study_area_ds = ogr.Open(ini['INPUTS']['study_area_path'], 0)
        study_area_lyr = study_area_ds.GetLayer()
        ini['EXPORT']['output_osr'] = osr.SpatialReference()
        ini['EXPORT']['output_osr'] = study_area_lyr.GetSpatialRef()
        ini['EXPORT']['output_crs'] = str(
            ini['EXPORT']['output_osr'].ExportToWkt())
        study_area_ds = None
        del study_area_lyr, study_area_ds
        logging.debug('\n  {:16s} {}'.format(
            'Output crs:', ini['EXPORT']['output_crs']))


    logging.debug('\nInitializing Earth Engine')
    ee.Initialize()

    # Get current running tasks
    tasks = utils.get_ee_tasks()

    # Get list of existing images/files
    if ini['EXPORT']['export_dest'] == 'ASSET':
        logging.debug('\nGetting EE asset list')
        try:
            asset_list = subprocess.check_output(
                ['earthengine', 'ls', ini['EXPORT']['output_ws']],
                universal_newlines=True)
            asset_list = [x.strip() for x in asset_list.split('\n') if x]
            # logging.debug(asset_list)
        except ValueError as e:
            logging.info('  Collection doesn\'t exist')
            logging.debug('  {}'.format(str(e)))
            asset_list = []
        except Exception as e:
            logging.error('\n  Unknown error, returning False')
            logging.error(e)
            return False
    elif ini['EXPORT']['export_dest'] == 'CLOUD':
        logging.debug('\nGetting cloud storage file list')
        cloud_list = utils.get_bucket_files(
            ini['EXPORT']['project_name'], ini['EXPORT']['output_ws'])
        # It may be necessary to remove image tile notation
    elif ini['EXPORT']['export_dest'] == 'GDRIVE':
        logging.debug('\nGetting Google drive file list')
        gdrive_list = [
            os.path.join(ini['EXPORT']['output_ws'], x)
            for x in os.listdir(ini['EXPORT']['output_ws'])]
        # Very large tiles may get split up automatically by EE
        # Strip the EE tile notation data from the image list
        gdrive_list = list(set([
            re.sub('-\d{10}-\d{10}.tif', '.tif', x) for x in gdrive_list]))
        # logging.debug(gdrive_list)


    # Get output image based on the study area
    logging.debug('\nBuilding export list')
    export_list = list(image_export_generator(
        ini['INPUTS']['study_area_path'],
        wrs2_coll=ini['INPUTS']['wrs2_coll'],
        cell_size=ini['EXPORT']['cell_size'],
        output_crs=ini['EXPORT']['output_crs'],
        output_osr=ini['EXPORT']['output_osr'],
        wrs2_tile_list=ini['INPUTS']['wrs2_tiles'],
        wrs2_tile_field=ini['INPUTS']['wrs2_tile_field'],
        snap_x=ini['EXPORT']['snap_x'],
        snap_y=ini['EXPORT']['snap_y'],
        wrs2_buffer=ini['INPUTS']['wrs2_buffer']))
    if not export_list:
        logging.error('\nEmpty export list, exiting')
        return False

    # Save export list to json
    with open('export_image.json', 'w') as json_f:
        json.dump(export_list, json_f)


    # Process each image separately (there should only be one)
    logging.info('\nImage Exports')
    for export_i, export_info in enumerate(export_list):
        logging.info('Image: 1')
        logging.debug('  Shape:      {}'.format(export_info['shape']))
        logging.debug('  Transform:  {}'.format(export_info['geo']))
        logging.debug('  Extent:     {}'.format(export_info['extent']))
        logging.debug('  MaxPixels:  {}'.format(export_info['maxpixels']))
        logging.debug('  WRS2 tiles: {}'.format(
            ', '.join(export_info['wrs2_tiles'])))


        if ini['INPUTS']['et_model'] == 'EEFLUX':
            # Get the Landsat collection
            landsat_coll = landsat.get_landsat_coll(
                wrs2_tile_list=export_info['wrs2_tiles'],
                cloud_cover=ini['INPUTS']['cloud_cover'],
                start_date=ini['INTERPOLATE']['start_date'],
                end_date=ini['INTERPOLATE']['end_date'],
                landsat5_flag=ini['INPUTS']['landsat5_flag'],
                landsat7_flag=ini['INPUTS']['landsat7_flag'],
                landsat8_flag=ini['INPUTS']['landsat8_flag'],
                landsat_type='RAD')

            # Compute ETf for each Landsat scene
            # The 'BQA' band is also being returned by the etrf method
            def apply_et_fraction(image):
                etrf_obj = eeflux.EEFlux(ee.Image(image)).etrf
                etrf_img = ee.Image(etrf_obj.select(['etrf'], ['etf'])) \
                    .clamp(-1, 2)
                cloud_mask = landsat.landsat_bqa_cloud_mask_func(
                    ee.Image(etrf_obj. select(['BQA'])))
                return etrf_img.updateMask(cloud_mask) \
                    .copyProperties(image, ['system:time_start'])
            scene_et_fraction_coll = ee.ImageCollection(
                landsat_coll.map(apply_et_fraction))

        else:
            logging.error('\nInvalid/unsupported ET Model: {}'.format(
                ini['INPUTS']['et_model']))
            return False


        # Daily reference ET collection
        # DEADBEEF - Hard coding to GRIDMET for now
        # Should this be retrieved from the model?
        daily_et_reference_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
            .filterDate(ini['INPUTS']['start_date'], ini['INPUTS']['end_date']) \
            .select(['etr'], ['et_reference'])

        # Compute composite/mosaic images for each image date
        daily_et_fraction_coll = ee.ImageCollection(interpolate.aggregate_daily(
            image_coll=scene_et_fraction_coll,
            start_date=ini['INTERPOLATE']['start_date'],
            end_date=ini['INTERPOLATE']['end_date']))

        # Interpolate daily ETf, multiply by daily ETr, and sum to ET
        daily_et_actual_coll = ee.ImageCollection(interpolate.interp_et_coll(
            et_reference_coll=daily_et_reference_coll,
            et_fraction_coll=daily_et_fraction_coll,
            interp_days=ini['INTERPOLATE']['interp_days'],
            interp_type=ini['INTERPOLATE']['interp_type']))

        # Export products
        for product in ini['EXPORT']['products']:
            logging.debug('\n  Product:   {}'.format(product))
            export_id = ini['EXPORT']['export_id_fmt'] \
                .replace('_{index}', '') \
                .format(
                    model=ini['INPUTS']['et_model'].lower(),
                    product=product.lower(),
                    study_area=ini['INPUTS']['study_area_name'],
                    start=ini['INPUTS']['start_date'],
                    end=ini['INPUTS']['end_date'],
                    export=ini['EXPORT']['export_dest'].lower())
            export_id = export_id.replace('-', '')
            logging.debug('    Export ID: {}'.format(export_id))

            if product == 'scene_id':
                # Export the scene list CSV to Google Drive
                if ini['EXPORT']['export_dest'] in ['ASSET', 'GDRIVE']:
                    export_path = os.path.join(
                        ini['EXPORT']['output_ws'], export_id + '.csv')
                elif ini['EXPORT']['export_dest'] == 'CLOUD':
                    export_path = '{}/{}'.format(
                        ini['EXPORT']['output_ws'], export_id + '.csv')
            elif ini['EXPORT']['export_dest'] == 'ASSET':
                export_path = ini['EXPORT']['output_ws'] + '/' + export_id
            elif ini['EXPORT']['export_dest'] == 'CLOUD':
                # Write each product to a separate folder
                export_path = '{}/{}'.format(
                    ini['EXPORT']['output_ws'], export_id + '.tif')
            elif ini['EXPORT']['export_dest'] == 'GDRIVE':
                export_path = os.path.join(
                    ini['EXPORT']['output_ws'], export_id + '.tif')
            logging.debug('    Export folder: {}'.format(
                os.path.dirname(export_path)))
            logging.debug('    Export file: {}'.format(
                os.path.basename(export_path)))

            if overwrite_flag:
                if export_id in tasks.keys():
                    logging.debug('    Task already submitted, cancelling')
                    ee.data.cancelTask(tasks[export_id])

                # This is intentionally not an "elif" so that a task can be
                # cancelled and an existing image/file/asset can be removed
                if (ini['EXPORT']['export_dest'] == 'ASSET' and
                        export_path in asset_list):
                    logging.debug('    Existing asset will be overwritten')
                    # logging.debug('    Asset already exists, removing')
                    # input('Press ENTER to remove asset')
                    # subprocess.call(['earthengine', 'rm', export_path])
                elif (ini['EXPORT']['export_dest'] == 'CLOUD' and
                        export_path in cloud_list):
                    logging.debug('    Export image already exists')
                    # Files in cloud storage are easily overwritten
                    #   so it is unneccesary to manually remove them
                    # Automatically generated image tiles may not get overwritten
                    subprocess.call([
                        'gsutil', 'rm', export_path.replace('.tif', '*.tif')])
                elif (ini['EXPORT']['export_dest'] == 'GDRIVE' and
                        export_path in gdrive_list):
                    logging.debug('    Export image already exists, removing')
                    # Remove automatically generated image tiles
                    for f in glob.glob(export_path.replace('.tif', '*.tif')):
                        os.remove(f)
                    # os.remove(export_path)
            else:
                if export_id in tasks.keys():
                    logging.debug('    Task already submitted, skipping')
                    continue
                elif (ini['EXPORT']['export_dest'] == 'ASSET' and
                        export_path in asset_list):
                    logging.debug('    Asset already exists, skipping')
                    continue
                elif (ini['EXPORT']['export_dest'] == 'CLOUD' and
                        export_path in cloud_list):
                    logging.debug('    Export file already exists, skipping')
                    continue
                elif (ini['EXPORT']['export_dest'] == 'GDRIVE' and
                        os.path.isfile(export_path)):
                    logging.debug('    Export file already exists, skipping')
                    continue

            # Compute target product
            if product == 'scene_id':
                def scene_id_extract(image):
                    return ee.Feature(None).setMulti({
                        'SCENE_ID': ee.String(image.get('SCENE_ID'))})
                scene_id_coll = ee.FeatureCollection(
                    scene_et_fraction_coll.map(scene_id_extract)).sort('SCENE_ID')
            elif product == 'et_actual':
                # Sum daily ET to total ET
                output_image = ee.Image(daily_et_actual_coll.sum())
            elif product == 'et_reference':
                # Sum daily reference ET to total reference ET
                output_image = ee.Image(daily_et_reference_coll.sum())
            elif product == 'et_fraction':
                # Compute mean ETf (ET / ETr)
                output_image = ee.Image(daily_et_actual_coll.sum()) \
                    .divide(ee.Image(daily_et_reference_coll.sum()))
            elif product == 'count':
                # Filter count date range to same period as reference ET
                output_image = ee.Image(daily_et_fraction_coll.filterDate(
                    ini['INPUTS']['start_dt'],
                    ini['INPUTS']['end_dt'] + datetime.timedelta(days=1)).count())
            # elif product == 'count_monthly':
            #     output_image = interpolate.aggregate_monthly(
            #         composite_etf_coll.filterDate(
            #             ini['INPUTS']['start_dt'],
            #             ini['INPUTS']['end_dt'] + datetime.timedelta(days=1)))
            else:
                logging.warning('  Unsupported product type, skipping')
                continue

            # Convert data types for export to Google Drive or Cloud Storage
            if (product in ['et_actual', 'et_reference', 'et_fraction'] and
                    ini['EXPORT']['export_dest'] in ['CLOUD', 'GDRIVE']):
                output_image = output_image.unmask(-9999, False).toFloat()
            # elif (product in ['count', 'count_monthly'] and
            elif (product in ['count'] and
                    ini['EXPORT']['export_dest'] in ['CLOUD', 'GDRIVE']):
                output_image = output_image.unmask(255, False).toUint8()
            elif ini['EXPORT']['export_dest'] in ['ASSET']:
                pass

            # Build export tasks
            if product == 'scene_id':
                if ini['EXPORT']['export_dest'] == 'CLOUD':
                    task = ee.batch.Export.table.toCloudStorage(
                        scene_id_coll,
                        description=export_id,
                        bucket=ini['EXPORT']['bucket_name'],
                        fileNamePrefix=export_id,
                        fileFormat='CSV')
                elif ini['EXPORT']['export_dest'] in ['ASSET', 'GDRIVE']:
                    # Export the scene list CSV to Google Drive
                    task = ee.batch.Export.table.toDrive(
                        scene_id_coll,
                        description=export_id,
                        folder=os.path.basename(ini['EXPORT']['output_ws']),
                        fileNamePrefix=export_id,
                        fileFormat='CSV')
            elif ini['EXPORT']['export_dest'] == 'ASSET':
                # Export the image as an asset
                task = ee.batch.Export.image.toAsset(
                    output_image,
                    description=export_id,
                    assetId=export_path,
                    dimensions=export_info['shape'],
                    crs=export_info['crs'],
                    crsTransform=export_info['geo'],
                    maxPixels=export_info['maxpixels'])
            elif ini['EXPORT']['export_dest'] == 'CLOUD':
                # Export the image to cloud storage
                task = ee.batch.Export.image.toCloudStorage(
                    output_image,
                    description=export_id,
                    bucket=ini['EXPORT']['bucket_name'],
                    fileNamePrefix=export_id,
                    dimensions=export_info['shape'],
                    crs=export_info['crs'],
                    crsTransform=export_info['geo'],
                    # shardSize=,
                    # fileDimensions=,
                    maxPixels=export_info['maxpixels'])
            elif ini['EXPORT']['export_dest'] == 'GDRIVE':
                # Export the images to your Google Drive
                task = ee.batch.Export.image.toDrive(
                    output_image,
                    description=export_id,
                    folder=os.path.basename(ini['EXPORT']['output_ws']),
                    fileNamePrefix=export_id,
                    dimensions=export_info['shape'],
                    crs=export_info['crs'],
                    crsTransform=export_info['geo'],
                    maxPixels=export_info['maxpixels'])
            else:
                logging.debug('  Export task not built, skipping')
                continue

            # Try to start the export task a few times
            logging.debug('  Starting export task')
            for i in range(1, 10):
                try:
                    task.start()
                    break
                except Exception as e:
                    logging.error(
                        '    Error: {}\n    Retrying ({}/10)'.format(e, i))
                    time.sleep(i ** 2)
                    i += 1
            # logging.debug('    Active: {}'.format(task.active()))
            # logging.debug('    Status: {}'.format(task.status()))


def image_export_generator(study_area_path, wrs2_coll,
                           cell_size=30, output_crs=None, output_osr=None,
                           wrs2_tile_list=[], wrs2_tile_field='WRS2_TILE',
                           snap_x=15, snap_y=15, wrs2_buffer=0,
                           n_max=1000, simplify_buffer=1000):
    """Generate image metadata for the study area geometry

    Args:
        study_area_path (str): File path of the study area shapefile
        wrs2_coll (str): WRS2 Landsat footprint asset ID.
            (should default to "projects/eeflux/wrs2_descending_custom")
        cell_size (float): Cell size [m].  Defaults to 30.
        output_crs (str): Output CRS (for setting 'crs' parameter in EE calls).
            Defaults to None.
        output_osr (osr.SpatialReference): Output coordinate system.
            Defaults to None.
        wrs2_tile_field (str): WRS2 tile field name in the fusion table
            Defaults to 'WRS2_TILE'
        wrs2_tile_list (list): User defined WRS2 tile subset
        snap_x (float): X snap coordinate [m].  Defaults to 0.
        snap_y (float): Y snap coordinate [m].  Defaults to 0.
        wrs2_buffer (float): WRS2 footprint buffer distance [m].
            Defaults to 10000.
        n_max (int): Maximum number of WRS2 tiles to join to feature.
            Defaults to 1000.
        simplify_buffer (float): Study area buffer/simplify distance [m].
            Defaults to 1000.

    Yields:
        dict: export information
    """
    logging.info('\nReading study area shapefile')
    logging.info('  {}'.format(study_area_path))
    study_area_ds = ogr.Open(study_area_path, 0)
    study_area_lyr = study_area_ds.GetLayer()
    study_area_osr = study_area_lyr.GetSpatialRef()
    study_area_proj = study_area_osr.ExportToWkt()
    # study_area_osr = study_area_osr.ExportToProj4()
    # logging.debug('  Projection: {}'.format(study_area_proj))
    # Convert WKT to EE WKT
    # study_area_crs = re.sub(
    #     '\s+', '', ee.Projection(study_area_proj).wkt().getInfo())
    study_area_crs = str(study_area_proj)
    logging.debug('  Study area projection: {}'.format(study_area_crs))

    # Get the dissolved/unioned geometry of the study area
    output_geom = ogr.Geometry(ogr.wkbMultiPolygon)
    # shape_list = []
    for study_area_ftr in study_area_lyr:
        # Union each feature
        output_geom = output_geom.Union(
            study_area_ftr.GetGeometryRef())
    study_area_ds = None

    # Project the study area geometry to the output coordinate system
    output_tx = osr.CoordinateTransformation(study_area_osr, output_osr)
    output_geom.Transform(output_tx)

    # Get the output extent from the projected geometry
    output_extent = list(output_geom.GetEnvelope())
    # OGR extents are swapped from GDAL extents
    output_extent[1], output_extent[2] = output_extent[2], output_extent[1]

    # Adjust extent to the cell size
    adjust_size = 2 * cell_size
    output_extent[0] = math.floor((
        output_extent[0] - snap_x) / adjust_size) * adjust_size + snap_x
    output_extent[1] = math.floor((
        output_extent[1] - snap_y) / adjust_size) * adjust_size + snap_y
    output_extent[2] = math.ceil((
        output_extent[2] - snap_x) / adjust_size) * adjust_size + snap_x
    output_extent[3] = math.ceil((
        output_extent[3] - snap_y) / adjust_size) * adjust_size + snap_y

    output_geo = [
        cell_size, 0, output_extent[0], 0, -cell_size, output_extent[3]]

    output_shape = '{0}x{1}'.format(
        int(abs(output_extent[2] - output_extent[0]) / cell_size),
        int(abs(output_extent[3] - output_extent[1]) / cell_size))

    max_pixels = 2 * reduce(mul, map(int, output_shape.split('x')))

    # Create simplified geometries to reduce pixel count
    # output_hull = output_geom.ConvexHull()

    # Buffer/simplify values are assuming the geometry units are in meters
    # DEADBEEF - Need to check that output_osr has units of meters
    output_simplify = output_geom.Buffer(simplify_buffer) \
        .SimplifyPreserveTopology(simplify_buffer)

    # Build EE geometry object from GeoJSON
    output_ee_geom = ee.Geometry(
        json.loads(output_simplify.ExportToJson()), output_crs, False)

    # Pre-filter the WRS2 descending collection
    #   with the buffered study area geometry
    # Then buffer the WRS2 descending collection
    if wrs2_buffer:
        wrs2_coll = ee.FeatureCollection(wrs2_coll) \
            .filterBounds(output_ee_geom.buffer(wrs2_buffer, 1)) \
            .map(lambda ftr: ftr.buffer(wrs2_buffer, 1))
    else:
        wrs2_coll = ee.FeatureCollection(wrs2_coll) \
            .filterBounds(output_ee_geom)

    #  Join intersecting geometries
    output_coll = ee.Join.saveAll(matchesKey='scenes').apply(
        ee.FeatureCollection([ee.Feature(output_ee_geom)]), wrs2_coll,
        ee.Filter.intersects(leftField='.geo', rightField='.geo', maxError=10))

    # Extract WRS2 tile values from joined collection
    def ftr_property(ftr):
        # Calling ".toList()" allows the map to return the WRS2 tiles as a list
        scenes = ee.FeatureCollection(ee.List(ee.Feature(ftr).get('scenes'))) \
            .toList(n_max).map(lambda pr: ee.Feature(pr).get(wrs2_tile_field))
        return ee.Feature(None, {'wrs2_tiles': scenes})
    output_wrs2_info = ee.FeatureCollection(output_coll)\
        .map(ftr_property).getInfo()

    # Pull the WRS2 tile list for the image
    output_wrs2_tiles = map(str, sorted(
        output_wrs2_info['features'][0]['properties']['wrs2_tiles']))

    # Limit the WRS2 tile list
    if wrs2_tile_list:
        output_wrs2_tiles = sorted(list(
            set(wrs2_tile_list) & set(output_wrs2_tiles)))

    yield {
        'crs': output_crs,
        'extent': output_extent,
        'geo': output_geo,
        # 'geojson': json.loads(output_simplify.ExportToJson()),
        # 'index':
        'maxpixels': max_pixels,
        'wrs2_tiles': output_wrs2_tiles,
        'shape': output_shape
    }


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Export ET/ETrF/ETr/count images',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', type=utils.arg_valid_file,
        help='Input file', metavar='FILE')
    parser.add_argument(
        '-o', '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    # Prompt user to select an INI file if not set at command line
    # if not args.ini:
    #     args.ini = utils.get_ini_path(os.getcwd())
    return args


if __name__ == '__main__':
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{}'.format('#' * 80))
    logging.info('{:<20s} {}'.format(
        'Run Time Stamp:', datetime.datetime.now().isoformat(' ')))
    logging.info('{:<20s} {}'.format('Current Directory:', os.getcwd()))
    logging.info('{:<20s} {}'.format(
        'Script:', os.path.basename(sys.argv[0])))

    main(ini_path=args.ini, overwrite_flag=args.overwrite)
