#--------------------------------
# Name:         interpolate_tile_export.py
# Purpose:      Export interpolated ETa/ETf/ETo/count image tiles
#--------------------------------

import argparse
from builtins import input
import datetime
from functools import reduce
import json
import logging
import math
from operator import mul
import os
import pprint
import sys
import time

import ee
from osgeo import ogr, osr

# Import interpolation functions
# This is an awful way of getting the parent folder into the path
# We really should package this up as a module with a setup.py
open_et_test_path = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(__file__)))))
sys.path.insert(0, os.path.join(open_et_test_path, 'interp-ee', 'interp'))
import inputs
import interpolate
import landsat
import utils

# Import ET Models
sys.path.insert(0, os.path.join(open_et_test_path, 'eeflux-ee'))
import eeflux


def main(ini_path=None, overwrite_flag=False,
         tile_i='', tile_j=''):
    """Export annual ET/ETrF/ETr/count image tiles

    Args:
        ini_path (str): Input file path
        overwrite_flag (bool): if True, overwrite existing files
        tile_i (str): comma separated list and/or range of tile row indices
        tile_j (str): comma separated list and/or range of tile columns indices

    Returns:
        None
    """
    logging.info('\nExport annual ET/ETrF/ETr/count image tiles')

    # Read config file
    ini = inputs.read(ini_path)
    inputs.parse_section(ini, section='INPUTS')
    inputs.parse_section(ini, section='INTERPOLATE')
    inputs.parse_section(ini, section='EXPORT')
    inputs.parse_section(ini, section='TILE')
    inputs.parse_section(ini, section=ini['INPUTS']['et_model'])

    if ini['EXPORT']['export_dest'] == 'ASSET':
        logging.error('\nERROR: ASSET tile export is not supported')
        return False

    # Limit tile ranges from command line
    # Eventually move to config file?
    try:
        tile_i_list = list(utils.parse_int_set(tile_i))
    except:
        tile_i_list = []
    try:
        tile_j_list = list(utils.parse_int_set(tile_j))
    except:
        tile_j_list = []

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
    if ini['EXPORT']['export_dest'] == 'CLOUD':
        logging.debug('\nGetting cloud storage file list')
        cloud_list = utils.get_bucket_files(
            ini['EXPORT']['project_name'], ini['EXPORT']['output_ws'])
        # It may be necessary to remove image tile notation
    elif ini['EXPORT']['export_dest'] == 'GDRIVE':
        logging.debug('\nGetting Google drive file list')
        gdrive_list = [
            os.path.join(ini['EXPORT']['output_ws'], x)
            for x in os.listdir(ini['EXPORT']['output_ws'])]
        # It may be necessary to remove image tile notation
        # Very large tiles may get split up automatically by EE
        # Strip the EE tile notation data from the image list
        # gdrive_list = list(set([
        #     re.sub('-\d{10}-\d{10}.tif', '.tif', x)
        #     for x in os.listdir(ini['EXPORT']['output_ws'])]))
        # logging.debug(gdrive_list)


    # Get list of tiles that intersect the study area
    logging.debug('\nBuilding export list')
    export_list = list(tile_export_generator(
        ini['INPUTS']['study_area_path'],
        wrs2_coll=ini['INPUTS']['wrs2_coll'],
        cell_size=ini['EXPORT']['cell_size'],
        output_crs=ini['EXPORT']['output_crs'],
        output_osr=ini['EXPORT']['output_osr'],
        wrs2_tile_list=ini['INPUTS']['wrs2_tiles'],
        wrs2_tile_field=ini['INPUTS']['wrs2_tile_field'],
        snap_x=ini['EXPORT']['snap_x'],
        snap_y=ini['EXPORT']['snap_y'],
        tile_cells=ini['TILE']['tile_cells'],
        wrs2_buffer=ini['INPUTS']['wrs2_buffer']))
    if not export_list:
        logging.error('\nEmpty export list, exiting')
        return False

    # Save export list to json
    with open('export_tiles.json', 'w') as json_f:
        json.dump(export_list, json_f)


    # Process each tile separately
    logging.info('\nImage Exports')
    for export_n, export_info in enumerate(export_list):
        logging.info('Tile: {}  ({}/{})'.format(
            export_info['index'], export_n + 1, len(export_list)))

        tile_i, tile_j = map(int, export_info['index'].split('_'))
        if tile_i_list and int(tile_i) not in tile_i_list:
            logging.debug('  Skipping tile')
            continue
        elif tile_j_list and int(tile_j) not in tile_j_list:
            logging.debug('  Skipping tile')
            continue
        # # DEADBEEF - Only test 1 block
        # if [int(tile_i), int(tile_j)] != [3, 3]:
        #     continue

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
            export_id = ini['EXPORT']['export_id_fmt'].format(
                model=ini['INPUTS']['et_model'].lower(),
                product=product.lower(),
                study_area=ini['INPUTS']['study_area_name'],
                index=export_info['index'],
                start=ini['INPUTS']['start_date'],
                end=ini['INPUTS']['end_date'],
                export=ini['EXPORT']['export_dest'].lower())
            export_id = export_id.replace('-', '')
            logging.debug('    Export ID: {}'.format(export_id))

            if product == 'scene_id':
                # Export the scene list CSV to Google Drive
                if ini['EXPORT']['export_dest'] == 'GDRIVE':
                    export_path = os.path.join(
                        ini['EXPORT']['output_ws'], export_id + '.csv')
                elif ini['EXPORT']['export_dest'] == 'CLOUD':
                    export_path = '{}/{}/{}'.format(
                        ini['EXPORT']['output_ws'], product, export_id + '.csv')
            elif ini['EXPORT']['export_dest'] == 'CLOUD':
                # Write each product to a separate folder
                export_path = '{}/{}/{}'.format(
                    ini['EXPORT']['output_ws'], product, export_id + '.tif')
            elif ini['EXPORT']['export_dest'] == 'GDRIVE':
                export_path = os.path.join(
                    ini['EXPORT']['output_ws'], export_id + '.tif')
            else:
                logging.warning('  Unsupported product type, skipping')
                continue
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
                if (ini['EXPORT']['export_dest'] == 'CLOUD' and
                        export_path in cloud_list):
                    logging.debug('    Export image already exists')
                    # Files in cloud storage are easily overwritten
                    #   so it is unneccesary to manually remove them
                    # # This would remove an existing file
                    # subprocess.call(['gsutil', 'rm', export_path])
                elif (ini['EXPORT']['export_dest'] == 'GDRIVE' and
                        export_path in gdrive_list):
                    logging.debug('    Export image already exists, removing')
                    os.remove(export_path)
                    # Remove automatically generated image tiles
                    # for f in glob.glob(export_path.replace('.tif', '*.tif')):
                    #     os.remove(f)
            else:
                if export_id in tasks.keys():
                    logging.debug('    Task already submitted, skipping')
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
                        fileNamePrefix='{}/{}'.format(product, export_id),
                        fileFormat='CSV')
                elif ini['EXPORT']['export_dest'] == 'GDRIVE':
                    # Export the scene list CSV to Google Drive
                    task = ee.batch.Export.table.toDrive(
                        scene_id_coll,
                        description=export_id,
                        folder=os.path.basename(ini['EXPORT']['output_ws']),
                        fileNamePrefix=export_id,
                        fileFormat='CSV')
            elif ini['EXPORT']['export_dest'] == 'CLOUD':
                # Export the image to cloud storage
                task = ee.batch.Export.image.toCloudStorage(
                    output_image,
                    description=export_id,
                    bucket=ini['EXPORT']['bucket_name'],
                    fileNamePrefix='{}/{}'.format(product, export_id),
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


def tile_export_generator(study_area_path, wrs2_coll,
                          cell_size=30, output_crs=None, output_osr=None,
                          wrs2_tile_list=[], wrs2_tile_field='WRS2_TILE',
                          snap_x=15, snap_y=15, tile_cells=2000,
                          wrs2_buffer=0, n_max=1000, simplify_buffer=240):
    """Generate tiles and metadata for the study area geometry

    This function is only a generator in order to match the image and WRS2 tile
        export functions.

    Args:
        study_area_path (str): File path of the study area shapefile
        wrs2_coll (str): WRS2 Landsat footprint asset ID.
            (should default to "projects/ssebop-gee/wrs2_descending_custom")
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
        tile_cells (int): Tile width and height [pixels]. Defaults to 2000.
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
    # tile_info['crs'] = re.sub(
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
    logging.debug('  {:16s} {}'.format('Output Extent:', output_extent))

    # Compute tile size (in meters)
    tile_size = float(tile_cells) * cell_size

    # Expand extent to fully include tiles
    output_extent[0] = math.floor(
        (output_extent[0] - snap_x) / tile_size) * tile_size + snap_x
    output_extent[1] = math.floor(
        (output_extent[1] - snap_y) / tile_size) * tile_size + snap_y
    output_extent[2] = math.ceil(
        (output_extent[2] - snap_x) / tile_size) * tile_size + snap_x
    output_extent[3] = math.ceil(
        (output_extent[3] - snap_y) / tile_size) * tile_size + snap_y
    logging.debug('  {:16s} {}'.format('Adjusted Extent:', output_extent))

    # Create simplified geometries to speed up checking tile intersections
    output_hull = output_geom.ConvexHull()

    # Buffer/simplify values are assuming the geometry units are in meters
    # DEADBEEF - Need to check that output_osr has units of meters
    output_simplify = output_geom.Buffer(simplify_buffer) \
        .SimplifyPreserveTopology(simplify_buffer)

    # Generate an EE feature
    output_ee_geom = ee.Geometry(
        json.loads(output_simplify.ExportToJson()), output_crs, False)

    # # DEBUG code for checking geometries
    # def point_count(geom):
    #     count = 0
    #     if geom.GetGeometryName() == 'POLYGON':
    #         for ring in geom:
    #             count += ring.GetPointCount()
    #     elif geom.GetGeometryName() == 'MULTIPOLYGON':
    #         for poly in geom:
    #             count += point_count(poly)
    #     return count
    # print(point_count(output_geom))
    # print(point_count(output_hull))
    # print(point_count(study_area_simplify))

    # # shape_list = [shape(json.loads(output_hull.ExportToJson()))]
    # shape_list = [shape(json.loads(study_area_simplify.ExportToJson()))]

    # # Draw study area geometries
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.set_xlim(output_extent[0], output_extent[2])
    # ax.set_ylim(output_extent[1], output_extent[3])
    # patches = [
    #     PolygonPatch(p, fc='none', ec='0.1', zorder=1)
    #     for p in MultiPolygon(shape_list)]
    # ax.add_collection(PatchCollection(patches, match_original=True))
    # ax.set_aspect(1.0)
    # plt.grid()
    # plt.show()
    # plt.close()
    # END DEBUG code

    logging.info('\nBuilding tiles')
    tile_cols = int(round(abs((
        output_extent[0] - output_extent[2]) / tile_size), 0))
    tile_rows = int(round(abs((
        output_extent[3] - output_extent[1]) / tile_size), 0))
    output_shape = '{0}x{0}'.format(tile_cells)

    max_pixels = 2 * reduce(mul, map(int, output_shape.split('x')))

    logging.debug('  {:16s} {} {}'.format('Cols/Rows:', tile_cols, tile_rows))
    logging.debug('  {:16s} {}'.format('Shape:', output_shape))
    logging.debug('  {:16s} {}'.format('Max Pixels:', max_pixels))

    # Build a list of tiles object that a function might return
    # We need a list of tile indices (or geometries) in order to build
    #   the tile collection below
    export_list = []
    for tile_i in range(0, tile_rows):
        for tile_j in range(0, tile_cols):
            tile_geo = [
                cell_size, 0, output_extent[0] + tile_j * tile_size,
                0, -cell_size, output_extent[3] - tile_i * tile_size]
            tile_extent = [
                tile_geo[2], tile_geo[5] - tile_size,
                tile_geo[2] + tile_size, tile_geo[5]]
            tile_geom = extent_geom(tile_extent)
            if not tile_geom.Intersects(output_hull):
                continue
            elif not tile_geom.Intersects(output_simplify):
                continue
            # elif not tile_geom.Intersects(study_area_geom):
            #     continue
            export_list.append({
                'crs': output_crs,
                'extent': tile_extent,
                'geo': tile_geo,
                'index': '{}_{}'.format(tile_i, tile_j),
                'maxpixels': max_pixels,
                'shape': output_shape,
            })

    # Build feature collection of tiles
    # Assume all tiles have the same CRS
    def tile_features(index):
        row = ee.Number.parse(ee.String(index).split('_').get(0))
        col = ee.Number.parse(ee.String(index).split('_').get(1))
        xmin = ee.Number(output_extent[0]).add(col.multiply(tile_size))
        ymax = ee.Number(output_extent[3]).subtract(row.multiply(tile_size))
        return ee.Feature(
            ee.Geometry.Rectangle(
                [xmin, ymax.subtract(tile_size), xmin.add(tile_size), ymax],
                ee.Projection(output_crs), False),
            {'index': row.format('%d').cat('_').cat(col.format('%d'))})
    index_list = [export_info['index'] for export_info in export_list]
    # If Internal server errors lead here, try making the tile_size smaller or
    #   cell size smaller.  It seems like there might be a problem when the
    #   rectangle is too big for the coordinate system?
    tiles_coll = ee.FeatureCollection(ee.List(index_list).map(tile_features))

    # Pre-filter the WRS2 descending collection
    #   with the buffered tile geometry
    # Then buffer the WRS2 descending collection
    if wrs2_buffer:
        wrs2_coll = ee.FeatureCollection(wrs2_coll) \
            .filterBounds(output_ee_geom.buffer(wrs2_buffer, 1)) \
            .map(lambda ftr: ftr.buffer(wrs2_buffer, 1))
    else:
        wrs2_coll = ee.FeatureCollection(wrs2_coll) \
            .filterBounds(output_ee_geom)

    #  Join intersecting geometries
    tiles_coll = ee.Join.saveAll(matchesKey='scenes').apply(
        tiles_coll, wrs2_coll,
        ee.Filter.intersects(leftField='.geo', rightField='.geo', maxError=10))

    def tile_scenes(tile):
        # Calling ".toList()" allows the map to return the WRS2 tiles as a list
        scenes = ee.FeatureCollection(ee.List(ee.Feature(tile).get('scenes'))) \
            .toList(n_max).map(lambda ftr: ee.Feature(ftr).get(wrs2_tile_field))
        return ee.Feature(None, {
            'index': tile.get('index'),
            'wrs2_tiles': scenes})
    tile_wrs2_info = ee.FeatureCollection(tiles_coll.map(tile_scenes)).getInfo()
    tile_wrs2_dict = {
        str(t['properties']['index']): map(str, t['properties']['wrs2_tiles'])
        for t in tile_wrs2_info['features']}

    # Pull the WRS2 tile list for each tile
    # Only yield exports that have intersecting WRS2 tiles
    for export_info in export_list:
        try:
            export_info['wrs2_tiles'] = sorted(
                tile_wrs2_dict[export_info['index']])
            yield export_info
        except KeyError:
            pass
        # t_index = export_info['index']
        # try:
        #     # export_list[i]['wrs2_tiles'] = tile_pr_dict[t_index]
        # except KeyError:
        #     # logging.debug('  Tile {} - no WRS2 tiles'.format(t_index))
        #     # export_list[i]['wrs2_tiles'] = []


def extent_geom(extent):
    """GDAL geometry object of the extent list

    Args:
        extent (list):  GDAL style extent (xmin, ymin, xmax, ymax)

    Returns:
        ogr.geometry
    """
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(extent[0], extent[3])
    ring.AddPoint(extent[2], extent[3])
    ring.AddPoint(extent[2], extent[1])
    ring.AddPoint(extent[0], extent[1])
    ring.CloseRings()
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    return polygon


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Export ET/ETrF/ETr/count image tiles',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-i', '--ini', type=utils.arg_valid_file,
        help='Input file', metavar='FILE')
    parser.add_argument(
        '--rows', default='',
        help='Rows to process (comma separate or range)')
    parser.add_argument(
        '--cols', default='',
        help='Columns to process (comma separate or range)')
    parser.add_argument(
        '-o', '--overwrite', default=False, action='store_true',
        help='Force overwrite of existing files')
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    # Prompt user to select an INI file if not set at command line
    if not args.ini:
        args.ini = utils.get_ini_path(os.getcwd())
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

    main(ini_path=args.ini, overwrite_flag=args.overwrite,
         tile_i=args.rows, tile_j=args.cols)
