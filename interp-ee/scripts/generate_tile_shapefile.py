#--------------------------------
# Name:         generate_tile_shapefile.py
# Purpose:      Generate a shapefile showing the image tile locations
#--------------------------------

import argparse
from builtins import input
import datetime
import logging
import math
import os
import sys

from osgeo import ogr, osr

# Import interpolation functions
# This is an awful way of getting the parent folder into the path
# We really should package this up as a module with a setup.py
open_et_test_path = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(__file__)))))
sys.path.insert(0, os.path.join(open_et_test_path, 'interp-ee', 'interp'))
import inputs
import utils


def main(ini_path=None, overwrite_flag=False, tile_i='', tile_j=''):
    """Export annual ET/ETrF/ETr/count image tiles

    Parameters
    ----------
    ini_path : str
        Input file path.
    overwrite_flag : bool, optional
        If True, overwrite existing files (the default is False).
    tile_i : str
        Comma separated list and/or range of tile row indices.
    tile_j : str
        Comma separated list and/or range of tile columns indices.

    Returns
    -------
    None

    """
    logging.info('\nGenerate tile shapefile')

    # Read config file
    ini = inputs.read(ini_path)
    inputs.parse_section(ini, section='INPUTS')
    inputs.parse_section(ini, section='EXPORT')

    output_path = ini['INPUTS']['study_area_path'].replace('.shp', '_tiles.shp')
    if os.path.isfile(output_path) and not overwrite_flag:
        logging.info(
            '\nOutput shapefile already exists and overwrite_flag is False\n')
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


    # Get list of tiles that intersect the study area
    logging.debug('\nBuilding export list')
    export_list = list(tile_export_generator(
        ini['INPUTS']['study_area_path'],
        cell_size=ini['EXPORT']['cell_size'],
        output_osr=ini['EXPORT']['output_osr'],
        snap_x=ini['EXPORT']['snap_x'],
        snap_y=ini['EXPORT']['snap_y'],
        tile_cells=ini['TILE']['tile_cells']))
    if not export_list:
        logging.error('\nEmpty export list, exiting')
        return False


    # Build the output shapefile
    # Write the scene Tcorr values to the shapefile
    logging.info('\nWriting tiles to the shapefile')
    logging.debug('  {}'.format(output_path))
    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    output_ds = shp_driver.CreateDataSource(output_path)
    output_lyr = output_ds.CreateLayer(
        output_path, ini['EXPORT']['output_osr'], ogr.wkbPolygon)

    field_name = ogr.FieldDefn('COL', ogr.OFTInteger)
    field_name.SetWidth(3)
    output_lyr.CreateField(field_name)
    field_name = ogr.FieldDefn('ROW', ogr.OFTInteger)
    field_name.SetWidth(3)
    output_lyr.CreateField(field_name)

    # Write each tile separately
    for export_n, export_info in enumerate(export_list):
        logging.info('Tile: {}  ({}/{})'.format(
            export_info['index'], export_n + 1, len(export_list)))
        logging.debug('  Extent: {}'.format(export_info['extent']))

        tile_i, tile_j = map(int, export_info['index'].split('_'))
        if tile_i_list and int(tile_i) not in tile_i_list:
            logging.debug('  Skipping tile')
            continue
        elif tile_j_list and int(tile_j) not in tile_j_list:
            logging.debug('  Skipping tile')
            continue

        feature = ogr.Feature(output_lyr.GetLayerDefn())
        feature.SetField('COL', tile_j)
        feature.SetField('ROW', tile_i)
        polygon = ogr.CreateGeometryFromWkt(
            "POLYGON(({0} {1}, {2} {1}, {2} {3}, {0} {3}, {0} {1}))".format(
                *export_info['extent']))
        feature.SetGeometry(polygon)
        output_lyr.CreateFeature(feature)
        feature = None
    output_ds = None


def tile_export_generator(study_area_path, cell_size=30, output_osr=None,
                          snap_x=15, snap_y=15, tile_cells=2000,
                          simplify_buffer=240):
    """Generate tiles and metadata for the study area geometry

    This function is only a generator in order to match the image and WRS2 tile
        export functions.

    Args:
        study_area_path (str): File path of the study area shapefile
        cell_size (float): Cell size [m].  Defaults to 30.
        output_osr (osr.SpatialReference): Output coordinate system.
            Defaults to None.
        snap_x (float): X snap coordinate [m].  Defaults to 0.
        snap_y (float): Y snap coordinate [m].  Defaults to 0.
        tile_cells (int): Tile width and height [pixels]. Defaults to 2000.
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

    logging.info('\nBuilding tiles')
    tile_cols = int(round(abs((
        output_extent[0] - output_extent[2]) / tile_size), 0))
    tile_rows = int(round(abs((
        output_extent[3] - output_extent[1]) / tile_size), 0))

    logging.debug('  {:16s} {} {}'.format('Cols/Rows:', tile_cols, tile_rows))

    # Build a list of tiles object that a function might return
    export_list = []
    for tile_j in range(0, tile_cols):
        for tile_i in range(0, tile_rows):
            tile_geo = [
                cell_size, 0, output_extent[0] + tile_j * tile_size,
                0, -cell_size, output_extent[3] - tile_i * tile_size]
            tile_extent = [
                tile_geo[2], tile_geo[5] - tile_size,
                tile_geo[2] + tile_size, tile_geo[5]]
            yield {
                'index': '{}_{}'.format(tile_i, tile_j),
                'extent': tile_extent}


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='Generate tile shapefile',
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

    main(ini_path=args.ini, overwrite_flag=args.overwrite,
         tile_i=args.rows, tile_j=args.cols)
