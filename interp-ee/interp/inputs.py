# import ConfigParser
from builtins import input
import datetime
import logging
import os
import re
import sys

import configparser
# from backports import configparser
import ee

# import openet.utils as utils
import utils


def read(ini_path):
    logging.debug('\nReading Input File')
    # Open config file
    config = configparser.ConfigParser()
    try:
        config.read(ini_path)
    except Exception as e:
        logging.error(
            '\nERROR: Input file could not be read, '
            'is not an input file, or does not exist\n'
            '  ini_path={}\n\nException: {}'.format(ini_path, e))
        sys.exit()

    # Force conversion of unicode to strings
    ini = dict()
    for section in config.keys():
        ini[str(section)] = {}
        for k, v in config[section].items():
            ini[str(section)][str(k)] = v
    return ini


def parse_section(ini, section):
    logging.debug('\nParsing/checking {} section'.format(section))
    if section not in ini.keys():
        logging.error(
            '\nERROR: Input file does not have an {} section'.format(section))
        sys.exit()

    if section == 'INPUTS':
        parse_inputs(ini)
    elif section == 'EXPORT':
        parse_export(ini)
    elif section == 'INTERPOLATE':
        parse_interpolate(ini)
    elif section == 'TILE':
        parse_tile(ini)
    # ET model specific sections
    elif section == 'NDVI':
        parse_ndvi(ini)
    elif section == 'EEFLUX':
        parse_eeflux(ini)


def get_param(ini, section, input_name, output_name, get_type,
              default='MANDATORY'):
    """Get INI parameters by type and set default values

    Args:
        ini (dict): Nested dictionary of INI file keys/values
        section (str): Section name
        input_name (str): Parameter name in INI file
        output_name (str): Parameter name in code
        get_type (): Python type
        default (): Default value to use if parameter was not set.
            Defaults to "MANDATORY".
            "MANDATORY" will cause script to exit if key does not exist.
    """

    try:
        if get_type is bool:
            ini[section][output_name] = (
                ini[section][input_name].lower() == "true")
            # ini[section][output_name] = distutils.util.strtobool(
            #     ini[section][input_name])
            # ini[section][output_name] = ini.getboolean(section, input_name)
            # ini[section][output_name] = ini[section].getboolean(input_name)
        elif get_type is int:
            ini[section][output_name] = int(ini[section][input_name])
        elif get_type is float:
            ini[section][output_name] = float(ini[section][input_name])
        elif get_type is list:
            ini[section][output_name] = str(ini[section][input_name])
        else:
            ini[section][output_name] = str(ini[section][input_name])
            # Convert 'None' (strings) to None
            if ini[section][output_name].lower() == 'none':
                ini[section][output_name] = None
    except (KeyError, configparser.NoOptionError):
        if default == 'MANDATORY':
            logging.error(
                '\nERROR: {} was not set in the INI, exiting\n'.format(
                    input_name))
            sys.exit()
        else:
            ini[section][output_name] = default
            logging.debug('  Setting {} = {}'.format(
                input_name, ini[section][output_name]))
    except ValueError:
        logging.error('\nERROR: Invalid value for "{}"'.format(
            input_name))
        sys.exit()
    except Exception as e:
        logging.error('\nERROR: Unhandled error\n  {}'.format(e))
        sys.exit()

    # If the parameter is renamed, remove the old name/parameter
    if input_name != output_name:
        del ini[section][input_name]


def parse_inputs(ini, section='INPUTS'):
    """Parse INPUTS section of INI file"""

    # MANDATORY PARAMETERS
    # section, input_name, output_name, get_type
    param_list = [
        ['et_model', 'et_model', str],
        # Move to INTERPOLATION section eventually
        ['start_date', 'start_date', str],
        ['end_date', 'end_date', str],
        ['study_area_path', 'study_area_path', str]
    ]
    for input_name, output_name, get_type in param_list:
        get_param(ini, section, input_name, output_name, get_type)

    ini[section]['et_model'] = ini[section]['et_model'].upper()

    # Check ET Model
    et_model_list = ['NDVI', 'EEFLUX']
    if not utils.is_option(ini[section]['et_model'], et_model_list):
        logging.error(
            '\nERROR: Invalid ET Model type: {}\n  Must be {}'.format(
                ini[section]['et_model'], ', '.join(et_model_list)))
        sys.exit()

    # Build and check file paths
    if not os.path.isfile(ini[section]['study_area_path']):
        logging.error(
            '\nERROR: The study area shapefile does not exist, '
            'exiting\n  {}'.format(ini[section]['study_area_path']))
        sys.exit()
    ini[section]['study_area_name'] = os.path.splitext(
        os.path.basename(ini[section]['study_area_path']))[0]

    # Check dates
    if not utils.is_date(ini[section]['start_date']):
        logging.error('\nERROR: Invalid start date')
        sys.exit()
    elif not utils.is_date(ini[section]['end_date']):
        logging.error('\nERROR: Invalid end date')
        sys.exit()
    ini[section]['start_dt'] = datetime.datetime.strptime(
        ini['INPUTS']['start_date'], '%Y-%m-%d')
    ini[section]['end_dt'] = datetime.datetime.strptime(
        ini['INPUTS']['end_date'], '%Y-%m-%d')

    # OPTIONAL PARAMETERS
    # param_section, input_name, output_name, get_type, default
    param_list = [
        # Control which Landsat images are used
        ['landsat5_flag', 'landsat5_flag', bool, True],
        ['landsat7_flag', 'landsat7_flag', bool, True],
        ['landsat8_flag', 'landsat8_flag', bool, True],
        # Average cloud cover percentage from metadata
        ['cloud_cover', 'cloud_cover', int, 100],
        # Path/row filtering
        ['wrs2_tiles', 'wrs2_tiles', list, []],
        ['wrs2_tile_fmt', 'wrs2_tile_fmt', str, 'p{:03d}r{:03d}'],
        ['wrs2_coll', 'wrs2_coll', str, 'projects/eeflux/wrs2_descending_custom'],
        ['wrs2_buffer', 'wrs2_buffer', int, 0],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    if ini[section]['cloud_cover'] < 0 or ini[section]['cloud_cover'] > 100:
        logging.error('\nERROR: cloud_cover must be in the range 0-100\n')
        sys.exit()

    # I'm not sure this needs to be checked every time
    try:
        ee.FeatureCollection(ini[section]['wrs2_coll']).first().getInfo()
    except Exception as e:
        logging.error(
            '\nUnhandled exception attempting to read the WRS2 collection\n  '
            'wrs2_coll: {}\n{}\n'.format(ini[section]['wrs2_coll'], e))
        sys.exit()
    if ini[section]['wrs2_buffer'] < 0:
        logging.error(
            '\nERROR: wrs2_buffer must be greater than or equal 0\n')
        sys.exit()

    # Convert WRS2 tile ranges to list
    if ini[section]['wrs2_tiles']:
        ini[section]['wrs2_tiles'] = sorted([
            wrs2_tile.strip()
            for wrs2_tile in ini[section]['wrs2_tiles'].split(',')])

        # Filter user WRS2 tiles based on pre-defined WRS2 tile format
        wrs2_tile_re = re.compile('p(\d{1,3})r(\d{1,3})')
        # logging.debug('\n  {:16s} {}'.format(
        #     'Input WRS2 tiles:', ', '.join(ini['INPUTS']['wrs2_tiles'])))
        ini[section]['wrs2_tiles'] = [
            ini['INPUTS']['wrs2_tile_fmt'].format(
                *map(int, wrs2_tile_re.findall(wrs2_tile)[0]))
            for wrs2_tile in ini[section]['wrs2_tiles']]
        logging.debug('\n  {:16s} {}'.format(
            'WRS2 tiles:', ', '.join(ini[section]['wrs2_tiles'])))

    # Inputs not yet set in INI file
    ini[section]['wrs2_tile_field'] = 'WRS2_TILE'

    # Read WRS2 tiles from file
    # Eventually make this a parameter?
    # pr_input_file = os.path.join('tcorr_conus', 'wrs2_tile_conus.txt')
    # pr_input_file = os.path.join(
    #     'tcorr_western_us', 'wrs2_tile_western_us.txt')
    # pr_input_file = os.path.join(
    #     'tcorr_red_river_basin', 'wrs2_tile_red_river_basin.txt')
    # pr_input_path = os.path.join(output_ws, pr_input_file)
    # logging.debug('\nReading WRS2 tile list from file:\n  {}'.format(
    #     pr_input_path))
    # with open(pr_input_path) as input_f:
    #     user_wrs2_tiles = [x.strip() for x in input_f.readlines()]
    # logging.debug('  {} WRS2 tiles'.format(len(wrs2_tiles)))

    logging.info('\n  {:16s} {}'.format(
        'Start Date:', ini[section]['start_date']))
    logging.info('  {:16s} {}'.format(
        'End Date:', ini[section]['end_date']))


def parse_export(ini, section='EXPORT'):
    """Parse EXPORT section of INI file

    Eventually add check on keywords in export_id_fmt
    """
    # DEADBEEF - This should be at the top but was moved to help debug an
    # error with Fiona
    from osgeo import osr

    # MANDATORY PARAMETERS
    # section, input_name, output_name, get_type
    param_list = [
        ['export_dest', 'export_dest', str],
        ['export_id_fmt', 'export_id_fmt', str],
    ]
    for input_name, output_name, get_type in param_list:
        get_param(ini, section, input_name, output_name, get_type)

    export_dest_list = ['GDRIVE', 'CLOUD', 'ASSET']
    if not utils.is_option(ini[section]['export_dest'], export_dest_list):
        logging.error(
            '\nERROR: Invalid Export Destination: {}\n  Must be {}'.format(
                ini[section]['export_dest'], ', '.join(export_dest_list)))
        sys.exit()

    # DEADBEEF - This might be better in an export module or separate function
    # Export destination specific options
    if ini[section]['export_dest'] == 'GDRIVE':
        logging.info('  Google Drive Export')
        get_param(ini, section, 'export_path', 'output_ws', str)
        ini[section]['output_ws'] = ini[section]['output_ws'].rstrip(os.sep)
        logging.debug('  {:16s} {}'.format(
            'GDrive Workspace:', ini[section]['output_ws']))
        if not os.path.isdir(ini[section]['output_ws']):
            os.makedirs(ini[section]['output_ws'])

    elif ini[section]['export_dest'] == 'ASSET':
        logging.info('  Asset Export')
        get_param(ini, section, 'export_path', 'output_ws', str)
        logging.debug('  {:16s} {}'.format(
            'Asset Worskpace:', ini[section]['output_ws']))

    elif ini[section]['export_dest'] == 'CLOUD':
        logging.info('  Cloud Storage')
        get_param(ini, section, 'project_name', 'project_name', str, 'steel-melody-531')
        get_param(ini, section, 'bucket_name', 'bucket_name', str, None)
        get_param(ini, section, 'bucket_folder', 'bucket_folder', str, None)

        if not ini[section]['project_name']:
            logging.error('\nERROR: {} must be set in INI, exiting\n'.format(
                ini[section]['project_name']))
        # DEADBEEF
        elif ini[section]['project_name'] not in ['eeflux', 'steel-melody-531']:
            logging.error(
                '\nERROR: When exporting to Cloud Storage, the '
                'project_name parameter sets the project name.'
                '  This parameter must be set to "steel-melody-531" for now')
            sys.exit()
        if ini[section]['bucket_name']:
            ini[section]['output_ws'] = 'gs://{}'.format(
                ini[section]['bucket_name'])
        else:
            logging.error('\nERROR: {} must be set in INI, exiting\n'.format(
                ini[section]['bucket_name']))
            sys.exit()
        if ini[section]['bucket_folder']:
            # Add folder to the path if it was set
            ini[section]['output_ws'] = '{}/{}'.format(
                ini[section]['output_ws'], ini[section]['bucket_folder'])
        logging.debug('  {:16s} {}'.format(
            'Project:', ini[section]['project_name']))
        logging.debug('  {:16s} {}'.format(
            'Bucket:', ini[section]['bucket_name']))
        logging.debug('  {:16s} {}'.format(
            'Output Workspace:', ini[section]['output_ws']))

        # DEADBEEF - inputs.py doesn't have a function get_buckets
        # bucket_list = utils.get_buckets(ini[section]['project_name'])
        # if ini[section]['bucket_name'] not in bucket_list:
        #     logging.error(
        #         '\nERROR: The bucket "{}" does not exist, exiting'.format(
        #             ini[section]['bucket_name']))
        #     return False
        #     # Try creating the storage bucket if it doesn't exist using gsutil
        #     # For now, I think it is better to make the user go do this
        #     # subprocess.check_output([
        #     #     'gsutil', 'mb', '-p', ini[section]['project_name'],
        #     #     'gs://{}-{}'.format(
        #     #         ini[section]['project_name'],
        #     #         ini[section]['bucket_name'])])

    # OPTIONAL PARAMETERS
    # param_section, input_name, output_name, get_type, default
    param_list = [
        #
        ['cell_size', 'cell_size', int, 30],
        ['snap_x', 'snap_x', int, 0],
        ['snap_y', 'snap_y', int, 0],
        #
        ['output_epsg', 'output_epsg', str, None],
        ['output_proj4', 'output_proj4', str, None],
        ['output_wkt', 'output_wkt', str, None],
        # Control output products
        ['et_actual_flag', 'et_actual_flag', bool, True],
        ['et_reference_flag', 'et_reference_flag', bool, True],
        ['et_fraction_flag', 'et_fraction_flag', bool, True],
        ['count_flag', 'count_flag', bool, False],
        ['scene_id_flag', 'scene_id_flag', bool, False],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    # Check optional inputs
    if ini[section]['cell_size'] < 0:
        logging.error('\nERROR: cell_size must be greater than 0\n')
        sys.exit()
    elif ini[section]['snap_x'] < 0 or ini[section]['snap_y'] < 0:
        logging.error(
            '\nERROR: snap_x/snap_y must be greater than 0\n')
        sys.exit()

    # Determine output coordinate system
    ini[section]['output_osr'] = osr.SpatialReference()
    if (ini[section]['output_epsg'] and
            utils.is_number(ini[section]['output_epsg'])):
        ini[section]['output_osr'].ImportFromEPSG(
            int(ini[section]['output_epsg']))
        ini[section]['output_crs'] = 'EPSG:' + ini[section]['output_epsg']
    elif (ini[section]['output_epsg'] and
            ini[section]['output_epsg'].startswith('EPSG:')):
        ini[section]['output_osr'].ImportFromEPSG(
            int(ini[section]['output_epsg'].replace('EPSG:', '')))
        ini[section]['output_crs'] = ini[section]['output_epsg']
    elif ini[section]['output_proj4']:
        ini[section]['output_osr'].ImportFromProj4(
            ini[section]['output_proj4'])
        ini[section]['output_crs'] = str(
            ini[section]['output_osr'].ExportToWkt())
    elif ini[section]['output_wkt']:
        ini[section]['output_osr'].ImportFromWkt(
            ini[section]['output_wkt'])
        ini[section]['output_crs'] = str(
            ini[section]['output_osr'].ExportToWkt())
    else:
        ini[section]['output_osr'] = None
        ini[section]['output_crs'] = None
    if ini[section]['output_crs']:
        logging.debug('\n  {:16s} {}'.format(
            'Output crs:', ini[section]['output_crs']))

    # Only keep products that are going to be computed
    product_dict = {
        'et_actual': ini[section]['et_actual_flag'],
        'et_fraction': ini[section]['et_fraction_flag'],
        'et_reference': ini[section]['et_reference_flag'],
        'count': ini[section]['count_flag'],
        # 'count_monthly': ini[section]['count_monthly_flag'],
        # 'mask': ini[section]['mask_flag'],
        'scene_id': ini[section]['scene_id_flag']
    }
    ini[section]['products'] = [k for k, v in product_dict.items() if v]

    logging.debug('  {:16s} {}'.format('ETa:', ini[section]['et_actual_flag']))
    logging.debug('  {:16s} {}'.format('ETf:', ini[section]['et_fraction_flag']))
    logging.debug('  {:16s} {}'.format('ETr:', ini[section]['et_reference_flag']))
    logging.debug('  {:16s} {}'.format('Count:', ini[section]['count_flag']))
    # logging.debug('  {:16s} {}'.format(
    #     'Monthly Count:', ini[section]['count_monthly_flag']))
    # logging.debug('  {:16s} {}'.format('Mask:', ini[section]['mask_flag']))
    logging.debug('  {:16s} {}'.format(
        'Scene ID CSV:', ini[section]['scene_id_flag']))


def parse_tile(ini, section='TILE'):
    """Parse TILE section of INI file"""

    # OPTIONAL PARAMETERS
    # section, input_name, output_name, get_type, default
    param_list = [
        ['tile_cells', 'tile_cells', int, 2000],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    if ini[section]['tile_cells'] < 1:
        logging.error(
            '\nERROR: Output tile size must be a positive integer\n')
        sys.exit()


def parse_interpolate(ini, section='INTERPOLATE'):
    """Parse INTERPOLATE section of INI file"""

    # OPTIONAL PARAMETERS
    # section, input_name, output_name, get_type, default
    param_list = [
        ['interp_days', 'interp_days', int, 64],
        ['interp_type', 'interp_type', str, 'LINEAR'],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    if ini[section]['interp_days'] < 1 or ini[section]['interp_days'] > 128:
        logging.error(
            '\nERROR: interp_days must be in the range 1-128\n')
        sys.exit()

    interp_type_list = ['LINEAR']
    ini[section]['interp_type'] = ini[section]['interp_type'].upper()
    if not utils.is_option(ini[section]['interp_type'], interp_type_list):
        logging.error(
            '\nERROR: Invalid interpolation type: {}\n  Must be {}'.format(
                ini[section]['interp_type'], ', '.join(interp_type_list)))
        sys.exit()

    # This assumes start/end date were already parsed
    # Should sections be able to interact with values from other sections?
    ini[section]['start_dt'] = (
        ini['INPUTS']['start_dt'] -
        datetime.timedelta(days=ini[section]['interp_days']))
    ini[section]['end_dt'] = (
        ini['INPUTS']['end_dt'] +
        datetime.timedelta(days=ini[section]['interp_days']))
    ini[section]['start_date'] = ini[section]['start_dt'].date().isoformat()
    ini[section]['end_date'] = ini[section]['end_dt'].date().isoformat()

    logging.info('\n  {:16s} {}'.format(
        'Interp. Start:', ini[section]['start_date']))
    logging.info('  {:16s} {}'.format(
        'Interp. End:', ini[section]['end_date']))
    logging.info('  {:16s} {}'.format(
        'Interp. Days:', ini[section]['interp_days']))


def parse_ndvi(ini, section='NDVI'):
    """"""

    # OPTIONAL PARAMETERS
    # section, input_name, output_name, description, get_type, default
    param_list = [
        ['m', 'm', float, 1.25],
        ['b', 'b', float, 0.0],
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)


def parse_eeflux(ini, section='EEFLUX'):
    """Parse EEFLUX section of INI file"""

    # MANDATORY PARAMETERS
    # section, input_name, output_name, get_type
    # param_list = [
    #     ['eto_source', 'eto_source', str],
    # ]
    # for input_name, output_name, get_type in param_list:
    #     get_param(ini, section, input_name, output_name, get_type)

    # OPTIONAL PARAMETERS
    # section, input_name, output_name, get_type, default
    param_list = [
        ['refet_source', 'refet_source', str, 'GRIDMET'],
        ['refet_type', 'refet_type', str, 'ETR']
        # ['refet_factor', 'refet_factor', float, 1.0]
    ]
    for input_name, output_name, get_type, default in param_list:
        get_param(ini, section, input_name, output_name, get_type, default)

    ini[section]['refet_source'] = ini[section]['refet_source'].upper()
    ini[section]['refet_type'] = ini[section]['refet_type'].upper()
    ini[section]['refet_factor'] = 1.0

    refet_source_list = ['GRIDMET']
    if (ini[section]['refet_source'] and
            not utils.is_option(ini[section]['refet_source'], refet_source_list)):
        logging.error(
            '\nERROR: Invalid Reference ET Source type: {}\n  Must be {}'.format(
                ini[section]['refet_source'], ', '.join(refet_source_list)))
        sys.exit()

    refet_type_list = ['ETO', 'ETR']
    if (ini[section]['refet_type'] and
            not utils.is_option(ini[section]['refet_type'], refet_type_list)):
        logging.error(
            '\nERROR: Invalid Reference ET type: {}\n  Must be {}'.format(
                ini[section]['refet_type'], ', '.join(refet_type_list)))
        sys.exit()
