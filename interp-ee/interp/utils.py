from builtins import int
import argparse
# import calendar
import datetime
import logging
import os
import subprocess
import sys

import ee

ee.Initialize()


def arg_valid_date(date_str):
    """Argparse specific function for validating date strings"""
    try:
        datetime.datetime.strptime(date_str, "%Y-%m-%d")
        return date_str
    except ValueError:
        raise argparse.ArgumentTypeError(
            '{} is an invalid date'.format(date_str))


def arg_valid_file(file_path):
    """Argparse specific function for testing if file exists

    Convert relative paths to absolute paths
    """
    if os.path.isfile(os.path.abspath(os.path.realpath(file_path))):
        return os.path.abspath(os.path.realpath(file_path))
        # return file_path
    else:
        raise argparse.ArgumentTypeError('{} does not exist'.format(file_path))


def get_ee_assets(asset_id, shell_flag=False):
    """Return Google Earth Engine assets

    Parameters
    ----------
    asset_id : str
        A folder or image collection ID.
    shell_flag : bool
        If True, execute the command through the shell (the default is True).

    Returns
    -------
    list of asset names

    """
    try:
        asset_list = subprocess.check_output(
            ['earthengine', 'ls', asset_id],
            universal_newlines=True, shell=shell_flag)
        asset_list = [x.strip() for x in asset_list.split('\n') if x]
        # logging.debug(asset_list)
    except ValueError as e:
        logging.info('  Collection doesn\'t exist')
        logging.debug('  {}'.format(str(e)))
        asset_list = []
    except Exception as e:
        logging.error('\n  Unknown error, returning False')
        logging.error(e)
        sys.exit()
    return asset_list


def get_bucket_files(project_name, bucket_name, shell_flag=True):
    """Return Google Cloud Storage buckets associated with project

    Parameters
    ----------
    project_name : str
        AppEngine project name.
    bucket_name : str
        Google Storage bucket name.
    shell_flag : bool
        If True, ? (the default is True).

    Returns
    -------
    list of file names

    """
    try:
        file_list = subprocess.check_output(
            ['gsutil', 'ls', '-r', '-p', project_name, bucket_name],
            universal_newlines=True, shell=shell_flag)
        # file_list = [x.strip() for x in file_list.split('\n') if x]
    except Exception as e:
        logging.error(
            '\nERROR: There was a problem getting the bucket file list ' +
            'using gsutil, exiting')
        logging.error('  Exception: {}'.format(e))
        sys.exit()
    return file_list


def get_ee_tasks(states=['RUNNING', 'READY']):
    """Return current active tasks

    Parameters
    ----------
    states : list

    Returns
    -------
    dict of task descriptions (key) and task IDs (value)

    """

    logging.debug('\nActive Tasks')
    try:
        task_list = ee.data.getTaskList()
    except Exception as e:
        logging.warning(
            '  Exception retrieving task list'
            '  {}'.format(e))
        return {}

    task_list = sorted([
        [t['state'], t['description'], t['id']] for t in task_list
        if t['state'] in states])
    if task_list:
        logging.debug('  {:8s} {}'.format('STATE', 'DESCRIPTION'))
        logging.debug('  {:8s} {}'.format('=====', '==========='))
    else:
        logging.debug('  None')

    tasks = {}
    for t_state, t_desc, t_id in task_list:
        logging.debug('  {:8s} {}'.format(t_state, t_desc))
        tasks[t_desc] = t_id
        # tasks[t_id] = t_desc
    return tasks


def is_date(input_date):
    """Check that a date string is ISO format (YYYY-MM-DD)

    DEADBEEF - It would probably make more sense to have this function
      parse the date using dateutil parser (http://labix.org/python-dateutil)
      and return the ISO format string
    """
    try:
        datetime.datetime.strptime(input_date, "%Y-%m-%d")
        return True
    except ValueError:
        return False


def is_number(x):
    try:
        float(x)
        return True
    except:
        return False


def is_option(value, options):
    """Test if a value is in a list of options"""
    if value in options:
        return True
    else:
        return False
        # logging.error('\nERROR: {} must be: {}\n'.format(
        #     value, ', '.join(options)))
        # raise ValueError('\nERROR: {} must be: {}\n'.format(
        #     value, ', '.join(options)))


# def millis(input_dt):
#     """Convert datetime to milliseconds since epoch"""
#     # Python 3 (or 2 with future module)
#     return 1000 * int(calendar.timegm(input_dt.timetuple()))
#     # Python 2
#     # return 1000 * long(calendar.timegm(input_dt.timetuple()))
#     # return 1000 * long(time.mktime(input_dt.timetuple()))


def parse_int_set(nputstr=""):
    """Return list of numbers given a string of ranges

    http://thoughtsbyclayg.blogspot.com/2008/10/parsing-list-of-numbers-in-python.html
    """
    selection = set()
    invalid = set()
    # tokens are comma seperated values
    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        try:
            # typically tokens are plain old integers
            selection.add(int(i))
        except:
            # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items seperated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token) - 1]
                    for x in range(first, last + 1):
                        selection.add(x)
            except:
                # not an int and not a range...
                invalid.add(i)
    # Report invalid tokens before returning valid selection
    # print "Invalid set: " + str(invalid)
    return selection