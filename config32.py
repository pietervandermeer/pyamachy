from __future__ import print_function, division

"""
Handles configuration of SDMF 3.2 parameters.
"""

import ConfigParser

def get_floats(tags, parser, part):
    """
    parse floating point variables from config and put them in a dict under the desired tags 
    """
    d = {}
    for tag in tags:
        d[tag] = float(parser.get(part, tag))
    return d

def get_floatlists(tags, parser, part):
    """
    parse comma-separated lists of floats from config and put them in a dict under the desired tags 
    """
    d = {}
    for tag in tags:
        string = parser.get(part, tag)
        d[tag] = [float(s) for s in string.split(',')]
    return d

def load(config_file):
    """ 
    load configuration from file 

    Parameters
    ----------

    config_file : string
        name of the configuration file

    Returns
    -------
    get_config : dict
        contains all relevant configuration settings as readily-usable types (int, boolean, float, arrays instead of strings)
    """
    import string

    parser=ConfigParser.SafeConfigParser()
    parser.readfp(config_file)
    d = {}

    #
    # parse file names
    #

    d['db_dir'] = parser.get('Global','masterdirectory')
    d['extract_fname'] = parser.get('Global','extract_file')
    d['dark_long_fname'] = parser.get('Global','dark_long_file')
    d['dark_short_fname'] = parser.get('Global','dark_short_file')
    d['pixelmask_fname'] = parser.get('Global','pixelmask_file')
    d['statedarkch6p_fname'] = parser.get('Global','statedarkch6p_file')
    d['statedarkch8_fname'] = parser.get('Global','statedarkch8_file')
    d['darklimbch6_fname'] = parser.get('Global','darklimbch6_file')
    d['darklimbch8_fname'] = parser.get('Global','darklimbch8_file')
    d['darklimbch8_fname'] = parser.get('Global','darklimbch8_file')

    #
    # parse lists of floats.. useful for multi-channel configurations (although this class is single-channel.. it's left in for future)
    #

    floatlist_tags = ['dc_sat_time','dc_err_thres','res_thres','max_signal','deadflaggingstates']
    proc_floatlists = get_floatlists(floatlist_tags, parser, 'Processor')

    #
    # parse floats (mostly thresholds and weights) used for channel 8 processing
    #

    float_tags = ['w_err','w_res','w_noise','w_sun','w_wls','thresh_sun','thresh_wls','max_noise']
    proc_floats = get_floats(float_tags, parser, 'Processor')
    return dict(d.items() + proc_floatlists.items() + proc_floats.items())
