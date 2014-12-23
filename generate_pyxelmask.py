#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generates orbital quality mask for SCIA channel 8. 
Individual criteria are floating point [0.0 .. 1.0]. 
These are combined and are thresholded, and result in the boolean `combinedFlag'.
"""

from __future__ import division, print_function

import numpy as np
from pixelquality import PixelQuality 
import sciamachy

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import logging
    from datetime import datetime
    from os.path import basename, isfile

    from envisat import parseOrbitList

    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output', dest='output_fname', type=str,
                        help="output file name")
    parser.add_argument('--log', dest='loglevel', type=str, 
                        choices=("DEBUG", "INFO", "WARN", "ERROR", "FATAL"), 
                        default="INFO", help="logging level")
    parser.add_argument('-c', action='store_true', dest="sdmf30_compat", 
                        help="compatibility, approximate SDMF3.0 pixel mask")
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', 
                        version='%(prog)s 0.1')
    parser.add_argument('-P', '--path', dest='path')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--orbitrange', default='all', 
                        help='sets orbit range f.e. "43000-44000", "all"', type=parseOrbitList)
    parser.add_argument('-p', '--pixnr', action='store', type=int, default=None,
                        dest='pixnr', help="pixel number to be examined [0..1023]")
    args = parser.parse_args()

    if args.orbitrange is not None:
        orbit_range = args.orbitrange
    else:
        orbit_range = 4151,53200
    if args.path is not None:
        path = args.path
    else:
        path = "."
    if args.output_fname is not None:
        output_fname = args.output_fname

    print("sdmf30_compat =", args.sdmf30_compat)
    print("path =", path)
    print("output_fname =", output_fname)
    print("orbit_range =", orbit_range)
    print("log level =", args.loglevel)

    #
    # setup log level
    #

    numeric_level = logging.getLevelName(args.loglevel.upper())

    #
    # create log with name with form generate_pyxelmask_YY-MMM-DD_N.log 
    #

    timestamp = datetime.now().strftime("%Y-%m-%d")
    i = 0
    while True:
        postfix = "_"+str(i) if i>0 else ""
        logname = basename(__file__)+"_"+timestamp+postfix+".log"
        if not isfile(logname):
            break
        i += 1

    #
    # open log
    #

    logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()

    fileHandler = logging.FileHandler(logname)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logger = rootLogger 
    rootLogger.setLevel(numeric_level)

    #
    # log parameters
    #

    logger.info("sdmf30_compat ="+str(args.sdmf30_compat))
    logger.info("path ="+str(path))
    logger.info("output_fname ="+str(output_fname))
    logger.info("orbit_range ="+str(orbit_range))
    logger.info("log level ="+str(args.loglevel))

    #
    # generate the masks
    #

    p = PixelQuality(sdmf30_compat=args.sdmf30_compat)

    for orbit in range(orbit_range[0],orbit_range[1]):
        try:
            p.calculate(orbit)
        except sciamachy.SdmfDataMissingError as e:
            logger.exception(e)
            logger.warning("calculation failed for orbit %d!" % orbit)
            continue

        logger.info("orbit %d computed" % orbit)
        p.write(directory=path, fname=output_fname)
        logger.info("Pixel mask data for orbit %d written to db." % orbit)
